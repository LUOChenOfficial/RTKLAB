/*------------------------------------------------------------------------------
* rtk_gls.c : current-epoch GLS solution separation
*
* notes  : Float protective subsets share the full-state prior.
*          Gaussian bounds exclude prior bias and RANSAC selection calibration.
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "rtk_int_internal.h"
#include <limits.h>

static int excluded(int sat, const unsigned char *mask)
{
    return sat>0&&sat<=MAXSAT&&mask[sat-1];
}

/* A DD row is ref-target. If ref is excluded, subtract ref-anchor to obtain
 * anchor-target. The same sparse transform must be applied to v, H and Q. */
/* construct satellite-exclusion DD transform --------------------------------
* args   : int *flags      I   DD row flags (nv x 1)
*          int nv          I   number of input DD rows
*          int sat1,sat2   I   excluded satellites (0: none)
*          int *rows       O   retained row indices
*          int *ixref      O   subtracted reference row indices (-1: none)
* return : number of output rows
* notes  : Apply the same transform to residuals, design and covariance.
*-----------------------------------------------------------------------------*/
EXPORT int rtkglsddmask(const int *flags, int nv, const unsigned char *mask,
                     int *rows, int *ixref)
{
    int i,j,n=0;
    for (i=0;i<nv;i++) {
        int ref=(flags[i]>>16)&255,target=(flags[i]>>8)&255,
            anchor= -1;
        if (ref<=0||target<=0||excluded(target,mask)) continue;
        if (excluded(ref,mask)) {
            for (j=0;j<nv;j++) {
                int r=(flags[j]>>16)&255,t=(flags[j]>>8)&255;
                if (r==ref&&(flags[j]&255)==(flags[i]&255)&&
                    t>0&&!excluded(t,mask)) {
                    anchor=j;
                    break;
                }
            }
            if (anchor<0||anchor==i) continue;
        }
        rows[n]=i;
        ixref[n++]=anchor;
    }
    return n;
}

EXPORT int rtkglsddmat(const int *flags, int nv, int sat1, int sat2,
                       int *rows, int *ixref)
{
    unsigned char mask[MAXSAT]={0};
    if (sat1>0&&sat1<=MAXSAT) mask[sat1-1]=1;
    if (sat2>0&&sat2<=MAXSAT) mask[sat2-1]=1;
    return rtkglsddmask(flags,nv,mask,rows,ixref);
}

static int nextcomb(int *ix, int order, int ns)
{
    int i,j;
    for (i=order-1;i>=0;i--) if (ix[i]<ns-order+i) {
        ix[i]++;
        for (j=i+1;j<order;j++) ix[j]=ix[j-1]+1;
        return 1;
    }
    return 0;
}

static int psd3(const double *q, double tol)
{
    double l[9]={0};
    int i,j,k;

    for (i=0;i<3;i++)
        for (j=0;j<=i;j++) {
            /* Test Q + tol I. A rank-one separation covariance can have a
         * tiny pivot and a larger off-diagonal; skipping that pivot falsely
         * rejects valid PSD matrices. */
            double s=0.5*(q[i+j*3]+q[j+i*3])+(i==j?tol:0.0);
            if (!isfinite(s)) return 0;
            for (k=0;k<j;k++) s-=l[i+k*3]*l[j+k*3];
            if (i==j) {
                if (s<=0.0) return 0;
                l[i+j*3]=sqrt(s);
            } else l[i+j*3]=s/l[j+j*3];
        }
    return 1;
}

/* Exact position block of the full-state KF update, without allocating a
 * second state/filter: F = first 3 rows of P H; Q = H' P H + R. */
/* projected full-state KF position update -----------------------------------
* args   : double *Q       I   innovation covariance (nv x nv, m^2)
*          double *F       I   position block of P*H (3 x nv)
*          double *v       I   innovations (nv x 1, m)
*          double *x0      I   prior ECEF position (3 x 1, m)
*          double *P0      I   prior ECEF position covariance (3 x 3, m^2)
*          int nv,nr       I   input and transformed observation counts
*          int *rows       I   retained row indices (nr x 1)
*          int *ixref      I   subtracted row indices (nr x 1, -1: none)
*          double *pos     O   updated ECEF position (3 x 1, m)
*          double *cov     O   updated covariance (3 x 3, m^2)
* return : status (1:ok,0:error)
* notes  : All matrices are column-major. Output is valid only on success.
*-----------------------------------------------------------------------------*/
EXPORT int rtkglsupdate(const double *Q, const double *F, const double *v,
                        const double *x0, const double *P0, int nv,
                        const int *rows, const int *ixref, int nr, double *pos,
                        double *cov)
{
    double *qs=NULL,*fs=NULL,*vs=NULL,*scale=NULL;
    int i,j,k,ok=0;

    if (nv<=0||nr<=0||nr>nv) return 0;
    qs=mat(nr,nr);
    fs=mat(3,nr);
    vs=mat(nr,1);
    scale=mat(nr,1);
    if (!qs||!fs||!vs||!scale) goto done;
    for (i=0;i<nr;i++) {
        int r=rows[i],a=ixref[i];
        if (r<0||r>=nv||a>=nv) goto done;
        vs[i]=v[r]-(a>=0?v[a]:0.0);
        for (k=0;k<3;k++) fs[k+i*3]=F[k+r*3]-(a>=0?F[k+a*3]:0.0);
        for (j=0;j<nr;j++) {
            int s=rows[j],b=ixref[j];
            if (s<0||s>=nv||b>=nv) goto done;
            qs[i+j*nr]=Q[r+s*nv]-(a>=0?Q[a+s*nv]:0.0)-
                             (b>=0?Q[r+b*nv]:0.0)+
                             (a>=0&&b>=0?Q[a+b*nv]:0.0);
        }
    }
    /* Scale and factor the innovation covariance. Explicit inversion loses
     * enough precision in large ambiguity priors to make rank-one separation
     * covariances appear indefinite. */
    for (i=0;i<nr;i++) {
        if (!(qs[i+i*nr]>0.0)||!isfinite(qs[i+i*nr])) goto done;
        scale[i]=sqrt(qs[i+i*nr]);
    }
    for (i=0;i<nr;i++)
        for (j=0;j<=i;j++) {
            double s=0.5*(qs[i+j*nr]+qs[j+i*nr])/
                       (scale[i]*scale[j]);
            for (k=0;k<j;k++) s-=qs[i+k*nr]*qs[j+k*nr];
            if (i==j) {
                if (!(s>0.0)) goto done;
                qs[i+j*nr]=sqrt(s);
            } else qs[i+j*nr]=s/qs[j+j*nr];
        }
    for (i=0;i<nr;i++) {
        double s=vs[i]/scale[i];
        for (j=0;j<i;j++) s-=qs[i+j*nr]*vs[j];
        vs[i]=s/qs[i+i*nr];
        for (k=0;k<3;k++) {
            s=fs[k+i*3]/scale[i];
            for (j=0;j<i;j++) s-=qs[i+j*nr]*fs[k+j*3];
            fs[k+i*3]=s/qs[i+i*nr];
        }
    }
    matcpy(pos,x0,3,1);
    matcpy(cov,P0,3,3);
    matmul("NN",3,1,nr,1.0,fs,vs,1.0,pos);
    matmul("NT",3,3,nr,-1.0,fs,fs,1.0,cov);
    for (i=0;i<3;i++) {
        if (!isfinite(pos[i])) goto done;
        for (j=0;j<i;j++)
            cov[i+j*3]=cov[j+i*3]=
                0.5*(cov[i+j*3]+cov[j+i*3]);
    }
    if (!psd3(cov,1e-7)) goto done;
    ok=1;
done:
    free(qs);
    free(fs);
    free(vs);
    free(scale);
    return ok;
}

/* Minimized augmented-GLS objective (real residual + prior correction).
 * With a full-rank prior for every active state, redundancy is nv, not
 * nv-3 nor nv-number_of_states. Factor S rather than explicitly invert it. */
EXPORT int rtkglsvarscale(const double *S, const double *v, int nv,
                          double *vtpv)
{
    double *L=NULL,*z=NULL,*s=NULL,total=0.0;
    int i,j,k,ok=0;
    if (!S||!v||!vtpv||nv<=0) return 0;
    L=zeros(nv,nv);z=mat(nv,1);s=mat(nv,1);
    if (!L||!z||!s) goto done;
    for (i=0;i<nv;i++) {
        if (!(S[i+i*nv]>0.0)||!isfinite(S[i+i*nv])||!isfinite(v[i])) goto done;
        s[i]=sqrt(S[i+i*nv]);
    }
    for (i=0;i<nv;i++) {
        for (j=0;j<=i;j++) {
            double a=0.5*(S[i+j*nv]+S[j+i*nv])/(s[i]*s[j]);
            if (!isfinite(a)) goto done;
            for (k=0;k<j;k++) a-=L[i+k*nv]*L[j+k*nv];
            if (i==j) {
                if (!(a>0.0)) goto done;
                L[i+j*nv]=sqrt(a);
            } else L[i+j*nv]=a/L[j+j*nv];
        }
        z[i]=v[i]/s[i];
        for (j=0;j<i;j++) z[i]-=L[i+j*nv]*z[j];
        z[i]/=L[i+i*nv];total+=z[i]*z[i];
    }
    if (!isfinite(total)) goto done;
    *vtpv=total;ok=1;
done:
    free(L);free(z);free(s);return ok;
}

EXPORT void rtkimglsprior(rtk_t *rtk, const double *x, const double *P,
                          const double *H, const double *R, int nx, int nv)
{
    rtkim_t *mon=&rtk->intg;
    double *pc=NULL,*hc=NULL,*f=NULL;
    int *ix=NULL,i,j,k,n=0;
    mon->glsready=0;
    rtk->glsassessed=0;
    rtk->sol.post_varfactor=1.0;
    rtk->sol.post_vtpv=0.0;rtk->sol.post_dof=0;
    if (!rtk->opt.intmon||rtk->intchild||rtk->opt.niter!=1||
        rtk->opt.mode==PMODE_MOVEB||rtk->opt.modear==ARMODE_FIXHOLD||
        rtk->opt.glomodear==GLO_ARMODE_FIXHOLD||nv<=0||
        nv>RTKIM_MAX_BIAS_ROWS||mon->bnv!=nv||mon->bnx!=nx)
        return;
    ix=imat(nx,1);
    if (!ix) return;
    for (i=0;i<nx;i++)
        if (x[i]!=0.0&&P[i+i*nx]>0.0) ix[n++]=i;
    if (n<3||ix[0]!=0||ix[1]!=1||ix[2]!=2) goto done;
    pc=mat(n,n);
    hc=mat(n,nv);
    f=mat(n,nv);
    if (!pc||!hc||!f) goto done;
    free(mon->glsQ);
    free(mon->glsF);
    mon->glsQ=mat(nv,nv);
    mon->glsF=mat(3,nv);
    if (!mon->glsQ||!mon->glsF) goto done;
    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) pc[i+j*n]=P[ix[i]+ix[j]*nx];
        for (j=0;j<nv;j++) hc[i+j*n]=H[ix[i]+j*nx];
    }
    matmul("NN",n,nv,n,1.0,pc,hc,0.0,f);
    matcpy(mon->glsQ,R,nv,nv);
    matmul("TN",nv,nv,n,1.0,hc,f,1.0,mon->glsQ);
#if RTK_GLS_POST_SCALE
    if (!rtkglsvarscale(mon->glsQ,mon->bv,nv,&rtk->sol.post_vtpv))
        goto done;
    rtk->sol.post_dof=nv;
    /* Tiny numerical floor only; do not force sigma0 >= 1. */
    rtk->sol.post_varfactor=fmax(1e-12,rtk->sol.post_vtpv/nv);
#endif
    for (i=0;i<3;i++) {
        mon->glspos[i]=x[i];
        for (j=0;j<3;j++) mon->glscov[i+j*3]=P[i+j*nx];
        for (k=0;k<nv;k++) mon->glsF[i+k*3]=f[i+k*n];
    }
    mon->glsready=1;
done:
    free(ix);
    free(pc);
    free(hc);
    free(f);
}

static double normtail(double p)
{
    double lo=0.0,hi=12.0;
    int i;

    for (i=0;i<70;i++) {
        double mid=0.5*(lo+hi);
        if (0.5*erfc(mid/sqrt(2.0))>p) lo=mid;
        else hi=mid;
    }
    return 0.5*(lo+hi);
}

static void covlocal(const double *E, const double *cov, double *q)
{
    double tmp[9];

    matmul("NN",3,3,3,1.0,E,cov,0.0,tmp);
    matmul("NT",3,3,3,1.0,tmp,E,0.0,q);
}

/* Require a full-rank current code geometry, not merely a finite prior. */
static int chkgeom(const rtkim_t *mon, const int *rows, const int *ixref,
                   int nr)
{
    double normal[9]={0};
    int i,j,k,n=0;
    for (i=0;i<nr;i++) {
        double h[3];
        int r=rows[i],a=ixref[i];
        if (((mon->bflg[r]>>4)&15)!=1) continue;
        n++;
        for (j=0;j<3;j++)
            h[j]=mon->bH[j+r*mon->bnx]-
                   (a>=0?mon->bH[j+a*mon->bnx]:0.0);
        for (j=0;j<3;j++)
            for (k=0;k<3;k++) normal[j+k*3]+=h[j]*h[k];
    }
    return n>=3&&matinv(normal,3)==0;
}

static int boundsub(rtk_t *rtk, const double *pos, const double *cov,
                     const double *outcov, const double *E, double kfa,
                     double kn, int mode, int sat1, int sat2,
                     int id)
{
    rtkim_t *mon=&rtk->intg;
    rtkim_pl_t *pl=&mon->pl;
    double q[9],delta[9],qss[9],d[3],sep[3],sigma[3],sss[3],axis[3],h;
    int i;
    for (i=0;i<9;i++) delta[i]=cov[i]-outcov[i];
    if (!psd3(delta,1e-7)) return 0;
    covlocal(E,cov,q);
    covlocal(E,delta,qss);
    for (i=0;i<3;i++) d[i]=rtk->sol.rr[i]-pos[i];
    matmul("NN",3,1,3,1.0,E,d,0.0,sep);
    for (i=0;i<3;i++) {
        sigma[i]=sqrt(fmax(q[i+i*3],0.0));
        sss[i]=sqrt(fmax(qss[i+i*3],0.0));
    }
    if (mode!=RTKIM_F_NONE) {
        rtkimglsout(rtk->sol.time,mode,sat1,sat2,sep,sigma);
        /* Observation-only and joint modes reuse a conservative float bound,
         * but must have separate, disjoint probability masses offline. */
        if (mode==4||mode==5||mode==7)
            rtkimglsout(rtk->sol.time,mode==4?1:mode==5?3:6,sat1,sat2,sep,
                        sigma);
    }
    /* Realized separation is a detection statistic, not a PL term. */
    for (i=0;i<3;i++) axis[i]=kfa*sss[i]+kn*sigma[i];
    h=hypot(axis[0],axis[1]);
    if (h>pl->hpl) {
        pl->hpl=h;
        pl->hsrc=id;
        pl->hmode=mode;
        pl->hsat=sat1;
        mon->glshsat2=sat2;
        for (i=0;i<3;i++) {
            pl->hsig[i]=sqrt(fmax(q[i+i*3],0.0));
            pl->hsss[i]=sss[i];
            pl->hsep[i]=sep[i];
        }
    }
    if (axis[2]>pl->vpl) {
        pl->vpl=axis[2];
        pl->vsrc=id;
        pl->vmode=mode;
        pl->vsat=sat1;
        mon->glsvsat2=sat2;
        for (i=0;i<3;i++) {
            pl->vsig[i]=sqrt(fmax(q[i+i*3],0.0));
            pl->vsss[i]=sss[i];
            pl->vsep[i]=sep[i];
        }
    }
    return 1;
}

EXPORT void rtkimglsassess(rtk_t *rtk)
{
    rtkim_t *mon=&rtk->intg;
    rtkim_pl_t *pl=&mon->pl;
    int sats[MAXSAT],ns=0,nv=mon->bnv,*rows=NULL,*ixref=NULL;
    int i,j,k,order,more,ix[8],maxorder,id=0,nr;
    int fixed=rtk->sol.stat==SOLQ_FIX,obs_modes=0;
    double combinations=1.0;
    unsigned char mask[MAXSAT];
    double mainpos[3],maincov[9],outcov[9],pos[3],cov[9],E[9],llh[3],qmain[9],kn,kfa;
    double pfa=rtk->opt.intpfa,pmd=rtk->opt.intpmd;
    double fdescore=0.0;
#if defined(RTK_GLS_DIAG_TESTS) && RTK_GLS_DIAG_TESTS
    double diag_float=0.0,diag_joint=0.0,diag_fix_sep[3]={0},diag_fix_std[3]={0};
#endif
    double scale=rtk->sol.post_varfactor>0.0?rtk->sol.post_varfactor:1.0;
    int fdeorder=0;
    clock_t begin=clock();
    memset(pl,0,sizeof(*pl));
    memset(&mon->fde,0,sizeof(mon->fde));
    memset(&mon->act,0,sizeof(mon->act));
    rtk->intfdemode=rtk->intfdesat=0;
    rtk->glsassessed=1;
    memset(rtk->glsfdemask,0,sizeof(rtk->glsfdemask));
    mon->glsexp=mon->glsvalid=mon->glsmiss=mon->glsalert=0;
    mon->glshsat2=mon->glsvsat2=0;
    mon->glsmatch=INFINITY;
    mon->glsgeofail=mon->glssolfail=mon->glscovfail=0;
    mon->glsobsalert=mon->glsfixalert=0;
    mon->glsobsratio=mon->glsfixratio=mon->glsfloatratio=0.0;
    mon->glsncandidate=0;
    if (!mon->glsready||rtk->sol.stat==SOLQ_NONE||
        (fixed&&!rtk->opt.intamb)||!rtk->opt.intsat)
        goto unavailable;
    rows=imat(nv,1);
    ixref=imat(nv,1);
    if (!rows||!ixref) goto unavailable;
    for (i=0;i<nv;i++)
        for (j=0;j<2;j++) {
            int sat=(mon->bflg[i]>>(j?8:16))&255;
            if (sat<=0||sat>MAXSAT) goto unavailable;
            for (k=0;k<ns;k++)
                if (sats[k]==sat) break;
            if (k==ns) sats[ns++]=sat;
        }
    maxorder=ns<RTK_GLS_MAX_FAULTS?ns:RTK_GLS_MAX_FAULTS;
    for (order=1;order<=maxorder;order++) {
        combinations=combinations*(ns-order+1)/order;
        if (combinations>INT_MAX/2-1-obs_modes) goto unavailable;
        obs_modes+=(int)(combinations+0.5);
    }
    /* Fixed outputs also test joint fixed-to-subset separations. */
    mon->glsexp=obs_modes*(fixed?2:1)+(fixed?1:0);
    if (RTK_GLS_FDE_SEARCH&&rtk->opt.intfde&&rtk->intreproc<=0) {
        int capacity=ns+(maxorder>=2?ns*(ns-1)/2:0);
        if (capacity>mon->glscapacity) {
            rtkim_glscand_t *pool=(rtkim_glscand_t *)realloc(mon->glscandidates,
                                          capacity*sizeof(rtkim_glscand_t));
            if (pool) {mon->glscandidates=pool;mon->glscapacity=capacity;}
        }
    }
    if (!(pfa>0.0&&pfa<1.0)) pfa=RTK_GLS_PFA;
    if (!(pmd>0.0&&pmd<1.0)) pmd=RTK_GLS_PMD;
    kn=normtail(pmd/(6.0*(mon->glsexp+1)));
    kfa=normtail(pfa/(6.0*(mon->glsexp+1)));
    for (i=0;i<nv;i++) {
        rows[i]=i;
        ixref[i]= -1;
    }
    if (!chkgeom(mon,rows,ixref,nv)||
        !rtkglsupdate(mon->glsQ,mon->glsF,mon->bv,mon->glspos,mon->glscov,
                      nv,rows,ixref,nv,mainpos,maincov))
        goto unavailable;
    mon->glsmatch=0.0;
    for (i=0;i<3;i++)
        mon->glsmatch=fmax(mon->glsmatch,fabs(mainpos[i]-rtk->x[i]));
    if (mon->glsmatch>1e-4) goto unavailable;
    for (i=0;i<9;i++) maincov[i]*=scale;
    ecef2pos(rtk->sol.rr,llh);
    xyz2enu(llh,E);
    covlocal(E,maincov,qmain);
    for (i=0;i<3;i++) pl->msig[i]=sqrt(fmax(qmain[i+i*3],0.0));
    for (i=0;i<3;i++) for (j=0;j<3;j++)
        outcov[i+j*3]=fixed?scale*rtk->Pa[i+j*rtk->na]:maincov[i+j*3];
    {
        double qnom[9],nomcov[9],zero[3]={0.0},sigma[3];
        for (i=0;i<3;i++)
            for (j=0;j<3;j++)
                nomcov[i+j*3]=
                    fixed?scale*rtk->Pa[i+j*rtk->na]:maincov[i+j*3];
        covlocal(E,nomcov,qnom);
        for (i=0;i<3;i++) sigma[i]=sqrt(fmax(qnom[i+i*3],0.0));
        /* Formal correct-fix covariance only: selection calibration is still
         * required before interpreting this as a conditional error bound. */
        rtkimglsout(rtk->sol.time,0,0,0,zero,sigma);
    }
    if (!boundsub(rtk,mainpos,maincov,outcov,E,kfa,kn,
                  fixed?RTKIM_F_NFIX:RTKIM_F_NONE,0,0,++id))
        goto unavailable;
    if (fixed) {
        double qfix[9],delta[9],qe[9],diff[3],sep[3];
        mon->glsvalid++;
        for (i=0;i<3;i++)
            for (j=0;j<3;j++) {
                qfix[i+j*3]=scale*rtk->Pa[i+j*rtk->na];
                delta[i+j*3]=maincov[i+j*3]-qfix[i+j*3];
            }
        if (!psd3(delta,1e-7)) goto unavailable;
        covlocal(E,delta,qe);
        for (i=0;i<3;i++) diff[i]=rtk->sol.rr[i]-mainpos[i];
        matmul("NN",3,1,3,1.0,E,diff,0.0,sep);
        for (i=0;i<3;i++) {
            double ratio=fabs(sep[i])/(kfa*sqrt(fmax(qe[i+i*3],0.0))+1e-6);
            mon->glsfixratio=fmax(mon->glsfixratio,ratio);
#if defined(RTK_GLS_DIAG_TESTS) && RTK_GLS_DIAG_TESTS
            diag_fix_sep[i]=sep[i];diag_fix_std[i]=sqrt(fmax(qe[i+i*3],0.0));
#endif
            if (ratio>1.0) mon->glsalert=mon->glsfixalert=1;
        }
    }
    for (order=1;order<=maxorder;order++) {
        for (i=0;i<order;i++) ix[i]=i;
        for (more=1;more;more=nextcomb(ix,order,ns)) {
            double delta[9],qe[9],diff[3],sep[3];
            int sat1=sats[ix[0]],sat2=order>1?sats[ix[1]]:0;
            int mode=order>2?6:order==2?3:RTKIM_F_SAT;
            memset(mask,0,sizeof(mask));
            for (i=0;i<order;i++) mask[sats[ix[i]]-1]=1;
            nr=rtkglsddmask(mon->bflg,nv,mask,rows,ixref);
            if (!chkgeom(mon,rows,ixref,nr)) {
                mon->glsgeofail++;
                mon->glsmiss+=fixed?2:1;
                continue;
            }
            if (!rtkglsupdate(mon->glsQ,mon->glsF,mon->bv,mon->glspos,
                              mon->glscov,nv,rows,ixref,nr,pos,cov)) {
                mon->glssolfail++;
                mon->glsmiss+=fixed?2:1;
                continue;
            }
            for (i=0;i<9;i++) {cov[i]*=scale;delta[i]=cov[i]-maincov[i];}
            if (!psd3(delta,1e-7)) {
                mon->glscovfail++;
                mon->glsmiss+=fixed?2:1;
                continue;
            }
            mon->glsvalid+=fixed?2:1;
            covlocal(E,delta,qe);
            for (i=0;i<3;i++) diff[i]=mainpos[i]-pos[i];
            matmul("NN",3,1,3,1.0,E,diff,0.0,sep);
            {
            double score=0.0;
            for (i=0;i<3;i++) {
                double ratio=fabs(sep[i])/
                               (kfa*sqrt(fmax(qe[i+i*3],0.0))+1e-6);
                mon->glsobsratio=fmax(mon->glsobsratio,ratio);
                mon->glsfloatratio=fmax(mon->glsfloatratio,ratio);
#if defined(RTK_GLS_DIAG_TESTS) && RTK_GLS_DIAG_TESTS
                diag_float=fmax(diag_float,ratio);
#endif
                if (ratio>1.0) mon->glsalert=mon->glsobsalert=1;
                score=fmax(score,ratio);
            }
            if (RTK_GLS_FDE_SEARCH&&rtk->opt.intfde&&rtk->intreproc<=0&&order<=2&&
                mon->glsncandidate<mon->glscapacity) {
                rtkim_glscand_t *candidate=mon->glscandidates+mon->glsncandidate++;
                candidate->order=order;candidate->sat[0]=sat1;candidate->sat[1]=sat2;
                candidate->score=score;
            }
            if (fixed) {
                /* Match the joint PL threshold to its actual test: protective
                 * float subset minus output fixed solution. The covariance
                 * difference is conditional on the assumed integer candidate. */
                for (i=0;i<9;i++) delta[i]=cov[i]-outcov[i];
                if (!psd3(delta,1e-7)) {
                    mon->glscovfail++;
                    mon->glsvalid-=2;mon->glsmiss+=2;
                    continue;
                }
                covlocal(E,delta,qe);
                for (i=0;i<3;i++) diff[i]=rtk->sol.rr[i]-pos[i];
                matmul("NN",3,1,3,1.0,E,diff,0.0,sep);
                for (i=0;i<3;i++) {
                    double ratio=fabs(sep[i])/
                        (kfa*sqrt(fmax(qe[i+i*3],0.0))+1e-6);
                    if (ratio>1.0) mon->glsalert=mon->glsobsalert=1;
                    mon->glsobsratio=fmax(mon->glsobsratio,ratio);
#if defined(RTK_GLS_DIAG_TESTS) && RTK_GLS_DIAG_TESTS
                    diag_joint=fmax(diag_joint,ratio);
#endif
                    score=fmax(score,ratio);
                }
            }
            /* Experimental candidate: smallest alarmed exclusion order, then
             * largest separation score. This is not proven fault isolation.
             * A full retry must pass the monitor before finite PL is accepted. */
            if (score>1.0&&(!fdeorder||order<fdeorder||
                           (order==fdeorder&&score>fdescore))) {
                fdeorder=order; fdescore=score;
                memcpy(rtk->glsfdemask,mask,sizeof(mask));
            }
            }
            /* Detection threshold + protective float-subset noise. */
            if (!boundsub(rtk,pos,cov,outcov,E,kfa,kn,
                          fixed?(order>2?7:order==2?5:4):mode,sat1,sat2,++id)) {
                mon->glscovfail++;
                mon->glsvalid-=fixed?2:1;mon->glsmiss+=fixed?2:1;
                continue;
            }
            if (order>2&&(pl->hsrc==id||pl->vsrc==id)) {
                rtkimdbg("GLS_HIGH_SOURCE %s %d %d %d",time_str(rtk->sol.time,3),
                         id,pl->hsrc==id,pl->vsrc==id);
                for (i=0;i<order;i++) rtkimdbg(" %d",sats[ix[i]]);
                rtkimdbg("\n");
            }
        }
    }
    pl->hpl0=pl->hpl;
#if defined(RTK_GLS_DIAG_TESTS) && RTK_GLS_DIAG_TESTS
    rtkimdbg("GLS_TEST_SPLIT %s retry=%d float=%.9g joint=%.9g fix=%.9g kfa=%.9g sigma0=%.9g sep=%.9g,%.9g,%.9g std=%.9g,%.9g,%.9g\n",
        time_str(rtk->sol.time,3),rtk->intreproc,diag_float,diag_joint,mon->glsfixratio,kfa,sqrt(scale),
        diag_fix_sep[0],diag_fix_sep[1],diag_fix_sep[2],diag_fix_std[0],diag_fix_std[1],diag_fix_std[2]);
#endif
    pl->vpl0=pl->vpl;
    if (rtk->opt.intfde&&rtk->intreproc<=0&&mon->glsvalid==mon->glsexp) {
        if (RTK_GLS_FDE_FLOAT_FIRST&&fixed&&mon->glsfloatratio<=1.0&&
            (mon->glsobsalert||mon->glsfixalert)) {
            /* Float consistency passing does not prove a wrong integer, but
             * permits testing float recovery before attributing it to a sat. */
            mon->fde.act=RTKIM_A_FLOAT;
            mon->fde.mode=RTKIM_F_NFIX;
            mon->fde.score=fmax(mon->glsobsratio,mon->glsfixratio);
            memset(rtk->glsfdemask,0,sizeof(rtk->glsfdemask));
        } else if (mon->glsobsalert&&fdeorder) {
            mon->fde.act=RTKIM_A_SAT;
            mon->fde.mode=RTKIM_F_SAT;
            mon->fde.score=fdescore;
            for (i=0;i<MAXSAT;i++) if (rtk->glsfdemask[i]) {
                if (!mon->fde.sat) mon->fde.sat=i+1;
                rtkimdbg("GLS_FDE_EXCLUDE %s %d\n",time_str(rtk->sol.time,3),i+1);
            }
        } else if (!mon->glsobsalert&&mon->glsfixalert) {
            mon->fde.act=RTKIM_A_FLOAT;
            mon->fde.mode=RTKIM_F_NFIX;
            mon->fde.score=mon->glsfixratio;
        }
        mon->act=mon->fde;
        rtk->intfdemode=mon->fde.mode;
        rtk->intfdesat=mon->fde.sat;
    }
    if (mon->glsvalid!=mon->glsexp||mon->glsalert) goto unavailable;
    goto done;
unavailable:
    pl->hpl=pl->vpl=INFINITY;
done:
    mon->glstime=(clock()-begin)/(double)CLOCKS_PER_SEC;
    rtk->inthpl=pl->hpl;
    rtk->intvpl=pl->vpl;
    free(rows);
    free(ixref);
}
