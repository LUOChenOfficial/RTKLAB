/*------------------------------------------------------------------------------
* rtk_slope.c : innovation slope integrity monitoring
*
* notes  : Thresholds, KF prior projection and position protection levels.
*-----------------------------------------------------------------------------*/
#include "rtk_int_internal.h"

static double gammap(double a, double x)
{
    double sum,term,b,c,d,h,an,del,fac;
    int i;

    if (x<=0.0) return 0.0;
    fac=exp(a*log(x)-x-lgamma(a));
    if (x<a+1.0) {
        sum=term=1.0/a;
        for (i=1;i<10000;i++) {
            term*=x/(a+i);
            sum+=term;
            if (fabs(term)<fabs(sum) *1E-14) break;
        }
        return fmin(1.0,fmax(0.0,sum*fac));
    }
    b=x+1.0-a;
    c=1E300;
    d=1.0/b;
    h=d;
    for (i=1;i<10000;i++) {
        an= -i*(i-a);
        b+=2.0;
        d=an*d+b;
        if (fabs(d)<1E-300) d=1E-300;
        c=b+an/c;
        if (fabs(c)<1E-300) c=1E-300;
        d=1.0/d;
        del=d*c;
        h*=del;
        if (fabs(del-1.0)<1E-14) break;
    }
    return fmin(1.0,fmax(0.0,1.0-fac*h));
}

static double ncchicdf(double x, int df, double lam)
{
    double mu=lam*0.5,w0,w,sum;
    int k,k0;

    if (mu<=0.0) return gammap(df*0.5,x*0.5);
    k0=(int)floor(mu);
    w0=exp(-mu+k0*log(mu)-lgamma(k0+1.0));
    sum=w0*gammap(df*0.5+k0,x*0.5);
    w=w0;
    for (k=k0;k>0;k--) {
        w*=k/mu;
        sum+=w*gammap(df*0.5+k-1,x*0.5);
        if (w<1E-15) break;
    }
    w=w0;
    for (k=k0+1;k<k0+10000;k++) {
        w*=mu/k;
        sum+=w*gammap(df*0.5+k,x*0.5);
        if (w<1E-15) break;
    }
    return fmin(1.0,fmax(0.0,sum));
}

/* slope detector thresholds -------------------------------------------------
* args   : int dof         I   innovation degrees of freedom
*          double pfa,pmd  I   false-alarm and missed-detection probabilities
*          double *thres  O   central chi-square threshold
*          double *lam    O   noncentrality parameter
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
EXPORT int rtkslopethres(int dof, double pfa, double pmd, double *thres,
                         double *lam)
{
    double lo=0.0,hi,mid,tau;
    int i;
    if (!thres||!lam||dof<=0||dof>RTKIM_MAX_BIAS_ROWS||
        pfa<1E-12||pfa>=1.0||pmd<1E-12||pmd>=1.0-pfa)
        return 0;
    hi=dof;
    while (gammap(dof*0.5,hi*0.5)<1.0-pfa&&hi<1E6) hi*=2.0;
    if (hi>=1E6) return 0;
    for (i=0;i<70;i++) {
        mid=(lo+hi)*0.5;
        if (gammap(dof*0.5,mid*0.5)<1.0-pfa) lo=mid;
        else hi=mid;
    }
    tau=(lo+hi)*0.5;
    lo=0.0;
    hi=dof;
    while (ncchicdf(tau,dof,hi)>pmd&&hi<1E6) hi*=2.0;
    if (hi>=1E6) return 0;
    for (i=0;i<70;i++) {
        mid=(lo+hi)*0.5;
        if (ncchicdf(tau,dof,mid)>pmd) lo=mid;
        else hi=mid;
    }
    *thres=tau;
    *lam=(lo+hi)*0.5;
    return isfinite(*thres)&&isfinite(*lam);
}

#ifdef ENABLE_RTK_INTEGRITY

#define SQR_(x) ((x)*(x))
#define SQRT_(x) ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

#if ENABLE_RTK_RANSAC_SLOPE
/* Single-epoch float KF innovation slope; conditioned on the selected model.
 * A physical SD fault affects every DD row containing that satellite/band/type. */
static void calcplslkf(rtk_t *rtk)
{
    rtkim_t *mon=&rtk->intg;
    rtkim_pl_t *pl=&mon->pl;
    int keys[RTKIM_MAX_BIAS_ROWS*2],n=mon->bnv,nm=0,i,j,k,a,key,sat,f,code;
    double c[RTKIM_MAX_BIAS_ROWS],g[3],enu[3],pos[3],E[9],q[6];
    double tau=0.0,lam=0.0,test=0.0,hs=0.0,vs=0.0,pfa,pmd,kn;
    static int nv0=0,nm0=0;
    static double pfa0=0.0,pmd0=0.0,tau0=0.0,lam0=0.0;
    memset(pl,0,sizeof(*pl));
    pl->hpl=pl->vpl=INFINITY;
    if (!mon->slready||mon->slnv!=n||n<=0||!mon->bv||
        rtk->sol.stat!=SOLQ_FLOAT||rtk->opt.modear!=ARMODE_OFF||
        rtk->opt.niter!=1)
        goto done;
    for (i=0;i<n;i++) {
        int flag=mon->bflg[i],ref=(flag>>16)&255,
            target=(flag>>8)&255;
        f=flag&15;
        code=(flag>>4)&15;
        if (f>=NFREQ||code>1) continue;
        for (a=0;a<2;a++) {
            sat=a?target:ref;
            if (sat<=0||sat>MAXSAT) continue;
            key=(sat*NFREQ+f)*2+code;
            for (j=0;j<nm;j++)
                if (keys[j]==key) break;
            if (j==nm) keys[nm++]=key;
        }
    }
    if (nm<=0) goto done;
    pfa=rtkimprob(rtk->opt.intpfa,RTKINT_DEFAULT_PFA);
    pmd=rtkimprob(rtk->opt.intpmd,RTKINT_DEFAULT_PMD);
    /* Equal conservative allocation to single-fault modes and nominal noise. */
    if (nv0!=n||nm0!=nm||pfa0!=pfa||pmd0!=pmd) {
        if (!rtkslopethres(n,pfa,pmd/(nm+1),&tau,&lam)) goto done;
        nv0=n;
        nm0=nm;
        pfa0=pfa;
        pmd0=pmd;
        tau0=tau;
        lam0=lam;
    }
    tau=tau0;
    lam=lam0;
    for (i=0;i<n;i++)
        for (j=0;j<n;j++) test+=mon->bv[i]*mon->slQi[i+j*n]*mon->bv[j];
    if (!isfinite(test)||test>tau)
        goto done; /* detected epoch: unavailable PL, no FDE */
    q[0]=rtk->P[0];
    q[1]=rtk->P[1+rtk->nx];
    q[2]=rtk->P[2+2*rtk->nx];
    q[3]=rtk->P[1];
    q[4]=rtk->P[2+rtk->nx];
    q[5]=rtk->P[2];
    rtkimsigenu(rtk->x,q,pl->msig);
    ecef2pos(rtk->x,pos);
    xyz2enu(pos,E);
    for (a=0;a<nm;a++) {
        double d=0.0,h,v;
        key=keys[a];
        code=key%2;
        f=(key/2)%NFREQ;
        sat=key/(2*NFREQ);
        for (i=0;i<n;i++) {
            int flag=mon->bflg[i];
            c[i]=((flag&15)==f&&((flag>>4)&15)==code)
                       ?(((flag>>16)&255)==sat)-
                             (((flag>>8)&255)==sat)
                       :0.0;
        }
        for (k=0;k<3;k++) {
            g[k]=0.0;
            for (i=0;i<n;i++) g[k]+=mon->slK[k+i*3]*c[i];
        }
        matmul("NN",3,1,3,1.0,E,g,0.0,enu);
        for (i=0;i<n;i++)
            for (j=0;j<n;j++) d+=c[i]*mon->slQi[i+j*n]*c[j];
        if (!isfinite(d)||d<=0.0) goto done;
        h=hypot(enu[0],enu[1])/sqrt(d);
        v=fabs(enu[2])/sqrt(d);
        if (!isfinite(h)||!isfinite(v)) goto done;
        if (h>hs) {
            hs=h;
            pl->hsrc=a+1;
            pl->hsat=sat;
            pl->hmode=RTKIM_F_SAT;
        }
        if (v>vs) {
            vs=v;
            pl->vsrc=a+1;
            pl->vsat=sat;
            pl->vmode=RTKIM_F_SAT;
        }
    }
    kn=rtkimnorminv(1.0-pmd/(4.0*(nm+1)));
    pl->hpl0=hs*sqrt(lam);
    pl->vpl0=vs*sqrt(lam);
    pl->hadd=kn*hypot(pl->msig[0],pl->msig[1]);
    pl->vadd=kn*pl->msig[2];
    pl->hpl=pl->hpl0+pl->hadd;
    pl->vpl=pl->vpl0+pl->vadd;
done:
    rtk->inthpl=pl->hpl;
    rtk->intvpl=pl->vpl;
    rtkimdbg(
        "SLOPE_KF %s %d %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
        time_str(rtk->sol.time,3),isfinite(pl->hpl)&&isfinite(pl->vpl),n,
        nm,test,tau,lam,hs,vs,pl->hpl,pl->vpl);
}
#endif

void rtkimslassess(rtk_t *rtk)
{
#if ENABLE_RTK_RANSAC_SLOPE
    calcplslkf(rtk);
    return;
#endif
    rtkim_t *mon=&rtk->intg;
    rtkim_pl_t *pl=&mon->pl;
    double qmain[6],msig[3],kmain=0.0,ksub=0.0,kss=0.0;
    double P[9],pos[3],E[9],*F=NULL,*Q=NULL,*Qinv=NULL,
                               *Kecef=NULL,*Kenu=NULL,*A=NULL,
                               *S=NULL;
    double hsmax=0.0,vsmax=0.0,pbias=0.0,sigmab=0.0;
    int i,j,k,n=mon->bnv,ih= -1,iv= -1;

    memset(pl,0,sizeof(*pl));
    if (n<=0||!mon->bH||!mon->bR||!rtk->x||!rtk->P||rtk->nx<3) {
        rtk->inthpl=0.0;
        rtk->intvpl=0.0;
        return;
    }
    rtkimthres(&rtk->opt,n,&kmain,&ksub,&kss);
    qmain[0]=rtk->P[0];
    qmain[1]=rtk->P[1+rtk->nx];
    qmain[2]=rtk->P[2+2*rtk->nx];
    qmain[3]=rtk->P[1];
    qmain[4]=rtk->P[2+rtk->nx];
    qmain[5]=rtk->P[2];
    rtkimsigenu(rtk->x,qmain,msig);
    for (i=0;i<3;i++) pl->msig[i]=msig[i];

    P[0]=rtk->P[0];
    P[1]=rtk->P[1];
    P[2]=rtk->P[2];
    P[3]=rtk->P[rtk->nx];
    P[4]=rtk->P[1+rtk->nx];
    P[5]=rtk->P[2+rtk->nx];
    P[6]=rtk->P[2*rtk->nx];
    P[7]=rtk->P[1+2*rtk->nx];
    P[8]=rtk->P[2+2*rtk->nx];

    ecef2pos(rtk->x,pos);
    xyz2enu(pos,E);

    F=mat(3,n);
    Q=mat(n,n);
    Qinv=mat(n,n);
    Kecef=mat(3,n);
    Kenu=mat(3,n);
    A=eye(n);
    S=mat(n,n);
    if (!F||!Q||!Qinv||!Kecef||!Kenu||!A||!S) {
        free(F);
        free(Q);
        free(Qinv);
        free(Kecef);
        free(Kenu);
        free(A);
        free(S);
        rtk->inthpl=0.0;
        rtk->intvpl=0.0;
        return;
    }
    for (j=0;j<n;j++) {
        for (i=0;i<3;i++) {
            F[i+j*3]=0.0;
            for (k=0;k<3;k++) F[i+j*3]+=P[i+k*3]*mon->bH[k+j*mon->bnx];
        }
    }
    for (j=0;j<n;j++)
        for (i=0;i<n;i++) {
            double s=mon->bR[i+j*n];
            for (k=0;k<3;k++) s+=mon->bH[k+i*mon->bnx]*F[k+j*3];
            Q[i+j*n]=s;
        }
    matcpy(Qinv,Q,n,n);
    if (!matinv(Qinv,n)) {
        matmul("NN",3,n,n,1.0,F,Qinv,0.0,Kecef);
        for (j=0;j<n;j++) {
            for (i=0;i<3;i++) {
                Kenu[i+j*3]=0.0;
                for (k=0;k<3;k++) Kenu[i+j*3]+=E[i+k*3]*Kecef[k+j*3];
            }
        }
        for (j=0;j<n;j++)
            for (i=0;i<n;i++) {
                double hk=0.0;
                for (k=0;k<3;k++) hk+=mon->bH[k+i*mon->bnx]*Kecef[k+j*3];
                A[i+j*n]-=hk;
            }
        matmul("TN",n,n,n,1.0,A,A,0.0,S);
        for (j=0;j<n;j++) {
            double sii=S[j+j*n],hs,vs;
            if (sii<=0.0) continue;
            hs=SQRT_(SQR_(Kenu[0+j*3])+SQR_(Kenu[1+j*3]))/
                 sqrt(sii);
            vs=fabs(Kenu[2+j*3])/sqrt(sii);
            if (hs>hsmax) {
                hsmax=hs;
                ih=j;
            }
            if (vs>vsmax) {
                vsmax=vs;
                iv=j;
            }
        }
        if (ih>=0) {
            for (j=0;j<n;j++) sigmab+=mon->bR[ih+j*n];
            sigmab=SQRT_(fabs(sigmab));
            pbias=kss*sigmab;
            pl->hsrc=ih+1;
            pl->hsat=rtkimsat(mon->bflg[ih]);
            pl->hmode=RTKIM_F_SAT;
        }
        if (iv>=0) {
            pl->vsrc=iv+1;
            pl->vsat=rtkimsat(mon->bflg[iv]);
            pl->vmode=RTKIM_F_SAT;
        }
        pl->hpl0=hsmax*pbias;
        pl->vpl0=vsmax*pbias;
    }
    free(F);
    free(Q);
    free(Qinv);
    free(Kecef);
    free(Kenu);
    free(A);
    free(S);
    pl->pe=pl->hpl0;
    pl->pn=0.0;
    pl->pu=pl->vpl0;
    pl->hpl=pl->hpl0;
    pl->vpl=pl->vpl0;
    pl->hadd=pl->hpl-pl->hpl0;
    pl->vadd=pl->vpl-pl->vpl0;
    rtk->inthpl=pl->hpl;
    rtk->intvpl=pl->vpl;
}

/* Mirror filter()'s active-state selection before its measurement update.
 * Capturing here also avoids using R after the legacy MDB code inverts it. */
EXPORT void rtkimslprior(rtk_t *rtk, const double *x, const double *P,
                         const double *H, const double *R, int nx, int nv)
{
#if ENABLE_RTK_RANSAC_SLOPE
    rtkim_t *mon=&rtk->intg;
    int ix[MAXSAT*NFREQ*2],i,j,k,n=0;
    double *Pc=NULL,*Hc=NULL,*F=NULL,*Qi=NULL,*Kp=NULL;
    mon->slready=0;
    if (!rtk->opt.intmon||rtk->intchild||!x||!P||!H||!R||nv<=0||
        nv>RTKIM_MAX_BIAS_ROWS||nx<=0||nx>MAXSAT*NFREQ*2||
        mon->bnv!=nv)
        return;
    for (i=0;i<nx;i++)
        if (x[i]!=0.0&&P[i+i*nx]>0.0) ix[n++]=i;
    if (n<3) return;
    Pc=mat(n,n);
    Hc=mat(n,nv);
    F=mat(n,nv);
    Qi=mat(nv,nv);
    Kp=zeros(3,nv);
    if (!Pc||!Hc||!F||!Qi||!Kp) goto done;
    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) Pc[i+j*n]=P[ix[i]+ix[j]*nx];
        for (j=0;j<nv;j++) Hc[i+j*n]=H[ix[i]+j*nx];
    }
    matmul("NN",n,nv,n,1.0,Pc,Hc,0.0,F);
    matcpy(Qi,R,nv,nv);
    matmul("TN",nv,nv,n,1.0,Hc,F,1.0,Qi);
    if (matinv(Qi,nv)) goto done;
    for (i=0;i<n;i++)
        if (ix[i]<3) {
            for (j=0;j<nv;j++)
                for (k=0;k<nv;k++) Kp[ix[i]+j*3]+=F[i+k*n]*Qi[k+j*nv];
        }
    for (i=0;i<3*nv;i++)
        if (!isfinite(Kp[i])) goto done;
    for (i=0;i<nv*nv;i++)
        if (!isfinite(Qi[i])) goto done;
    free(mon->slK);
    free(mon->slQi);
    mon->slK=Kp;
    mon->slQi=Qi;
    mon->slnv=nv;
    mon->slready=1;
    Kp=Qi=NULL;
done:
    free(Pc);
    free(Hc);
    free(F);
    free(Qi);
    free(Kp);
#else
    (void)rtk;
    (void)x;
    (void)P;
    (void)H;
    (void)R;
    (void)nx;
    (void)nv;
#endif
}

#endif
