#include "rtklib.h"

#ifdef ENABLE_RTK_INTEGRITY

#define SQR_(x)      ((x)*(x))
#define SQRT_(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define K_MAIN       5.33
#define K_SUB        5.33
#define K_SS         5.10
#define NF_(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP_(opt)     ((opt)->dynamics==0?3:9)
#define NI_(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT_(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL_(opt)     ((opt)->glomodear!=GLO_ARMODE_AUTOCAL?0:NFREQGLO)
#define NR_(opt)     (NP_(opt)+NI_(opt)+NT_(opt)+NL_(opt))
#define IB_(s,f,opt) (NR_(opt)+MAXSAT*(f)+(s)-1)

static FILE *fp_pint=NULL,*fp_pld=NULL,*fp_sub=NULL,*fp_rbias=NULL;
static int out_sub=0,out_rbias=0;

static void sigenu(const double *rr, const double *qr, double *sig)
{
    double pos[3],q[9],qe[9];
    ecef2pos(rr,pos);
    q[0]=qr[0]; q[4]=qr[1]; q[8]=qr[2];
    q[1]=q[3]=qr[3]; q[5]=q[7]=qr[4]; q[2]=q[6]=qr[5];
    covenu(pos,q,qe);
    sig[0]=SQRT_(qe[0]);
    sig[1]=SQRT_(qe[4]);
    sig[2]=SQRT_(qe[8]);
}

static void sepenu(const double *base, const double *pos, double *sep)
{
    double bl[3],llh[3],E[9];
    bl[0]=pos[0]-base[0]; bl[1]=pos[1]-base[1]; bl[2]=pos[2]-base[2];
    ecef2pos(base,llh);
    xyz2enu(llh,E);
    matmul("NN",3,1,3,1.0,E,bl,0.0,sep);
}

static void setpath(char *dst, const char *src, const char *ext)
{
    char *p;
    strcpy(dst,src);
    if ((p=strrchr(dst,'.'))) *p='\0';
    strcat(dst,ext);
}

static int hasobs(const obsd_t *obs, int n, int rcv, int sat)
{
    int i;
    for (i=0;i<n;i++) if (obs[i].rcv==rcv&&obs[i].sat==sat) return 1;
    return 0;
}

static int commonsat(const obsd_t *obs, int n, int *sat)
{
    int i,ns=0;
    for (i=0;i<n;i++) {
        if (obs[i].rcv!=1) continue;
        if (!hasobs(obs,n,2,obs[i].sat)) continue;
        sat[ns++]=obs[i].sat;
        if (ns>=MAXSAT) break;
    }
    return ns;
}

static void countobs(const obsd_t *obs, int n, int *nu, int *nr)
{
    int i;
    *nu=*nr=0;
    for (i=0;i<n;i++) {
        if (obs[i].rcv==1) (*nu)++;
        else if (obs[i].rcv==2) (*nr)++;
    }
}

static int finddef(rtkint_t *mon, int mode, int sat)
{
    int i;
    for (i=0;i<mon->ndef;i++) {
        if (mon->def[i].mode==mode&&mon->def[i].sat==sat) return i;
    }
    return -1;
}

static int adddef(rtk_t *rtk, int mode, int sat)
{
    rtkint_t *mon=&rtk->intg;
    prcopt_t *opt=&rtk->opt;
    int max=opt->rtk_integrity_max_subset_filters;
    int i=mon->ndef;
    if (max<=0) return -1;
    if (i>=max||i>=MAXSAT+1) return -1;
    mon->def[i].id=i+1;
    mon->def[i].mode=mode;
    mon->def[i].sat=sat;
    mon->def[i].fixed=mode==RTKINT_F_NFIX;
    mon->ndef++;
    return i;
}

static int getdef(rtk_t *rtk, int mode, int sat)
{
    int i=finddef(&rtk->intg,mode,sat);
    return i>=0?i:adddef(rtk,mode,sat);
}

static void freechild(rtkint_ss_t *ss)
{
    if (!ss->rtk) return;
    rtkfree(ss->rtk);
    free(ss->rtk);
    ss->rtk=NULL;
}

static void initchild(rtk_t *rtk, rtkint_ss_t *ss)
{
    prcopt_t opt=rtk->opt;
    if (ss->rtk) return;
    opt.enable_rtk_integrity_monitor=0;
    opt.modear=ARMODE_OFF;
    ss->rtk=(rtk_t *)malloc(sizeof(rtk_t));
    if (!ss->rtk) return;
    rtkinit(ss->rtk,&opt);
    ss->rtk->int_child=1;
}

static void filtobs(const obsd_t *obs, int n, int sat, obsd_t *out, int *m)
{
    int i;
    *m=0;
    for (i=0;i<n;i++) {
        if (sat>0&&obs[i].sat==sat) continue;
        out[(*m)++]=obs[i];
    }
}

static void setsolfromfilter(rtk_t *rtk)
{
    int i,nx=rtk->nx;
    for (i=0;i<3;i++) {
        rtk->sol.rr[i]=rtk->x[i];
        rtk->sol.qr[i]=(float)rtk->P[i+i*nx];
    }
    rtk->sol.qr[3]=(float)rtk->P[1];
    rtk->sol.qr[4]=(float)rtk->P[1+2*nx];
    rtk->sol.qr[5]=(float)rtk->P[2];
}

static void condamb(rtk_t *main, rtkint_def_t *def, rtk_t *child,
                    rtkint_ss_t *ss)
{
    int idx[MAXSAT*NFREQ];
    int i,f,j,k,l,nv=0,nx=child->nx,nf=NF_(&child->opt);

    ss->ambc=0;
    if (main->sol.stat!=SOLQ_FIX||def->mode==RTKINT_F_NFIX) return;
    if (main->na<=0||main->na!=child->na||main->nx!=child->nx) return;

    for (i=1;i<=MAXSAT&&nv<MAXSAT*NFREQ;i++) for (f=0;f<nf;f++) {
        j=IB_(i,f,&child->opt);
        if (j<0||j>=nx||j>=main->na) continue;
        if (def->mode==RTKINT_F_SAT&&i==def->sat) continue;
        if (main->xa[j]==0.0||child->x[j]==0.0) continue;
        if (main->Pa[j+j*main->na]<=0.0||child->P[j+j*nx]<=0.0) continue;
        idx[nv++]=j;
        if (nv>=MAXSAT*NFREQ) break;
    }
    if (nv<=0) return;

    for (i=0;i<nv;i++) {
        double v,r,s;
        double *kvec,*prow;
        j=idx[i];
        v=main->xa[j]-child->x[j];
        r=main->Pa[j+j*main->na]>1E-6?main->Pa[j+j*main->na]:1E-6;
        s=child->P[j+j*nx]+r;
        if (s<=0.0) continue;
        kvec=(double *)malloc(sizeof(double)*nx);
        prow=(double *)malloc(sizeof(double)*nx);
        if (!kvec||!prow) {
            free(kvec); free(prow);
            continue;
        }
        for (k=0;k<nx;k++) kvec[k]=child->P[k+j*nx]/s;
        for (l=0;l<nx;l++) prow[l]=child->P[j+l*nx];
        for (k=0;k<nx;k++) child->x[k]+=kvec[k]*v;
        for (k=0;k<nx;k++) for (l=0;l<nx;l++) {
            child->P[k+l*nx]-=kvec[k]*prow[l];
        }
        free(kvec); free(prow);
        ss->ambc++;
    }
    if (ss->ambc>0) setsolfromfilter(child);
}

static void runsub(rtk_t *rtk, int idx, const obsd_t *obs, int n,
                   const nav_t *nav, const sta_t *sta)
{
    rtkint_t *mon=&rtk->intg;
    rtkint_def_t *def=mon->def+idx;
    rtkint_ss_t *ss=mon->ss+idx;
    obsd_t ob[MAXOBS*2];
    int m=0,ok,i;

    ss->act=1;
    initchild(rtk,ss);
    if (!ss->rtk) return;
    ss->init=1;
    ss->rtk->opt=rtk->opt;
    ss->rtk->opt.enable_rtk_integrity_monitor=0;
    ss->rtk->opt.modear=ARMODE_OFF;
    ss->rtk->int_child=1;
    filtobs(obs,n,def->mode==RTKINT_F_SAT?def->sat:0,ob,&m);
    ok=rtkpos(ss->rtk,ob,m,nav,sta);
    ss->ready=1;
    ss->valid=ok&&ss->rtk->sol.stat!=SOLQ_NONE;
    if (!ss->valid) return;
    condamb(rtk,def,ss->rtk,ss);
    for (i=0;i<3;i++) ss->pos[i]=ss->rtk->sol.rr[i];
    for (i=0;i<6;i++) ss->q[i]=ss->rtk->sol.qr[i];
    sigenu(ss->pos,ss->q,ss->sig);
    sepenu(rtk->sol.rr,ss->pos,ss->sep);
    ss->qi=ss->rtk->sol.stat;
    ss->ratio=ss->rtk->sol.ratio;
    ss->pdop=ss->rtk->sol.dops[1];
    ss->hdop=ss->rtk->sol.dops[2];
    ss->vdop=ss->rtk->sol.dops[3];
    ss->sig0=ss->rtk->sol.Ftestvalue;
    ss->chisq=ss->rtk->sol.Ftestvalue;
    ss->dof=ss->rtk->sol.numofnv;
}

static void calcpl(rtk_t *rtk)
{
    rtkint_t *mon=&rtk->intg;
    rtkint_pl_t *pl=&mon->pl;
    double msig[3],pe,pn,pu,ssig[3],sss,subp,qmain[6];
    int i,k;

    memset(pl,0,sizeof(*pl));
    for (i=0;i<6;i++) qmain[i]=rtk->sol.qr[i];
    sigenu(rtk->sol.rr,qmain,msig);
    for (i=0;i<3;i++) pl->msig[i]=msig[i];
    pe=K_MAIN*msig[0]; pn=K_MAIN*msig[1]; pu=K_MAIN*msig[2];
    pl->hsrc=pl->vsrc=0;
    for (i=0;i<mon->ndef;i++) {
        rtkint_ss_t *ss=mon->ss+i;
        rtkint_def_t *def=mon->def+i;
        if (!ss->act||!ss->valid) continue;
        ssig[0]=ss->sig[0]; ssig[1]=ss->sig[1]; ssig[2]=ss->sig[2];
        for (k=0;k<3;k++) {
            sss=SQRT_(SQR_(ssig[k])-SQR_(msig[k]));
            subp=K_SUB*ssig[k]+K_SS*sss;
            if (k==0&&subp>pe) {
                pe=subp; pl->hsrc=def->id; pl->hmode=def->mode; pl->hsat=def->sat;
                pl->hsig[0]=ssig[0]; pl->hsig[1]=ssig[1]; pl->hsig[2]=ssig[2];
                pl->hsep[0]=ss->sep[0]; pl->hsep[1]=ss->sep[1]; pl->hsep[2]=ss->sep[2];
            }
            if (k==1&&subp>pn) {
                pn=subp; pl->hsrc=def->id; pl->hmode=def->mode; pl->hsat=def->sat;
                pl->hsig[0]=ssig[0]; pl->hsig[1]=ssig[1]; pl->hsig[2]=ssig[2];
                pl->hsep[0]=ss->sep[0]; pl->hsep[1]=ss->sep[1]; pl->hsep[2]=ss->sep[2];
            }
            if (k==2&&subp>pu) {
                pu=subp; pl->vsrc=def->id; pl->vmode=def->mode; pl->vsat=def->sat;
                pl->vsig[0]=ssig[0]; pl->vsig[1]=ssig[1]; pl->vsig[2]=ssig[2];
                pl->vsep[0]=ss->sep[0]; pl->vsep[1]=ss->sep[1]; pl->vsep[2]=ss->sep[2];
            }
        }
    }
    pl->pe=pe; pl->pn=pn; pl->pu=pu;
    pl->hpl0=SQRT_(SQR_(pe)+SQR_(pn));
    pl->vpl0=pu;
    pl->hpl=pl->hpl0;
    pl->vpl=pl->vpl0;
    rtk->int_hpl=pl->hpl;
    rtk->int_vpl=pl->vpl;
}

static void evalfde(rtk_t *rtk)
{
    rtkint_t *mon=&rtk->intg;
    rtkint_fde_t best={0};
    int i,k;
    if (!rtk->opt.enable_rtk_integrity_fde_recovery) {
        mon->fde=best;
        return;
    }
    if (rtk->sol.stat==SOLQ_NONE||rtk->sol.Ftestvalue<=4.0) {
        mon->fde=best;
        return;
    }
    for (i=0;i<mon->ndef;i++) {
        rtkint_ss_t *ss=mon->ss+i;
        rtkint_def_t *def=mon->def+i;
        double score=0.0;
        if (!ss->act||!ss->valid) continue;
        for (k=0;k<3;k++) {
            double th=K_SS*SQRT_(SQR_(ss->sig[k])-SQR_(mon->pl.msig[k]));
            double st=fabs(ss->sep[k]);
            double sc=th>0.0?st/th:0.0;
            if (sc>score) score=sc;
        }
        if (score>best.score) {
            best.mode=def->mode;
            best.id=def->id;
            best.sat=def->sat;
            best.score=score;
            best.act=def->mode==RTKINT_F_SAT?RTKINT_A_SAT:
                     def->mode==RTKINT_F_NFIX?RTKINT_A_FLOAT:RTKINT_A_NONE;
        }
    }
    if (best.score<=1.0) memset(&best,0,sizeof(best));
    mon->fde=best;
    mon->act=best;
    rtk->int_fde_mode=best.mode;
    rtk->int_fde_sat=best.sat;
}

static void subout(const rtk_t *rtk)
{
    const rtkint_t *mon=&rtk->intg;
    int i;
    if (!fp_sub||!rtk->opt.enable_rtk_integrity_subset_debug_output) return;
    for (i=0;i<mon->ndef;i++) {
        const rtkint_def_t *def=mon->def+i;
        const rtkint_ss_t *ss=mon->ss+i;
        if (!ss->act) continue;
        if (rtk->opt.rtk_integrity_debug_subset_id>0&&
            rtk->opt.rtk_integrity_debug_subset_id!=def->id) continue;
        if (rtk->opt.rtk_integrity_debug_satellite>0&&
            rtk->opt.rtk_integrity_debug_satellite!=def->sat) continue;
        fprintf(fp_sub,"%s %d %d %d %d %d %.4f %.4f %.4f %.4f %.4f %.4f\n",
            time_str(mon->ep.time,3),def->id,def->mode,def->sat,ss->valid,
            ss->ambc,ss->pos[0],ss->pos[1],ss->pos[2],ss->sep[0],ss->sep[1],
            ss->sep[2]);
    }
}

EXPORT void rtkint_init(rtk_t *rtk, const prcopt_t *opt)
{
    memset(&rtk->intg,0,sizeof(rtk->intg));
    rtk->intg.ena=opt->enable_rtk_integrity_monitor;
    rtk->int_reproc=0;
    rtk->int_fde_mode=0;
    rtk->int_fde_sat=0;
    rtk->int_hpl=rtk->int_vpl=0.0;
}

EXPORT void rtkint_free(rtk_t *rtk)
{
    int i;
    for (i=0;i<rtk->intg.ndef;i++) freechild(rtk->intg.ss+i);
    memset(&rtk->intg,0,sizeof(rtk->intg));
}

EXPORT void rtkint_update(rtk_t *rtk, const obsd_t *obs, int nobs,
                          const nav_t *nav, const sta_t *sta)
{
    rtkint_t *mon=&rtk->intg;
    int sat[MAXSAT],ns,i,nu,nr;
    if (!rtk->opt.enable_rtk_integrity_monitor||rtk->int_child) return;
    mon->ena=mon->init=1;
    for (i=0;i<mon->ndef;i++) mon->ss[i].act=0;
    ns=commonsat(obs,nobs,sat);
    countobs(obs,nobs,&nu,&nr);
    mon->ep.time=rtk->sol.time;
    mon->ep.nu=nu;
    mon->ep.nr=nr;
    mon->ep.ns=ns;
    mon->ep.qi=rtk->sol.stat;
    mon->ep.ratio=rtk->sol.ratio;
    mon->ep.fixed=rtk->sol.stat==SOLQ_FIX;
    mon->ep.fixed_upd=mon->ep.fixed;
    if (rtk->opt.enable_monitor_single_satellite_fault) {
        for (i=0;i<ns;i++) getdef(rtk,RTKINT_F_SAT,sat[i]);
    }
    if (mon->ep.fixed&&rtk->opt.enable_monitor_fixed_ambiguity_update_fault) {
        getdef(rtk,RTKINT_F_NFIX,0);
    }
    for (i=0;i<mon->ndef;i++) {
        if (mon->def[i].mode==RTKINT_F_SAT&&!hasobs(obs,nobs,1,mon->def[i].sat)) continue;
        if (mon->def[i].mode==RTKINT_F_NFIX&&!mon->ep.fixed) continue;
        runsub(rtk,i,obs,nobs,nav,sta);
    }
    mon->nact=0;
    for (i=0;i<mon->ndef;i++) if (mon->ss[i].act) mon->nact++;
    calcpl(rtk);
    evalfde(rtk);
}

EXPORT int rtkint_redo(const rtk_t *rtk, int *mode, int *sat)
{
    const rtkint_t *mon=&rtk->intg;
    if (!rtk->opt.enable_rtk_integrity_fde_recovery) return 0;
    if (mon->fde.act==RTKINT_A_NONE) return 0;
    *mode=mon->fde.act;
    *sat=mon->fde.sat;
    return 1;
}

EXPORT int rtkint_open(const char *outfile, const prcopt_t *opt)
{
    char path[1024];
    if (!opt->enable_rtk_integrity_monitor||!outfile||!*outfile) return 0;
    setpath(path,outfile,".pint");
    fp_pint=fopen(path,"w");
    setpath(path,outfile,".pldiag");
    fp_pld=fopen(path,"w");
    out_sub=opt->enable_rtk_integrity_subset_debug_output;
    out_rbias=opt->enable_rtk_integrity_rbias_export;
    if (out_sub) {
        setpath(path,outfile,".subdiag");
        fp_sub=fopen(path,"w");
    }
    if (out_rbias) {
        setpath(path,outfile,".rbias");
        fp_rbias=fopen(path,"w");
    }
    if (fp_pint) fprintf(fp_pint,"%% time x y z nsat QI ratio sigma0 subset_total subset_active HPL VPL FDEmode FDEsat\n");
    if (fp_pld) fprintf(fp_pld,"%% time nsat QI ratio subset_total subset_active HPL VPL bias_e bias_n bias_u HPL0 VPL0 Hsrc Hmode Hsat HsigE HsigN HsigU HsepE HsepN HsepU Vsrc Vmode Vsat VsigE VsigN VsigU VsepE VsepN VsepU\n");
    if (fp_sub) fprintf(fp_sub,"%% time subset mode sat valid amb_cond x y z sep_e sep_n sep_u\n");
    if (fp_rbias) fprintf(fp_rbias,"%% time sat res bias_e bias_n bias_u\n");
    return fp_pint||fp_pld;
}

EXPORT void rtkint_close(void)
{
    if (fp_pint) fclose(fp_pint); fp_pint=NULL;
    if (fp_pld) fclose(fp_pld); fp_pld=NULL;
    if (fp_sub) fclose(fp_sub); fp_sub=NULL;
    if (fp_rbias) fclose(fp_rbias); fp_rbias=NULL;
}

EXPORT void rtkint_out(const rtk_t *rtk)
{
    const rtkint_t *mon=&rtk->intg;
    const rtkint_pl_t *pl=&mon->pl;
    if (!rtk->opt.enable_rtk_integrity_monitor||!mon->init) return;
    if (fp_pint) {
        fprintf(fp_pint,"%s %.4f %.4f %.4f %d %d %.3f %.4f %d %d %.4f %.4f %d %d\n",
            time_str(mon->ep.time,3),rtk->sol.rr[0],rtk->sol.rr[1],rtk->sol.rr[2],
            rtk->sol.ns,rtk->sol.stat,rtk->sol.ratio,rtk->sol.Ftestvalue,
            mon->ndef,mon->nact,pl->hpl,pl->vpl,mon->fde.mode,mon->fde.sat);
    }
    if (fp_pld) {
        fprintf(fp_pld,"%s %d %d %.3f %d %d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %d %d %d %.4f %.4f %.4f %.4f %.4f %.4f %d %d %d %.4f %.4f %.4f %.4f %.4f %.4f\n",
            time_str(mon->ep.time,3),rtk->sol.ns,rtk->sol.stat,rtk->sol.ratio,
            mon->ndef,mon->nact,pl->hpl,pl->vpl,pl->be,pl->bn,pl->bu,pl->hpl0,
            pl->vpl0,pl->hsrc,pl->hmode,pl->hsat,pl->hsig[0],pl->hsig[1],
            pl->hsig[2],pl->hsep[0],pl->hsep[1],pl->hsep[2],pl->vsrc,
            pl->vmode,pl->vsat,pl->vsig[0],pl->vsig[1],pl->vsig[2],
            pl->vsep[0],pl->vsep[1],pl->vsep[2]);
    }
    subout(rtk);
}

#endif
