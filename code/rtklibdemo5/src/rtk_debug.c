/*------------------------------------------------------------------------------
* rtk_debug.c : RTK diagnostic output and ambiguity Monte Carlo
*
* notes  : Diagnostic streams are owned by this module.
*-----------------------------------------------------------------------------*/
#include "rtk_int_internal.h"

#ifdef ENABLE_RTK_INTEGRITY

#if ENABLE_RTK_DEBUG_OUTPUT
static FILE *fp_rtkdbg=NULL;
static FILE *fp_vtest=NULL;
#endif

static void setpath(char *dst, const char *src, const char *ext)
{
    char *p;

    strcpy(dst,src);
    if ((p=strrchr(dst,'.'))) *p='\0';
    strcat(dst,ext);
}

/* Private logging bridge: monitoring code never owns diagnostic streams. */
void rtkimdbg(const char *format, ...)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    va_list args;
    if (!fp_rtkdbg) return;
    va_start(args,format);
    vfprintf(fp_rtkdbg,format,args);
    va_end(args);
#else
    (void)format;
#endif
}

EXPORT void rtkdbgrobust(gtime_t time, int nv, int iterations, double change,
                        int status, const char *reason)
{
    rtkimdbg("IGG3_ITER %s %d %d %.9g %d %s\n",time_str(time,3),nv,
             iterations,change,status,reason);
}

EXPORT void rtkdbgweight(gtime_t time,int row,int flag,double weight)
{
    char ref[16],sat[16];
    satno2id((flag>>16)&255,ref);
    satno2id((flag>>8)&255,sat);
    rtkimdbg("IGG3_WEIGHT %s %d %s %s %d %d %.17g\n",
             time_str(time,3),row,ref,sat,(flag>>4)&15,(flag&15)+1,weight);
}

EXPORT int rtkdbgopen(const char *outfile)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    char path[1024];
    if (!outfile||!*outfile) return 0;
    setpath(path,outfile,".rtk_debug");
    fp_rtkdbg=fopen(path,"w");
    if (fp_rtkdbg) {
        fprintf(fp_rtkdbg,
                "%% type time stage action sat kept_obs total_obs stat ns "
                "ratio chi2 nx na rr0 rr1 rr2 drr0 drr1 drr2\n");
        fprintf(fp_rtkdbg,"%% RAWCOUNTS time stage n nu nr\n");
        fprintf(fp_rtkdbg,"%% RAWSATSET time stage count sat_list\n");
        fprintf(fp_rtkdbg,
                "%% RAWLINE time stage satid preset_sat decoded_sat ok\n");
        fprintf(fp_rtkdbg,"%% SATSET time stage count sat_list\n");
        fprintf(fp_rtkdbg,"%% RANSAC_SD time sat freq model_dim obs_count "
                           "consensus_count keep sd_m variance_m2\n");
        fprintf(fp_rtkdbg,"%% IGG3_ITER time nv iterations max_log_weight_change "
                           "status reason\n");
        fprintf(fp_rtkdbg,
                "%% SLOPE_KF time available nv fault_modes statistic "
                "threshold noncentrality h_slope v_slope HPL VPL\n");
        fprintf(
            fp_rtkdbg,
            "%% AMB_FIX time ratio threshold dd_index dd_count ref_sat sat "
            "freq float_dd_cycles integer_dd_cycles residual_cycles\n");
        fprintf(fp_rtkdbg,"%% AMB_TRUTH time ref_sat sat freq "
                           "truth_constrained_dd_cycles\n");
        fprintf(fp_rtkdbg,"%% AMB_FLOAT id time nb index float_dd_cycles\n");
        fprintf(fp_rtkdbg,"%% AMB_COV id time nb row col q_cycles2\n");
        fprintf(fp_rtkdbg,"%% MC_EPOCH id time nb samples ratio_threshold\n");
        fprintf(fp_rtkdbg,"%% MC_MODE id time nb k1 ... kn count probability "
                           "dx dy dz norm\n");
        fprintf(fp_rtkdbg,"%% MC_RESULT id time nb samples accepted correct "
                           "wrong acceptance_rate conditional_correct_rate\n");
    }
    setpath(path,outfile,".vtest");
    fp_vtest=fopen(path,"w");
    if (fp_vtest) {
        fprintf(
            fp_vtest,
            "%% VAL time stage solstat ratio nv np dof chi2 threshold pass\n");
        fprintf(fp_vtest,"%% TOP time stage rank sat_ref sat_sat meas freq "
                          "resid sigma contrib\n");
    }
    return fp_rtkdbg!=NULL||fp_vtest!=NULL;
#else
    (void)outfile;
    return 0;
#endif
}

EXPORT void rtkdbgclose(void)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    if (fp_rtkdbg) fclose(fp_rtkdbg);
    fp_rtkdbg=NULL;
    if (fp_vtest) fclose(fp_vtest);
    fp_vtest=NULL;
#endif
}

EXPORT void rtkdbgvaltest(gtime_t time, const char *stage, int solstat,
                          float ratio, int nv, int np, int dof, double chi2,
                          double chi2_thres, int pass, const double *v,
                          const double *R, const int *vflg, int topn)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    int i,j,k,n,imax[8]={0},used[8]={0};
    double c[8]={0},sigma;
    const char *stype;
    int sat1,sat2,type,freq;
    if (!fp_vtest||!stage) return;
    fprintf(fp_vtest,"VAL %s %s %d %.3f %d %d %d %.4f %.4f %d\n",
            time_str(time,3),stage,solstat,ratio,nv,np,dof,chi2,
            chi2_thres,pass);
    if (!v||!R||!vflg||nv<=0) {
        fflush(fp_vtest);
        return;
    }
    n=topn<0?0:topn;
    if (n>8) n=8;
    for (i=0;i<n;i++) {
        double best= -1.0;
        int bestj= -1;
        for (j=0;j<nv;j++) {
            int skip=0;
            double rii=R[j+j*nv];
            if (rii<=0.0) continue;
            for (k=0;k<i;k++) {
                if (used[k]==j) {
                    skip=1;
                    break;
                }
            }
            if (skip) continue;
            c[j%8]=v[j]*v[j]/rii;
            if (c[j%8]>best) {
                best=c[j%8];
                bestj=j;
            }
        }
        if (bestj<0) break;
        used[i]=bestj;
        imax[i]=bestj;
        sat1=(vflg[bestj]>>16)&0xFF;
        sat2=(vflg[bestj]>>8)&0xFF;
        type=(vflg[bestj]>>4)&0xF;
        freq=(vflg[bestj]&0xF)+1;
        stype=type==0?"L":(type==1?"P":"C");
        sigma=
            R[bestj+bestj*nv]>0.0?sqrt(R[bestj+bestj*nv]):0.0;
        fprintf(fp_vtest,"TOP %s %s %d %d %d %s %d %.4f %.4f %.4f\n",
                time_str(time,3),stage,i+1,sat1,sat2,stype,freq,
                v[bestj],sigma,v[bestj]*v[bestj]/R[bestj+bestj*nv]);
    }
    fflush(fp_vtest);
#else
    (void)time;
    (void)stage;
    (void)solstat;
    (void)ratio;
    (void)nv;
    (void)np;
    (void)dof;
    (void)chi2;
    (void)chi2_thres;
    (void)pass;
    (void)v;
    (void)R;
    (void)vflg;
    (void)topn;
#endif
}

EXPORT void rtkdbgepoch(const rtk_t *rtk, const char *stage, int action,
                        int sat, int nkeep, int nobs, int stat0, int ns0,
                        float ratio0, double test0, const double *rr0)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    const sol_t *sol;
    double d0=0.0,d1=0.0,d2=0.0;
    if (!fp_rtkdbg||!rtk||!stage) return;
    sol=&rtk->sol;
    if (rr0) {
        d0=sol->rr[0]-rr0[0];
        d1=sol->rr[1]-rr0[1];
        d2=sol->rr[2]-rr0[2];
    }
    fprintf(fp_rtkdbg,
            "EPOCH %s %s %d %d %d %d %d %d %.3f %.4f %d %d %.4f %.4f %.4f "
            "%.4f %.4f %.4f\n",
            time_str(sol->time,3),stage,action,sat,nkeep,nobs,
            !strcmp(stage,"before")?stat0:sol->stat,
            !strcmp(stage,"before")?ns0:sol->ns,
            !strcmp(stage,"before")?ratio0:sol->ratio,
            !strcmp(stage,"before")?test0:sol->chi2testvalue,rtk->nx,
            rtk->na,sol->rr[0],sol->rr[1],sol->rr[2],d0,d1,d2);
    fflush(fp_rtkdbg);
#else
    (void)rtk;
    (void)stage;
    (void)action;
    (void)sat;
    (void)nkeep;
    (void)nobs;
    (void)stat0;
    (void)ns0;
    (void)ratio0;
    (void)test0;
    (void)rr0;
#endif
}

EXPORT void rtkdbgrawcounts(const char *stage, gtime_t time, int nobs, int nu,
                            int nr)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    FILE *fp=fp_rtkdbg;
    if (!fp||!stage) return;
    fprintf(fp,"RAWCOUNTS %s %s %d %d %d\n",time_str(time,3),stage,nobs,
            nu,nr);
    fflush(fp);
#else
    (void)stage;
    (void)time;
    (void)nobs;
    (void)nu;
    (void)nr;
#endif
}

EXPORT void rtkdbgrawset(const char *stage, gtime_t time, const obsd_t *obs,
                         int nobs)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    int i;
    char id[32];
    FILE *fp=fp_rtkdbg;
    if (!fp||!stage||!obs||nobs<=0) return;
    fprintf(fp,"RAWSATSET %s %s %d",time_str(time,3),stage,nobs);
    for (i=0;i<nobs;i++) {
        satno2id(obs[i].sat,id);
        fprintf(fp," %s",id);
    }
    fputc('\n',fp);
    fflush(fp);
#else
    (void)stage;
    (void)time;
    (void)obs;
    (void)nobs;
#endif
}

EXPORT void rtkdbgrawline(const char *stage, gtime_t time, const char *satid,
                          int preset_sat, int decoded_sat, int ok)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    FILE *fp=fp_rtkdbg;
    if (!fp||!stage) return;
    fprintf(fp,"RAWLINE %s %s %s %d %d %d\n",time_str(time,3),stage,
            satid&& *satid?satid:"---",preset_sat,decoded_sat,ok);
    fflush(fp);
#else
    (void)stage;
    (void)time;
    (void)satid;
    (void)preset_sat;
    (void)decoded_sat;
    (void)ok;
#endif
}

EXPORT void rtkdbgcounts(const rtk_t *rtk, const char *stage, int nobs, int nu,
                         int nr)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    if (!fp_rtkdbg||!rtk||!stage) return;
    fprintf(fp_rtkdbg,"COUNTS %s %s %d %d %d\n",time_str(rtk->sol.time,3),
            stage,nobs,nu,nr);
    fflush(fp_rtkdbg);
#else
    (void)rtk;
    (void)stage;
    (void)nobs;
    (void)nu;
    (void)nr;
#endif
}

EXPORT void rtkdbgransac(gtime_t time, int sat, int freq, int np, int nv,
                         int kept, int keep, double sd, double var)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    char id[16];
    if (!fp_rtkdbg) return;
    satno2id(sat,id);
    fprintf(fp_rtkdbg,"RANSAC_SD %s %s %d %d %d %d %d %.17g %.17g\n",
            time_str(time,3),id,freq,np,nv,kept,keep,sd,var);
#else
    (void)time;
    (void)sat;
    (void)freq;
    (void)np;
    (void)nv;
    (void)kept;
    (void)keep;
    (void)sd;
    (void)var;
#endif
}

EXPORT void rtkdbgambfix(gtime_t time, float ratio, float thres, int idx,
                         int nb, const char *refsat, const char *sat, int freq,
                         double af, double az)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    if (!fp_rtkdbg) return;
    fprintf(fp_rtkdbg,"AMB_FIX %s %.3f %.3f %d %d %s %s %d %.6f %.0f %.6f\n",
            time_str(time,3),ratio,thres,idx,nb,
            refsat&& *refsat?refsat:"---",sat&& *sat?sat:"---",
            freq,af,az,af-az);
    fflush(fp_rtkdbg);
#else
    (void)time;
    (void)ratio;
    (void)thres;
    (void)idx;
    (void)nb;
    (void)refsat;
    (void)sat;
    (void)freq;
    (void)af;
    (void)az;
#endif
}

EXPORT void rtkdbgambtruth(gtime_t time, const char *refsat, const char *sat,
                           int freq, double truth_cycles)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    if (!fp_rtkdbg) return;
    fprintf(fp_rtkdbg,"AMB_TRUTH %s %s %s %d %.6f\n",time_str(time,3),
            refsat&& *refsat?refsat:"---",sat&& *sat?sat:"---",
            freq,truth_cycles);
    fflush(fp_rtkdbg);
#else
    (void)time;
    (void)refsat;
    (void)sat;
    (void)freq;
    (void)truth_cycles;
#endif
}

/* Export the exact float DD ambiguity vector and covariance passed to LAMBDA. */
EXPORT void rtkdbgfloatamb(gtime_t time, int nb, const double *amb,
                           const double *Q)
{
#if ENABLE_RTK_DEBUG_OUTPUT
    static unsigned long long seq=0;
    unsigned long long id=++seq;
    int i,j;
    if (!fp_rtkdbg||!amb||!Q||nb<=0) return;
    for (i=0;i<nb;i++) {
        fprintf(fp_rtkdbg,"AMB_FLOAT %llu %s %d %d %.17g\n",id,
                time_str(time,3),nb,i+1,amb[i]);
    }
    for (i=0;i<nb;i++)
        for (j=0;j<nb;j++) {
            fprintf(fp_rtkdbg,"AMB_COV %llu %s %d %d %d %.17g\n",id,
                    time_str(time,3),nb,i+1,j+1,Q[i+j*nb]);
        }
    fflush(fp_rtkdbg);
#else
    (void)time;
    (void)nb;
    (void)amb;
    (void)Q;
#endif
}

/* Monte Carlo diagnostic for the final ambiguity set of one epoch. */
static void rtkdbgmcepochlegacy(gtime_t time, int nb, const double *Q,
                                int samples, double thres)
{
#if ENABLE_RTK_DEBUG_OUTPUT&&ENABLE_RTK_ONLINE_MC
    double *L=NULL,*eta=NULL,*F=NULL,s[2],u1,u2,gauss;
    int i,j,m,info,accepted=0,correct=0;
    static unsigned long long seq=0;
    unsigned long long id=++seq;
    if (!fp_rtkdbg||!Q||nb<=0||nb>MAXSAT*NFREQ||samples<=0) return;
    L=zeros(nb,nb);
    eta=zeros(nb,1);
    F=zeros(nb,2);
    if (!L||!eta||!F) {
        free(L);
        free(eta);
        free(F);
        return;
    }
    matcpy(L,Q,nb,nb);
    /* lower-triangular Cholesky factor L, Q=L*L' */
    for (i=0;i<nb;i++) {
        for (j=0;j<=i;j++) {
            double v=L[i+j*nb];
            int k;
            for (k=0;k<j;k++) v-=L[i+k*nb]*L[j+k*nb];
            if (i==j) {
                if (v<=0.0) {
                    free(L);
                    free(eta);
                    free(F);
                    return;
                }
                L[i+j*nb]=sqrt(v);
            } else L[i+j*nb]=v/L[j+j*nb];
        }
        for (j=i+1;j<nb;j++) L[i+j*nb]=0.0;
    }
    fprintf(fp_rtkdbg,"MC_EPOCH %llu %s %d %d %.6g\n",id,time_str(time,3),
            nb,samples,thres);
    for (m=0;m<samples;m++) {
        for (i=0;i<nb;i++) {
            double z=0.0;
            for (j=0;j<=i;j++) {
                u1=((double)rand()+1.0)/((double)RAND_MAX+2.0);
                u2=((double)rand()+1.0)/((double)RAND_MAX+2.0);
                gauss=sqrt(-2.0*log(u1))*cos(2.0*PI*u2);
                z+=L[i+j*nb]*gauss;
            }
            eta[i]=z;
        }
        info=lambda(nb,2,eta,Q,F,s);
        if (info||s[0]<=0.0||s[1]/s[0]<thres) continue;
        accepted++;
        {
            int iszero=1;
            for (i=0;i<nb;i++)
                if (F[i]!=0.0) {
                    iszero=0;
                    break;
                }
            if (iszero) correct++;
        }
        fprintf(fp_rtkdbg,"MC_OFFSET %llu %d",id,m+1);
        for (i=0;i<nb;i++) fprintf(fp_rtkdbg," %.0f",F[i]);
        fprintf(fp_rtkdbg,"\n");
    }
    fprintf(fp_rtkdbg,"MC_RESULT %llu %s %d %d %d %d %d %.17g %.17g\n",id,
            time_str(time,3),nb,samples,accepted,correct,
            accepted-correct,accepted/(double)samples,
            accepted?correct/(double)accepted:0.0);
    fflush(fp_rtkdbg);
    free(L);
    free(eta);
    free(F);
#else
    (void)time;
    (void)nb;
    (void)Q;
    (void)samples;
    (void)thres;
#endif
}

EXPORT void rtkdbgmcepoch(gtime_t time, int na, int nb, const double *Qab,
                          const double *Q, int samples, double thres)
{
#if ENABLE_RTK_DEBUG_OUTPUT&&ENABLE_RTK_ONLINE_MC
    enum { MC_MAX_MODES=20000 };
    double *L=NULL,*eta=NULL,*F=NULL,*Qi=NULL,*dv=NULL,s[2];
    int *modes=NULL,*counts=NULL,mode_n=0,i,j,m,info,accepted=0,correct=0;
    static unsigned long long seq=0;
    unsigned long long id=++seq;
    if (!fp_rtkdbg||!Q||nb<=0||nb>MAXSAT*NFREQ||samples<=0) return;
    L=zeros(nb,nb);
    eta=zeros(nb,1);
    F=zeros(nb,2);
    Qi=mat(nb,nb);
    dv=zeros(na>0?na:1,1);
    modes=imat(nb,MC_MAX_MODES);
    counts=imat(1,MC_MAX_MODES);
    if (!L||!eta||!F||!Qi||!dv||!modes||!counts) goto done;
    matcpy(L,Q,nb,nb);
    for (i=0;i<nb;i++) {
        for (j=0;j<=i;j++) {
            double v=L[i+j*nb];
            int k;
            for (k=0;k<j;k++) v-=L[i+k*nb]*L[j+k*nb];
            if (i==j) {
                if (v<=0.0) goto done;
                L[i+j*nb]=sqrt(v);
            } else L[i+j*nb]=v/L[j+j*nb];
        }
        for (j=i+1;j<nb;j++) L[i+j*nb]=0.0;
    }
    matcpy(Qi,Q,nb,nb);
    if (matinv(Qi,nb)) goto done;
    fprintf(fp_rtkdbg,"MC_EPOCH %llu %s %d %d %.6g\n",id,time_str(time,3),
            nb,samples,thres);
    for (m=0;m<samples;m++) {
        for (i=0;i<nb;i++) {
            double z=0.0,u1,u2,g;
            for (j=0;j<=i;j++) {
                u1=((double)rand()+1.0)/((double)RAND_MAX+2.0);
                u2=((double)rand()+1.0)/((double)RAND_MAX+2.0);
                g=sqrt(-2.0*log(u1))*cos(2.0*PI*u2);
                z+=L[i+j*nb]*g;
            }
            eta[i]=z;
        }
        info=lambda(nb,2,eta,Q,F,s);
        if (info||s[0]<=0.0||s[1]/s[0]<thres) continue;
        accepted++;
        for (i=0;i<mode_n;i++) {
            int same=1;
            for (j=0;j<nb;j++)
                if (modes[j+i*nb]!=(int)F[j]) {
                    same=0;
                    break;
                }
            if (same) break;
        }
        if (i==mode_n&&mode_n<MC_MAX_MODES) {
            for (j=0;j<nb;j++) modes[j+mode_n*nb]=(int)F[j];
            counts[mode_n]=1;
            mode_n++;
        } else if (i<mode_n) counts[i]++;
        if (mode_n>0) {
            int zero=1;
            for (j=0;j<nb;j++)
                if (F[j]!=0.0) zero=0;
            if (zero) correct++;
        }
    }
    for (i=0;i<mode_n;i++) {
        double dx=0.0,dy=0.0,dz=0.0,norm=0.0;
        if (Qab&&na>0)
            for (j=0;j<na;j++) {
                int k;
                dv[j]=0.0;
                for (k=0;k<nb;k++)
                    for (m=0;m<nb;m++)
                        dv[j]-=Qab[j+k*na]*Qi[k+m*nb]*
                                 modes[m+i*nb];
                if (j==0) dx=dv[j];
                else if (j==1) dy=dv[j];
                else if (j==2) dz=dv[j];
            }
        norm=sqrt(dx*dx+dy*dy+dz*dz);
        fprintf(fp_rtkdbg,"MC_MODE %llu %s %d",id,time_str(time,3),nb);
        for (j=0;j<nb;j++) fprintf(fp_rtkdbg," %d",modes[j+i*nb]);
        fprintf(fp_rtkdbg," %d %.17g %.17g %.17g %.17g %.17g\n",counts[i],
                counts[i]/(double)samples,dx,dy,dz,norm);
    }
    fprintf(fp_rtkdbg,"MC_RESULT %llu %s %d %d %d %d %d %.17g %.17g\n",id,
            time_str(time,3),nb,samples,accepted,correct,
            accepted-correct,accepted/(double)samples,
            accepted?correct/(double)accepted:0.0);
    fflush(fp_rtkdbg);
done:
    free(L);
    free(eta);
    free(F);
    free(Qi);
    free(dv);
    free(modes);
    free(counts);
#else
    (void)time;
    (void)na;
    (void)nb;
    (void)Qab;
    (void)Q;
    (void)samples;
    (void)thres;
#endif
}

#endif
