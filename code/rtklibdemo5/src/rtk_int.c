/*------------------------------------------------------------------------------
* rtk_int.c : RTK integrity lifecycle and solution separation
*
* notes  : DD model capture, FILTER_SS/LS_SS, FDE and monitor output.
*-----------------------------------------------------------------------------*/
#include "rtk_int_internal.h"

#ifdef ENABLE_RTK_INTEGRITY

#define SQR_(x)      ((x)*(x))
#define SQRT_(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define NF_(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP_(opt)     ((opt)->dynamics==0?3:9)
#define NI_(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT_(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL_(opt)     ((opt)->glomodear!=GLO_ARMODE_AUTOCAL?0:NFREQGLO)
#define NR_(opt)     (NP_(opt)+NI_(opt)+NT_(opt)+NL_(opt))
#define IB_(s,f,opt) (NR_(opt)+MAXSAT*(f)+(s)-1)

static FILE *fp_pint=NULL,*fp_pld=NULL,*fp_sub=NULL,*fp_rbias=NULL;
static FILE *fp_gls=NULL;
static FILE *fp_glsmode=NULL;
static int outsub=0,outbias=0;
static char rbiaspath[1024]="";

#if ENABLE_RTK_ARAIM_PL_BIAS_TERM
static int calcbiasenu(rtk_t *rtk, double *be, double *bn, double *bu,
                       int *rows);
#endif

#if ENABLE_RTK_ARAIM_PL_BIAS_TERM
typedef struct {
    int sysidx,band,meas,sat;
    double mu,sigma,bound;
    int count;
} rbiasp_t;
static rbiasp_t *rbiasp=NULL;
static int nrbiasp=0,rbiasload=0,rbiastry=0;
#endif

static void sepenu(const double *base, const double *pos, double *sep)
{
    double bl[3],llh[3],E[9];

    bl[0]=pos[0]-base[0];
    bl[1]=pos[1]-base[1];
    bl[2]=pos[2]-base[2];
    ecef2pos(base,llh);
    xyz2enu(llh,E);
    matmul("NN",3,1,3,1.0,E,bl,0.0,sep);
}

static void qr2cov(const float *qr, double *q)
{
    q[0]=qr[0];
    q[1]=qr[1];
    q[2]=qr[2];
    q[3]=qr[3];
    q[4]=qr[4];
    q[5]=qr[5];
}

static int satinrow(int vflg, int sat)
{
    int su=(vflg>>8)&0xFF;
    int sf=(vflg>>16)&0xFF;

    return sat>0&&(su==sat||sf==sat);
}

static int lsqdd(const rtk_t *rtk, int exsat, double *Qpos, double *sig,
                 double *xenu, int *nrow)
{
    const rtkim_t *mon=&rtk->intg;
    double pos[3],E[9],qtmp[6],*Hpos=NULL,*Rsub=NULL,*W=NULL,
                                  *F=NULL,*N=NULL,*u=NULL,
                                  *vsub=NULL;
    int i,j,m=0,*rows=NULL,ok=0;

    if (nrow) *nrow=0;
    if (xenu) memset(xenu,0,sizeof(double) *3);
    if (sig) memset(sig,0,sizeof(double) *3);
    if (Qpos) memset(Qpos,0,sizeof(double) *9);
    if (!mon->bH||!mon->bR||!mon->bv||mon->bnv<=0||mon->bnx<3) return 0;

    rows=(int *)malloc(sizeof(int)*mon->bnv);
    if (!rows) return 0;
    for (i=0;i<mon->bnv;i++) {
        if (satinrow(mon->bflg[i],exsat)) continue;
        rows[m++]=i;
    }
    if (nrow) *nrow=m;
    if (m<3) {
        free(rows);
        return 0;
    }
    Hpos=mat(3,m);
    Rsub=mat(m,m);
    W=mat(m,m);
    F=mat(3,m);
    N=mat(3,3);
    u=mat(3,1);
    vsub=mat(m,1);
    if (!Hpos||!Rsub||!W||!F||!N||!u||!vsub) {
        free(rows);
        free(Hpos);
        free(Rsub);
        free(W);
        free(F);
        free(N);
        free(u);
        free(vsub);
        return 0;
    }
    for (j=0;j<m;j++) {
        int rj=rows[j];
        vsub[j]=mon->bv[rj];
        for (i=0;i<3;i++) Hpos[i+j*3]=mon->bH[i+rj*mon->bnx];
        for (i=0;i<m;i++) Rsub[i+j*m]=mon->bR[rows[i]+rj*mon->bnv];
    }
    matcpy(W,Rsub,m,m);
    if (matinv(W,m)) {
        goto done;
    }
    matmul("NN",3,m,m,1.0,Hpos,W,0.0,F);
    matmul("NT",3,3,m,1.0,F,Hpos,0.0,N);
    if (matinv(N,3)) {
        goto done;
    }
    for (i=0;i<3;i++) {
        u[i]=0.0;
        for (j=0;j<m;j++) u[i]+=F[i+j*3]*vsub[j];
    }
    if (xenu||Qpos||sig) {
        double xecef[3]={0};
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) xecef[i]+=N[i+j*3]*u[j];
        }
        if (Qpos) {
            for (i=0;i<9;i++) Qpos[i]=N[i];
        }
        qtmp[0]=N[0];
        qtmp[1]=N[4];
        qtmp[2]=N[8];
        qtmp[3]=N[1];
        qtmp[4]=N[5];
        qtmp[5]=N[2];
        if (sig) rtkimsigenu(rtk->sol.rr,qtmp,sig);
        if (xenu) {
            ecef2pos(rtk->sol.rr,pos);
            xyz2enu(pos,E);
            for (i=0;i<3;i++) {
                for (j=0;j<3;j++) xenu[i]+=E[i+j*3]*xecef[j];
            }
        }
    }
    ok=1;
done:
    free(rows);
    free(Hpos);
    free(Rsub);
    free(W);
    free(F);
    free(N);
    free(u);
    free(vsub);
    return ok;
}

static int lssubsets(rtk_t *rtk, double *sig, double *xenu)
{
    rtkim_t *mon=&rtk->intg;
    double msig[3]={0},enu[3]={0};
    int i,k,nrow=0,nvalid=0;

    if (sig) memset(sig,0,sizeof(double) *3);
    if (xenu) memset(xenu,0,sizeof(double) *3);
    mon->nact=0;
    for (i=0;i<mon->ndef;i++) memset(&mon->ss[i],0,sizeof(mon->ss[i]));
    if (!lsqdd(rtk,0,NULL,msig,enu,&nrow)) return 0;
    if (sig)
        for (i=0;i<3;i++) sig[i]=msig[i];
    if (xenu)
        for (i=0;i<3;i++) xenu[i]=enu[i];

    for (i=0;i<mon->ndef;i++) {
        rtkim_def_t *def=mon->def+i;
        rtkim_ss_t *ss=mon->ss+i;
        double ssig[3]={0},subenu[3]={0};
        if (def->mode!=RTKIM_F_SAT||def->sat<=0) continue;
        ss->act=1;
        if (!lsqdd(rtk,def->sat,NULL,ssig,subenu,&nrow)) continue;
        ss->ready=1;
        ss->valid=1;
        ss->qi=rtk->sol.stat;
        ss->ratio=rtk->sol.ratio;
        for (k=0;k<3;k++) {
            ss->sig[k]=ssig[k];
            ss->sep[k]=enu[k]-subenu[k];
        }
        nvalid++;
    }
    mon->nact=nvalid;
    return 1;
}

static void calcplls(rtk_t *rtk)
{
    rtkim_t *mon=&rtk->intg;
    rtkim_pl_t *pl=&mon->pl;
    double msig[3],ssig[3],sep[3];
    double pe,pn,pu,sss,subp,kmain=0.0,ksub=0.0,kss=0.0;
    int i,k;

    memset(pl,0,sizeof(*pl));
    if (rtk->sol.stat!=SOLQ_FIX) {
        rtk->inthpl=0.0;
        rtk->intvpl=0.0;
        return;
    }
    if (!lssubsets(rtk,msig,NULL)) {
        rtk->inthpl=0.0;
        rtk->intvpl=0.0;
        return;
    }
    rtkimthres(&rtk->opt,mon->ndef,&kmain,&ksub,&kss);
    for (i=0;i<3;i++) pl->msig[i]=msig[i];
    pe=kmain*msig[0];
    pn=kmain*msig[1];
    pu=kmain*msig[2];
    pl->hsrc=pl->vsrc=0;

    for (i=0;i<mon->ndef;i++) {
        rtkim_def_t *def=mon->def+i;
        rtkim_ss_t *ss=mon->ss+i;
        if (!ss->act||!ss->valid||def->mode!=RTKIM_F_SAT||
            def->sat<=0)
            continue;
        for (k=0;k<3;k++) {
            ssig[k]=ss->sig[k];
            sep[k]=ss->sep[k];
        }
        for (k=0;k<3;k++) {
            sss=SQRT_(fabs(SQR_(ssig[k])-SQR_(msig[k])));
            subp=ksub*ssig[k]+kss*sss;
            if (k==0&&subp>pe) {
                pe=subp;
                pl->hsrc=def->id;
                pl->hmode=def->mode;
                pl->hsat=def->sat;
                pl->hsig[0]=ssig[0];
                pl->hsig[1]=ssig[1];
                pl->hsig[2]=ssig[2];
                pl->hsss[0]=sss;
                pl->hsss[1]=0.0;
                pl->hsss[2]=0.0;
                pl->hsep[0]=sep[0];
                pl->hsep[1]=sep[1];
                pl->hsep[2]=sep[2];
            }
            if (k==1&&subp>pn) {
                pn=subp;
                pl->hsrc=def->id;
                pl->hmode=def->mode;
                pl->hsat=def->sat;
                pl->hsig[0]=ssig[0];
                pl->hsig[1]=ssig[1];
                pl->hsig[2]=ssig[2];
                pl->hsss[0]=0.0;
                pl->hsss[1]=sss;
                pl->hsss[2]=0.0;
                pl->hsep[0]=sep[0];
                pl->hsep[1]=sep[1];
                pl->hsep[2]=sep[2];
            }
            if (k==2&&subp>pu) {
                pu=subp;
                pl->vsrc=def->id;
                pl->vmode=def->mode;
                pl->vsat=def->sat;
                pl->vsig[0]=ssig[0];
                pl->vsig[1]=ssig[1];
                pl->vsig[2]=ssig[2];
                pl->vsss[0]=0.0;
                pl->vsss[1]=0.0;
                pl->vsss[2]=sss;
                pl->vsep[0]=sep[0];
                pl->vsep[1]=sep[1];
                pl->vsep[2]=sep[2];
            }
        }
    }
    pl->pe=pe;
    pl->pn=pn;
    pl->pu=pu;
    pl->hpl0=SQRT_(SQR_(pe)+SQR_(pn));
    pl->vpl0=pu;
    pl->hpl=pl->hpl0;
    pl->vpl=pl->vpl0;
    pl->hadd=0.0;
    pl->vadd=0.0;
    rtk->inthpl=pl->hpl;
    rtk->intvpl=pl->vpl;
}

static void setpath(char *dst, const char *src, const char *ext)
{
    char *p;

    strcpy(dst,src);
    if ((p=strrchr(dst,'.'))) *p='\0';
    strcat(dst,ext);
}

static void setbasename(char *dst, const char *src, const char *name)
{
    char *p;

    strcpy(dst,src);
    if ((p=strrchr(dst,'\\'))) *(p+1)='\0';
    else if ((p=strrchr(dst,'/')))*(p+1)='\0';
    else dst[0]='\0';
    strcat(dst,name);
}

static void printcenter(FILE *fp, int width, const char *text)
{
    int len,left,right,i;

    if (!fp||width<=0) return;
    if (!text) text="";
    len=(int)strlen(text);
    if (len>=width) {
        fprintf(fp,"%.*s",width,text);
        return;
    }
    left=(width-len)/2;
    right=width-len-left;
    for (i=0;i<left;i++) fputc(' ',fp);
    fputs(text,fp);
    for (i=0;i<right;i++) fputc(' ',fp);
}

static int sysidx(int sys)
{
    if (sys==SYS_GPS) return 0;
    if (sys==SYS_GLO) return 1;
    if (sys==SYS_GAL) return 2;
    if (sys==SYS_CMP) return 3;
    if (sys==SYS_QZS) return 4;
    return -1;
}

static const char *systxt(int sys)
{
    if (sys==SYS_GPS) return "G";
    if (sys==SYS_GLO) return "R";
    if (sys==SYS_GAL) return "E";
    if (sys==SYS_CMP) return "C";
    if (sys==SYS_QZS) return "J";
    return "?";
}

#if ENABLE_RTK_ARAIM_PL_BIAS_TERM
static int satfromid(const char *id)
{
    int sat=satid2no(id);

    if (sat>0) return sat;
    return atoi(id);
}

static int measidx(const char *s)
{
    if (!strcmp(s,"CODE")||!strcmp(s,"code")||!strcmp(s,"P")) return 1;
    if (!strcmp(s,"PHASE")||!strcmp(s,"phase")||!strcmp(s,"L")) return 0;
    return -1;
}

static int splitcsv(char *line, char **cols, int maxcols)
{
    int n=0;
    char *p=line;

    while (n<maxcols) {
        cols[n++]=p;
        if (!(p=strchr(p,','))) break;
        *p++='\0';
    }
    return n;
}

static int loadrbiasfile(const char *path)
{
    FILE *fp=fopen(path,"r");
    char buff[1024],*cols[16];

    rbiasp_t *tmp=NULL;
    int n=0,cap=0;
    if (!fp) return 0;
    while (fgets(buff,sizeof(buff),fp)) {
        rbiasp_t p={0};
        int nc;
        char *q=strchr(buff,'\n');
        if (q) *q='\0';
        nc=splitcsv(buff,cols,16);
        if (nc<9||!strcmp(cols[0],"sysidx")) continue;
        p.sysidx=atoi(cols[0]);
        p.band=atoi(cols[2]);
        p.meas=measidx(cols[3]);
        p.sat=satfromid(cols[4]);
        p.mu=atof(cols[5]);
        p.sigma=atof(cols[6]);
        p.bound=atof(cols[7]);
        p.count=atoi(cols[8]);
        if (p.sysidx<0||p.band<=0||p.meas<0||p.sat<=0) continue;
        if (n>=cap) {
            cap=cap<=0?128:cap*2;
            tmp=(rbiasp_t *)realloc(tmp,sizeof(rbiasp_t)*cap);
            if (!tmp) {
                fclose(fp);
                return 0;
            }
        }
        tmp[n++]=p;
    }
    fclose(fp);
    if (n<=0) {
        free(tmp);
        return 0;
    }
    free(rbiasp);
    rbiasp=tmp;
    nrbiasp=n;
    return 1;
}

static int ensurerbias(void)
{
    const char *paths[]={
        rbiaspath,
        "result\\Open\\rbias_params.csv",
        "results\\bias_envelope\\output\\rbias_params.csv",
        "result\\bias_envelope\\output\\rbias_params.csv",
        "..\\result\\Open\\rbias_params.csv",
        "..\\results\\bias_envelope\\output\\rbias_params.csv",
        "..\\result\\bias_envelope\\output\\rbias_params.csv",
        "..\\..\\result\\Open\\rbias_params.csv",
        "..\\..\\results\\bias_envelope\\output\\rbias_params.csv",
        "..\\..\\result\\bias_envelope\\output\\rbias_params.csv",
        "bias_envelope\\output\\rbias_params.csv"};
    int i;

    if (rbiasload) return 1;
    if (rbiastry) return 0;
    rbiastry=1;
    for (i=0;i<(int)(sizeof(paths)/sizeof(paths[0]));i++) {
        if (!paths[i]||!paths[i][0]) continue;
        if (loadrbiasfile(paths[i])) {
            rbiasload=1;
            return 1;
        }
    }
    return 0;
}

static int findrbias(int sysi, int band, int meas, int sat)
{
    int i;

    for (i=0;i<nrbiasp;i++) {
        if (rbiasp[i].sysidx==sysi&&rbiasp[i].band==band&&
            rbiasp[i].meas==meas&&rbiasp[i].sat==sat)
            return i;
    }
    return -1;
}

static int calcbiasenu(rtk_t *rtk, double *be, double *bn, double *bu,
                       int *rows)
{
    rtkim_t *mon=&rtk->intg;
    double pos[3],E[9],dxyz[3]={0},denu[3],*mu,*R,*W,*tmp,*dx;
    int i,j,k,n=0,nx=rtk->nx,sys,su,sf,frq,meas,idxu,idxf,sysi;
    int sel[RTKIM_MAX_BIAS_ROWS];
    *be=*bn=*bu=0.0;
    if (rows) *rows=0;
    if (!ensurerbias()||mon->bnv<=0||mon->bnx!=nx||!mon->bH||
        !mon->bR)
        return 0;
    for (i=0;i<mon->bnv&&n<RTKIM_MAX_BIAS_ROWS;i++) {
        sf=(mon->bflg[i]>>16)&0xFF;
        su=(mon->bflg[i]>>8)&0xFF;
        meas=(mon->bflg[i]>>4)&0xF;
        frq=(mon->bflg[i]&0xF)+1;
        sys=satsys(su,NULL);
        sysi=sysidx(sys);
        if (sf<=0||su<=0||sysi<0) continue;
        idxu=findrbias(sysi,frq,meas,su);
        idxf=findrbias(sysi,frq,meas,sf);
        if (idxu<0||idxf<0) continue;
        sel[n++]=i;
    }
    if (n<=0) return 0;
    mu=mat(n,1);
    R=mat(n,n);
    W=mat(n,n);
    tmp=mat(nx,1);
    dx=mat(nx,1);
    if (!mu||!R||!W||!tmp||!dx) {
        free(mu);
        free(R);
        free(W);
        free(tmp);
        free(dx);
        return 0;
    }
    for (i=0;i<n;i++) {
        int ri=sel[i];
        sf=(mon->bflg[ri]>>16)&0xFF;
        su=(mon->bflg[ri]>>8)&0xFF;
        meas=(mon->bflg[ri]>>4)&0xF;
        frq=(mon->bflg[ri]&0xF)+1;
        sysi=sysidx(satsys(su,NULL));
        idxu=findrbias(sysi,frq,meas,su);
        idxf=findrbias(sysi,frq,meas,sf);
        mu[i]=rbiasp[idxu].mu-rbiasp[idxf].mu;
        for (j=0;j<n;j++) R[i+j*n]=mon->bR[ri+sel[j]*mon->bnv];
    }
    matcpy(W,R,n,n);
    if (matinv(W,n)) {
        free(mu);
        free(R);
        free(W);
        free(tmp);
        free(dx);
        return 0;
    }
    for (k=0;k<nx;k++) {
        double s=0.0;
        for (i=0;i<n;i++)
            for (j=0;j<n;j++) {
                s+=mon->bH[k+sel[i]*nx]*W[i+j*n]*mu[j];
            }
        tmp[k]=s;
    }
    matmul("NN",nx,1,nx,1.0,rtk->P,tmp,0.0,dx);
    dxyz[0]=dx[0];
    dxyz[1]=dx[1];
    dxyz[2]=dx[2];
    ecef2pos(rtk->sol.rr,pos);
    xyz2enu(pos,E);
    matmul("NN",3,1,3,1.0,E,dxyz,0.0,denu);
    *be=fabs(denu[0]);
    *bn=fabs(denu[1]);
    *bu=fabs(denu[2]);
    if (rows) *rows=n;
    free(mu);
    free(R);
    free(W);
    free(tmp);
    free(dx);
    return 1;
}
#endif

static int hasobs(const obsd_t *obs, int n, int rcv, int sat)
{
    int i;

    for (i=0;i<n;i++)
        if (obs[i].rcv==rcv&&obs[i].sat==sat) return 1;
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
        if (obs[i].rcv==1)(*nu)++;
        else if (obs[i].rcv==2)(*nr)++;
    }
}

static int finddef(rtkim_t *mon, int mode, int sat)
{
    int i;

    for (i=0;i<mon->ndef;i++) {
        if (mon->def[i].mode==mode&&mon->def[i].sat==sat) return i;
    }
    return -1;
}

static int adddef(rtk_t *rtk, int mode, int sat)
{
    rtkim_t *mon=&rtk->intg;
    prcopt_t *opt=&rtk->opt;
    int max=opt->maxintsub;
    int i=mon->ndef;
    if (max<=0) return -1;
    if (i>=max||i>=MAXSAT+1) return -1;
    mon->def[i].id=i+1;
    mon->def[i].mode=mode;
    mon->def[i].sat=sat;
    mon->def[i].fixed=mode==RTKIM_F_NFIX;
    mon->ndef++;
    return i;
}

static int getdef(rtk_t *rtk, int mode, int sat)
{
    int i=finddef(&rtk->intg,mode,sat);

    return i>=0?i:adddef(rtk,mode,sat);
}

static void freechild(rtkim_ss_t *ss)
{
    if (!ss->rtk) return;
    rtkfree(ss->rtk);
    free(ss->rtk);
    ss->rtk=NULL;
}

static void initchild(rtk_t *rtk, rtkim_ss_t *ss)
{
    prcopt_t opt=rtk->opt;
    if (ss->rtk) return;
    opt.intmon=0;
    opt.modear=ARMODE_OFF;
    ss->rtk=(rtk_t *)malloc(sizeof(rtk_t));
    if (!ss->rtk) return;
    rtkinit(ss->rtk,&opt);
    ss->rtk->intchild=1;
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

static void condamb(rtk_t *main, rtkim_def_t *def, rtk_t *child,
                    rtkim_ss_t *ss)
{
    int idx[MAXSAT*NFREQ];
    int i,f,j,k,l,nv=0,nx=child->nx,nf=NF_(&child->opt);

    ss->ambc=0;
    if (main->sol.stat!=SOLQ_FIX||def->mode==RTKIM_F_NFIX) return;
    if (main->na<=0||main->na!=child->na||main->nx!=child->nx) return;

    for (i=1;i<=MAXSAT&&nv<MAXSAT*NFREQ;i++)
        for (f=0;f<nf;f++) {
            j=IB_(i,f,&child->opt);
            if (j<0||j>=nx||j>=main->na) continue;
            if (def->mode==RTKIM_F_SAT&&i==def->sat) continue;
            if (main->xa[j]==0.0||child->x[j]==0.0) continue;
            if (main->Pa[j+j*main->na]<=0.0||
                child->P[j+j*nx]<=0.0)
                continue;
            idx[nv++]=j;
            if (nv>=MAXSAT*NFREQ) break;
        }
    if (nv<=0) return;

    for (i=0;i<nv;i++) {
        double v,r,s;
        double *kvec,*prow;
        j=idx[i];
        v=main->xa[j]-child->x[j];
        r=main->Pa[j+j*main->na]>1E-6?main->Pa[j+j*main->na]
                                              :1E-6;
        s=child->P[j+j*nx]+r;
        if (s<=0.0) continue;
        kvec=(double *)malloc(sizeof(double)*nx);
        prow=(double *)malloc(sizeof(double)*nx);
        if (!kvec||!prow) {
            free(kvec);
            free(prow);
            continue;
        }
        for (k=0;k<nx;k++) kvec[k]=child->P[k+j*nx]/s;
        for (l=0;l<nx;l++) prow[l]=child->P[j+l*nx];
        for (k=0;k<nx;k++) child->x[k]+=kvec[k]*v;
        for (k=0;k<nx;k++)
            for (l=0;l<nx;l++) {
                child->P[k+l*nx]-=kvec[k]*prow[l];
            }
        free(kvec);
        free(prow);
        ss->ambc++;
    }
    if (ss->ambc>0) setsolfromfilter(child);
}

static void runsub(rtk_t *rtk, int idx, const obsd_t *obs, int n,
                   const nav_t *nav, const sta_t *sta)
{
    rtkim_t *mon=&rtk->intg;
    rtkim_def_t *def=mon->def+idx;
    rtkim_ss_t *ss=mon->ss+idx;
    obsd_t ob[MAXOBS*2];
    int m=0,ok,i;

    ss->act=1;
    initchild(rtk,ss);
    if (!ss->rtk) return;
    ss->init=1;
    ss->rtk->opt=rtk->opt;
    ss->rtk->opt.intmon=0;
    ss->rtk->opt.modear=ARMODE_OFF;
    ss->rtk->intchild=1;
#if ENABLE_RTK_RANSAC_ARAIM
    memcpy(ss->rtk->ransacmask,rtk->ransacmask,sizeof(rtk->ransacmask));
#endif
    filtobs(obs,n,def->mode==RTKIM_F_SAT?def->sat:0,ob,&m);
    ok=rtkpos(ss->rtk,ob,m,nav,sta);
    ss->ready=1;
    ss->valid=ok&&ss->rtk->sol.stat!=SOLQ_NONE;
    if (!ss->valid) return;
    condamb(rtk,def,ss->rtk,ss);
    for (i=0;i<3;i++) ss->pos[i]=ss->rtk->sol.rr[i];
    for (i=0;i<6;i++) ss->q[i]=ss->rtk->sol.qr[i];
    rtkimsigenu(ss->pos,ss->q,ss->sig);
    sepenu(rtk->sol.rr,ss->pos,ss->sep);
    ss->qi=ss->rtk->sol.stat;
    ss->ratio=ss->rtk->sol.ratio;
    ss->pdop=ss->rtk->sol.dops[1];
    ss->hdop=ss->rtk->sol.dops[2];
    ss->vdop=ss->rtk->sol.dops[3];
    ss->sig0=ss->rtk->sol.chi2testvalue;
    ss->chisq=ss->rtk->sol.chi2testvalue;
    ss->dof=ss->rtk->sol.dof;
}

static void calcpl(rtk_t *rtk)
{
#if RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS
    if (!rtk->glsassessed) rtkimglsassess(rtk);
    return;
#elif RTK_INT_METHOD==RTK_INT_METHOD_SLOPE
    rtkimslassess(rtk);
    return;
#elif RTK_INT_METHOD==RTK_INT_METHOD_LS_SS
    calcplls(rtk);
    return;
#else
    rtkim_t *mon=&rtk->intg;
    rtkim_pl_t *pl=&mon->pl;
    double msig[3],pe,pn,pu,ssig[3],qmain[6];
    double kmain,ksub,kss;
    double be=0.0,bn=0.0,bu=0.0;
    double hpl0;
    int brows=0;
    int i;

    memset(pl,0,sizeof(*pl));
    if (rtk->sol.stat==SOLQ_NONE) {
        rtk->inthpl=0.0;
        rtk->intvpl=0.0;
        return;
    }
    for (i=0;i<6;i++) qmain[i]=rtk->sol.qr[i];
    rtkimthres(&rtk->opt,mon->ndef,&kmain,&ksub,&kss);
    rtkimsigenu(rtk->sol.rr,qmain,msig);
    for (i=0;i<3;i++) pl->msig[i]=msig[i];
    pe=kmain*msig[0];
    pn=kmain*msig[1];
    pu=kmain*msig[2];
    hpl0=SQRT_(SQR_(pe)+SQR_(pn));
    pl->hsrc=pl->vsrc=0;
    for (i=0;i<mon->ndef;i++) {
        rtkim_ss_t *ss=mon->ss+i;
        rtkim_def_t *def=mon->def+i;
        double sss_e,sss_n,sss_u;
        double sube,subn,subu,hsub;
        if (!ss->act||!ss->valid) continue;
        ssig[0]=ss->sig[0];
        ssig[1]=ss->sig[1];
        ssig[2]=ss->sig[2];
        sss_e=SQRT_(SQR_(ssig[0])-SQR_(msig[0]));
        sss_n=SQRT_(SQR_(ssig[1])-SQR_(msig[1]));
        sss_u=SQRT_(SQR_(ssig[2])-SQR_(msig[2]));
        sube=ksub*ssig[0]+kss*sss_e;
        subn=ksub*ssig[1]+kss*sss_n;
        subu=ksub*ssig[2]+kss*sss_u;
        hsub=SQRT_(SQR_(sube)+SQR_(subn));
        if (hsub>hpl0) {
            hpl0=hsub;
            pe=sube;
            pn=subn;
            pl->hsrc=def->id;
            pl->hmode=def->mode;
            pl->hsat=def->sat;
            pl->hsig[0]=ssig[0];
            pl->hsig[1]=ssig[1];
            pl->hsig[2]=ssig[2];
            pl->hsss[0]=sss_e;
            pl->hsss[1]=sss_n;
            pl->hsss[2]=0.0;
            pl->hsep[0]=ss->sep[0];
            pl->hsep[1]=ss->sep[1];
            pl->hsep[2]=ss->sep[2];
        }
        if (subu>pu) {
            pu=subu;
            pl->vsrc=def->id;
            pl->vmode=def->mode;
            pl->vsat=def->sat;
            pl->vsig[0]=ssig[0];
            pl->vsig[1]=ssig[1];
            pl->vsig[2]=ssig[2];
            pl->vsss[0]=0.0;
            pl->vsss[1]=0.0;
            pl->vsss[2]=sss_u;
            pl->vsep[0]=ss->sep[0];
            pl->vsep[1]=ss->sep[1];
            pl->vsep[2]=ss->sep[2];
        }
    }
    pl->pe=pe;
    pl->pn=pn;
    pl->pu=pu;
    pl->hpl0=hpl0;
    pl->vpl0=pu;
#if ENABLE_RTK_ARAIM_PL_BIAS_TERM
    if (calcbiasenu(rtk,&be,&bn,&bu,&brows)) {
        pl->biasload=1;
        pl->nbias=brows;
        pl->be=be;
        pl->bn=bn;
        pl->bu=bu;
    }
#endif
    pl->hpl=SQRT_(SQR_(pe+pl->be)+SQR_(pn+pl->bn));
    pl->vpl=pu+pl->bu;
    pl->hadd=pl->hpl-pl->hpl0;
    pl->vadd=pl->vpl-pl->vpl0;
    rtk->inthpl=pl->hpl;
    rtk->intvpl=pl->vpl;
#endif
}

static void evalfde(rtk_t *rtk)
{
    rtkim_t *mon=&rtk->intg;
    rtkim_fde_t best={0};
    double kss=0.0;
#if ENABLE_RTK_INTEGRITY_CHI2_GATE
    double chi2_thres=0.0;
#endif
    int i,k;
    if (!rtk->opt.intfde) {
        mon->fde=best;
        return;
    }
#if ENABLE_RTK_INTEGRITY_CHI2_GATE
    if (rtk->sol.dof>0) {
        chi2_thres=
            rtk->sol.dof<=100?chisqr[rtk->sol.dof-1]:chisqr[99];
    }
    if (rtk->sol.stat==SOLQ_NONE||rtk->sol.dof<=0||
        rtk->sol.chi2testvalue<=chi2_thres) {
        mon->fde=best;
        return;
    }
#else
    if (rtk->sol.stat==SOLQ_NONE) {
        mon->fde=best;
        return;
    }
#endif
    rtkimthres(&rtk->opt,mon->ndef,NULL,NULL,&kss);
    for (i=0;i<mon->ndef;i++) {
        rtkim_ss_t *ss=mon->ss+i;
        rtkim_def_t *def=mon->def+i;
        double score=0.0;
        if (!ss->act||!ss->valid) continue;
        for (k=0;k<3;k++) {
            double th=kss*SQRT_(SQR_(ss->sig[k])-SQR_(mon->pl.msig[k]));
            double st=fabs(ss->sep[k]);
            double sc=th>0.0?st/th:0.0;
            if (sc>score) score=sc;
        }
        trace(3,
              "fde score: subset=%d mode=%d sat=%d qi=%d ratio=%.3f "
              "sep=%.4f %.4f %.4f score=%.4f\n",
              def->id,def->mode,def->sat,ss->qi,ss->ratio,ss->sep[0],
              ss->sep[1],ss->sep[2],score);
        if (score>best.score) {
            best.mode=def->mode;
            best.id=def->id;
            best.sat=def->sat;
            best.score=score;
            best.act=def->mode==RTKIM_F_SAT?RTKIM_A_SAT
                       :def->mode==RTKIM_F_NFIX?RTKIM_A_FLOAT
                                                   :RTKIM_A_NONE;
        }
    }
    if (best.score<=1.0) memset(&best,0,sizeof(best));
    mon->fde=best;
    mon->act=best;
    rtk->intfdemode=best.mode;
    rtk->intfdesat=best.sat;
    trace(3,"fde best : subset=%d mode=%d sat=%d score=%.4f act=%d\n",
          best.id,best.mode,best.sat,best.score,best.act);
}

static void subout(const rtk_t *rtk)
{
    const rtkim_t *mon=&rtk->intg;
    int i;
    if (!fp_sub||!rtk->opt.intoutsub) return;
    for (i=0;i<mon->ndef;i++) {
        const rtkim_def_t *def=mon->def+i;
        const rtkim_ss_t *ss=mon->ss+i;
        if (!ss->act) continue;
        if (rtk->opt.intdbgid>0&&rtk->opt.intdbgid!=def->id) continue;
        if (rtk->opt.intdbgsat>0&&rtk->opt.intdbgsat!=def->sat) continue;
        fprintf(fp_sub,
                "%-23s %6d %6d %6d %6d %8d %12.4f %12.4f %12.4f %10.4f %10.4f "
                "%10.4f\n",
                time_str(mon->ep.time,3),def->id,def->mode,def->sat,
                ss->valid,ss->ambc,ss->pos[0],ss->pos[1],ss->pos[2],
                ss->sep[0],ss->sep[1],ss->sep[2]);
    }
}

static double phaserms(const rtk_t *rtk, int *count)
{
    double s=0.0;
    int i,f,n=0,nf=NF_(&rtk->opt);

    for (i=0;i<MAXSAT;i++)
        for (f=0;f<nf&&f<NFREQ;f++) {
            double r=rtk->ssat[i].resc[f];
            if (!rtk->ssat[i].vsat[f]||r==0.0) continue;
            s+=r*r;
            n++;
        }
    if (count) *count=n;
    return n>0?sqrt(s/n):0.0;
}

/* initialize RTK integrity monitor -------------------------------------------
* args   : rtk_t    *rtk   IO  RTK control and monitor state
*          prcopt_t *opt   I   processing options
* return : none
* notes  : Called from rtkinit(); release an existing monitor before reuse.
*-----------------------------------------------------------------------------*/
EXPORT void rtkiminit(rtk_t *rtk, const prcopt_t *opt)
{
    memset(&rtk->intg,0,sizeof(rtk->intg));
    rtk->intg.ena=opt->intmon;
    rtk->intreproc=0;
    rtk->intfdemode=0;
    rtk->intfdesat=0;
    rtk->intfdemode0=0;
    rtk->intfdesat0=0;
    rtk->intfdeact=0;
    rtk->inthpl=rtk->intvpl=0.0;
}

/* free RTK integrity monitor -------------------------------------------------
* args   : rtk_t *rtk      IO  RTK control and monitor state
* return : none
* notes  : Release subset filters and model buffers; clear monitor state.
*-----------------------------------------------------------------------------*/
EXPORT void rtkimfree(rtk_t *rtk)
{
    int i;

    for (i=0;i<rtk->intg.ndef;i++) freechild(rtk->intg.ss+i);
    free(rtk->intg.bH);
    rtk->intg.bH=NULL;
    free(rtk->intg.bR);
    rtk->intg.bR=NULL;
    free(rtk->intg.bv);
    rtk->intg.bv=NULL;
    free(rtk->intg.slK);
    rtk->intg.slK=NULL;
    free(rtk->intg.slQi);
    rtk->intg.slQi=NULL;
    free(rtk->intg.glsQ);
    rtk->intg.glsQ=NULL;
    free(rtk->intg.glsF);
    rtk->intg.glsF=NULL;
    free(rtk->intg.glscandidates);
    memset(&rtk->intg,0,sizeof(rtk->intg));
}

/* save current DD observation model ------------------------------------------
* args   : rtk_t  *rtk     IO  RTK control and monitor state
*          double *H       I   design matrix (nx x nv, column-major)
*          double *R       I   observation covariance (nv x nv, m^2)
*          double *v       I   residuals (nv x 1, m)
*          int    *vflg    I   DD row flags (nv x 1)
*          int    nx,nv    I   state and observation dimensions
* return : none
* notes  : Captures the supplied model for subsequent integrity assessment.
*-----------------------------------------------------------------------------*/
EXPORT void rtkimsaveddr(rtk_t *rtk, const double *H, const double *R,
                         const double *v, const int *vflg, int nx, int nv)
{
    rtkim_t *mon=&rtk->intg;
    int i,n=nv;
    if (!rtk->opt.intmon||rtk->intchild) return;
    if (!H||!R||!v||!vflg||nx<=0||nv<=0) return;
    if (n>RTKIM_MAX_BIAS_ROWS) n=RTKIM_MAX_BIAS_ROWS;
    if (mon->bnx!=nx||!mon->bH) {
        free(mon->bH);
        mon->bH=(double *)malloc(sizeof(double)*nx*RTKIM_MAX_BIAS_ROWS);
    }
    if (!mon->bR)
        mon->bR=(double *)malloc(sizeof(double)*RTKIM_MAX_BIAS_ROWS*
                                   RTKIM_MAX_BIAS_ROWS);
    if (!mon->bv) mon->bv=(double *)malloc(sizeof(double) *RTKIM_MAX_BIAS_ROWS);
    if (!mon->bH||!mon->bR||!mon->bv) return;
    mon->bnx=nx;
    mon->bnv=n;
    for (i=0;i<n;i++) {
        mon->bflg[i]=vflg[i];
        mon->bv[i]=v[i];
    }
    matcpy(mon->bH,H,nx,n);
    for (i=0;i<n*n;i++) mon->bR[i]=0.0;
    for (i=0;i<n;i++) {
        int j;
        for (j=0;j<n;j++) mon->bR[i+j*n]=R[i+j*nv];
    }
}

EXPORT void rtkimoutbias(const rtk_t *rtk, const double *v, const int *vflg,
                         int nv)
{
    double tow;
    int week,i,qi;
    if (!fp_rbias||!outbias||!rtk->opt.intoutbias) return;
    if (!v||!vflg||nv<=0) return;
    tow=time2gpst(rtk->sol.time,&week);
    qi=rtk->sol.stat>=SOLQ_FIX?rtk->sol.stat:SOLQ_FIX;
    for (i=0;i<nv;i++) {
        int sf=(vflg[i]>>16)&0xFF;
        int su=(vflg[i]>>8)&0xFF;
        int meas=(vflg[i]>>4)&0xF;
        int frq=(vflg[i]&0xF)+1;
        int sys=satsys(su,NULL),si=sysidx(sys);
        char idu[16],idf[16];
        if (sf<=0||su<=0||si<0) continue;
        satno2id(su,idu);
        satno2id(sf,idf);
        fprintf(fp_rbias,
                "$RBIAS,%s,%.3f,%d,%d,%.3f,%d,%s,%d,%s,%s,%s,%.6f,%.6f\n",
                time_str(rtk->sol.time,3),tow,rtk->sol.ns,qi,
                rtk->sol.ratio,si,systxt(sys),frq,meas?"CODE":"PHASE",
                idu,idf,v[i],fabs(v[i]));
    }
}

/* prepare fault modes and detect faults --------------------------------------
* args   : rtk_t  *rtk     IO  RTK control and monitor state
*          obsd_t *obs     I   rover and base observations
*          int    nobs     I   number of observations
*          nav_t  *nav     I   navigation data
*          sta_t  *sta     I   station data
* return : none
* notes  : Called after positioning; method selection controls subset/FDE use.
*-----------------------------------------------------------------------------*/
EXPORT void rtkimdetect(rtk_t *rtk, const obsd_t *obs, int nobs,
                        const nav_t *nav, const sta_t *sta)
{
    rtkim_t *mon=&rtk->intg;
    int sat[MAXSAT],ns,i,nu,nr;
    if (!rtk->opt.intmon||rtk->intchild) return;
    if (rtk->intreproc<=0) {
        rtk->intfdemode0=0;
        rtk->intfdesat0=0;
        rtk->intfdeact=0;
    }
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
    mon->ep.fixupd=mon->ep.fixed;
#if RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS
    mon->ndef=mon->nact=0;
    memset(&mon->fde,0,sizeof(mon->fde));
    memset(&mon->act,0,sizeof(mon->act));
    rtk->intfdemode=rtk->intfdesat=0;
    if (rtk->opt.intfde) rtkimglsassess(rtk);
    return;
#elif RTK_INT_METHOD==RTK_INT_METHOD_SLOPE
    /* The slope method supplies PL only; it must not alter AR/FDE state. */
    mon->ndef=0;
    mon->nact=0;
    memset(&mon->fde,0,sizeof(mon->fde));
    memset(&mon->act,0,sizeof(mon->act));
    rtk->intfdemode=0;
    rtk->intfdesat=0;
    return;
#elif RTK_INT_METHOD==RTK_INT_METHOD_LS_SS
    if (rtk->opt.intsat) {
        for (i=0;i<ns;i++) getdef(rtk,RTKIM_F_SAT,sat[i]);
    }
    if (lssubsets(rtk,mon->pl.msig,NULL)) {
        evalfde(rtk);
    } else {
        mon->nact=0;
        memset(&mon->fde,0,sizeof(mon->fde));
    }
    return;
#endif
    if (rtk->opt.intsat) {
        for (i=0;i<ns;i++) getdef(rtk,RTKIM_F_SAT,sat[i]);
    }
    if (mon->ep.fixed&&rtk->opt.intamb) {
        getdef(rtk,RTKIM_F_NFIX,0);
    }
    for (i=0;i<mon->ndef;i++) {
        if (mon->def[i].mode==RTKIM_F_SAT&&
            !hasobs(obs,nobs,1,mon->def[i].sat))
            continue;
        if (mon->def[i].mode==RTKIM_F_NFIX&&!mon->ep.fixed) continue;
        runsub(rtk,i,obs,nobs,nav,sta);
    }
    mon->nact=0;
    for (i=0;i<mon->ndef;i++)
        if (mon->ss[i].act) mon->nact++;
    evalfde(rtk);
}

/* assess protection levels ---------------------------------------------------
* args   : rtk_t *rtk      IO  RTK control and monitor state
* return : none
* notes  : Requires initialized epoch data from rtkimdetect().
*-----------------------------------------------------------------------------*/
EXPORT void rtkimassess(rtk_t *rtk)
{
    if (!rtk->opt.intmon||rtk->intchild) return;
    if (!rtk->intg.init) return;
    calcpl(rtk);
}

EXPORT int rtkimopen(const char *outfile, const prcopt_t *opt)
{
    char path[1024];

    if ((!opt->intmon&&!opt->intoutbias
#if!ENABLE_RTK_DEBUG_OUTPUT
         )
#else
         &&1)
#endif
        ||!outfile||!*outfile)
        return 0;
    setbasename(rbiaspath,outfile,"rbias_params.csv");
    if (opt->intmon) {
        setpath(path,outfile,".pint");
        fp_pint=fopen(path,"w");
        setpath(path,outfile,".pldiag");
        fp_pld=fopen(path,"w");
#if RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS
        setpath(path,outfile,".glsdiag");
        fp_gls=fopen(path,"w");
        if (fp_gls)
            fputs("% time nv expected valid missing alert fixed float_match_m "
                  "subset_seconds HPL VPL Hmode Hsat1 Hsat2 Vmode Vsat1 Vsat2 "
                  "geometry_failed solve_failed covariance_failed obs_alert "
                  "fixed_alert obs_ratio fixed_ratio candidate_HPL "
                  "candidate_VPL\n",
                  fp_gls);
        if (opt->intoutsub) {
            setpath(path,outfile,".glsmodes");
            fp_glsmode=fopen(path,"w");
            if (fp_glsmode)
                fputs("% time mode sat1 sat2 sepE sepN sepU sigE sigN sigU\n",
                      fp_glsmode);
        }
#endif
        outsub=opt->intoutsub;
        if (outsub) {
            setpath(path,outfile,".subdiag");
            fp_sub=fopen(path,"w");
        }
    } else {
        outsub=0;
    }
    outbias=opt->intoutbias;
    if (outbias) {
        setpath(path,outfile,".rbias");
        fp_rbias=fopen(path,"w");
    }
    if (fp_pint) {
        fputs("% ",fp_pint);
        /* "% " is a two-character comment prefix, so use 21 to match the
         * 23-character timestamp field in the data rows. */
        printcenter(fp_pint,21,"time");
        fputc(' ',fp_pint);
        printcenter(fp_pint,12,"x");
        fputc(' ',fp_pint);
        printcenter(fp_pint,12,"y");
        fputc(' ',fp_pint);
        printcenter(fp_pint,12,"z");
        fputc(' ',fp_pint);
        printcenter(fp_pint,4,"nsat");
        fputc(' ',fp_pint);
        printcenter(fp_pint,2,"QI");
        fputc(' ',fp_pint);
        printcenter(fp_pint,7,"ratio");
        fputc(' ',fp_pint);
        printcenter(fp_pint,9,"chi2");
        fputc(' ',fp_pint);
        printcenter(fp_pint,6,"subT");
        fputc(' ',fp_pint);
        printcenter(fp_pint,6,"subA");
        fputc(' ',fp_pint);
        printcenter(fp_pint,8,"HPL");
        fputc(' ',fp_pint);
        printcenter(fp_pint,8,"VPL");
        fputc(' ',fp_pint);
        printcenter(fp_pint,7,"FDEmode");
        fputc(' ',fp_pint);
        printcenter(fp_pint,6,"FDEsat");
        fputc(' ',fp_pint);
        printcenter(fp_pint,13,"phase_res_rms");
        fputc(' ',fp_pint);
        printcenter(fp_pint,6,"nPhase");
        fputc(' ',fp_pint);
        printcenter(fp_pint,10,"FDEpreMode");
        fputc(' ',fp_pint);
        printcenter(fp_pint,9,"FDEpreSat");
        fputc(' ',fp_pint);
        printcenter(fp_pint,9,"FDEact");
        fputc('\n',fp_pint);
    }
    if (fp_pld) {
        fputs("% ",fp_pld);
        printcenter(fp_pld,23,"time");
        fputc(' ',fp_pld);
        printcenter(fp_pld,4,"nsat");
        fputc(' ',fp_pld);
        printcenter(fp_pld,2,"QI");
        fputc(' ',fp_pld);
        printcenter(fp_pld,7,"ratio");
        fputc(' ',fp_pld);
        printcenter(fp_pld,6,"subT");
        fputc(' ',fp_pld);
        printcenter(fp_pld,6,"subA");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HPL");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VPL");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"bias_e");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"bias_n");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"bias_u");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HPL0");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VPL0");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HADD");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VADD");
        fputc(' ',fp_pld);
        printcenter(fp_pld,5,"bLd");
        fputc(' ',fp_pld);
        printcenter(fp_pld,5,"bRow");
        fputc(' ',fp_pld);
        printcenter(fp_pld,4,"Hsrc");
        fputc(' ',fp_pld);
        printcenter(fp_pld,5,"Hmode");
        fputc(' ',fp_pld);
        printcenter(fp_pld,4,"Hsat");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"MsigE");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"MsigN");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"MsigU");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HsigE");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HsigN");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HsigU");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HsssE");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HsssN");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HsssU");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HsepE");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HsepN");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"HsepU");
        fputc(' ',fp_pld);
        printcenter(fp_pld,4,"Vsrc");
        fputc(' ',fp_pld);
        printcenter(fp_pld,5,"Vmode");
        fputc(' ',fp_pld);
        printcenter(fp_pld,4,"Vsat");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VsigE");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VsigN");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VsigU");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VsssE");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VsssN");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VsssU");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VsepE");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VsepN");
        fputc(' ',fp_pld);
        printcenter(fp_pld,8,"VsepU");
        fputc(' ',fp_pld);
        printcenter(fp_pld,10,"FDEpreMode");
        fputc(' ',fp_pld);
        printcenter(fp_pld,9,"FDEpreSat");
        fputc(' ',fp_pld);
        printcenter(fp_pld,9,"FDEact");
        fputc('\n',fp_pld);
    }
    if (fp_sub) {
        fputs("% ",fp_sub);
        printcenter(fp_sub,23,"time");
        fputc(' ',fp_sub);
        printcenter(fp_sub,6,"subset");
        fputc(' ',fp_sub);
        printcenter(fp_sub,6,"mode");
        fputc(' ',fp_sub);
        printcenter(fp_sub,6,"sat");
        fputc(' ',fp_sub);
        printcenter(fp_sub,6,"valid");
        fputc(' ',fp_sub);
        printcenter(fp_sub,8,"amb_cond");
        fputc(' ',fp_sub);
        printcenter(fp_sub,12,"x");
        fputc(' ',fp_sub);
        printcenter(fp_sub,12,"y");
        fputc(' ',fp_sub);
        printcenter(fp_sub,12,"z");
        fputc(' ',fp_sub);
        printcenter(fp_sub,10,"sep_e");
        fputc(' ',fp_sub);
        printcenter(fp_sub,10,"sep_n");
        fputc(' ',fp_sub);
        printcenter(fp_sub,10,"sep_u");
        fputc('\n',fp_sub);
    }
    if (fp_rbias)
        fprintf(fp_rbias,"%% "
                          "$RBIAS,time,tow_s,nsat,qi,ratio,sysidx,sys,band,"
                          "meas,sat_u,sat_f,resid_m,absres_m\n");
    return fp_pint||fp_pld||fp_rbias
#if ENABLE_RTK_DEBUG_OUTPUT
           ||0
#endif
        ;
}

EXPORT void rtkimclose(void)
{
    if (fp_pint) fclose(fp_pint);
    fp_pint=NULL;
    if (fp_pld) fclose(fp_pld);
    fp_pld=NULL;
    if (fp_sub) fclose(fp_sub);
    fp_sub=NULL;
    if (fp_gls) fclose(fp_gls);
    fp_gls=NULL;
    if (fp_glsmode) fclose(fp_glsmode);
    fp_glsmode=NULL;
    if (fp_rbias) fclose(fp_rbias);
    fp_rbias=NULL;
    rbiaspath[0]='\0';
}

EXPORT void rtkimglsout(gtime_t time, int mode, int sat1, int sat2,
                        const double *sep, const double *sigma)
{
    if (fp_glsmode)
        fprintf(fp_glsmode,
                "%s %d %d %d %.17g %.17g %.17g %.17g %.17g %.17g\n",
                time_str(time,3),mode,sat1,sat2,sep[0],sep[1],sep[2],
                sigma[0],sigma[1],sigma[2]);
}

EXPORT void rtkimout(const rtk_t *rtk)
{
    const rtkim_t *mon=&rtk->intg;
    const rtkim_pl_t *pl=&mon->pl;
    int ndef=RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS?mon->glsexp:mon->ndef;
    int nact=RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS?mon->glsvalid:mon->nact;
    if (!rtk->opt.intmon||!mon->init) return;
    if (fp_pint) {
        int phase_count=0;
        double phase_r=phaserms(rtk,&phase_count);
        fprintf(fp_pint,
                "%-23s %12.4f %12.4f %12.4f %4d %2d %7.3f %9.4f %6d %6d %8.4f "
                "%8.4f %7d %6d %13.6f %6d %10d %9d %9d\n",
                time_str(mon->ep.time,3),rtk->sol.rr[0],rtk->sol.rr[1],
                rtk->sol.rr[2],rtk->sol.ns,rtk->sol.stat,rtk->sol.ratio,
                rtk->sol.chi2testvalue,ndef,nact,pl->hpl,pl->vpl,
                mon->fde.mode,mon->fde.sat,phase_r,phase_count,
                rtk->intfdemode0,rtk->intfdesat0,rtk->intfdeact);
    }
    if (fp_pld) {
        fprintf(fp_pld,
                "%-23s %4d %2d %7.3f %6d %6d %8.4f %8.4f %8.4f %8.4f %8.4f "
                "%8.4f %8.4f %8.4f %8.4f %5d %5d %4d %5d %4d %8.4f %8.4f "
                "%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f "
                "%4d %5d %4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f "
                "%8.4f %10d %9d %9d\n",
                time_str(mon->ep.time,3),rtk->sol.ns,rtk->sol.stat,
                rtk->sol.ratio,ndef,nact,pl->hpl,pl->vpl,pl->be,pl->bn,
                pl->bu,pl->hpl0,pl->vpl0,pl->hadd,pl->vadd,pl->biasload,
                pl->nbias,pl->hsrc,pl->hmode,pl->hsat,pl->msig[0],
                pl->msig[1],pl->msig[2],pl->hsig[0],pl->hsig[1],
                pl->hsig[2],pl->hsss[0],pl->hsss[1],pl->hsss[2],
                pl->hsep[0],pl->hsep[1],pl->hsep[2],pl->vsrc,pl->vmode,
                pl->vsat,pl->vsig[0],pl->vsig[1],pl->vsig[2],pl->vsss[0],
                pl->vsss[1],pl->vsss[2],pl->vsep[0],pl->vsep[1],
                pl->vsep[2],rtk->intfdemode0,rtk->intfdesat0,
                rtk->intfdeact);
    }
    subout(rtk);
    if (fp_gls)
        fprintf(fp_gls,
                "%s %d %d %d %d %d %d %.9g %.6f %.9g %.9g %d %d %d %d %d %d "
                "%d %d %d %d %d %.9g %.9g %.9g %.9g\n",
                time_str(mon->ep.time,3),mon->bnv,mon->glsexp,
                mon->glsvalid,mon->glsmiss,mon->glsalert,mon->ep.fixed,
                mon->glsmatch,mon->glstime,pl->hpl,pl->vpl,pl->hmode,
                pl->hsat,mon->glshsat2,pl->vmode,pl->vsat,mon->glsvsat2,
                mon->glsgeofail,mon->glssolfail,mon->glscovfail,
                mon->glsobsalert,mon->glsfixalert,mon->glsobsratio,
                mon->glsfixratio,pl->hpl0,pl->vpl0);
}

#endif
