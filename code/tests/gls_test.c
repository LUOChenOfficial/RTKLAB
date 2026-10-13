#include "rtklib.h"

static int flag(int ref, int target) {return (ref<<16)|(target<<8)|(1<<4);}

int main(void)
{
    const int nx=5,nv=5;
    double x[5]={-2186035.8,5181367.8,2999261.2,2.0,0.0};
    double P[25]={0},H[25]={0},R[25]={0},v[5]={.1,-.2,.15,.05,-.1};
    double directions[15]={1,0,0,0,1,0,0,0,1,1,1,0,0,1,1};
    int flags[5],rows[5],ixref[5],nr,i,j,k,a,b,checks=0;
    rtk_t *rtk=(rtk_t *)calloc(1,sizeof(rtk_t));
    if (!rtk) return 2;
    rtk->nx=nx;rtk->na=3;rtk->opt=prcopt_default;
    rtk->opt.niter=1;rtk->opt.intmon=1;
    rtk->opt.modear=ARMODE_CONT;rtk->opt.glomodear=GLO_ARMODE_OFF;
    rtk->opt.intsat=1;rtk->opt.intamb=1;
    rtk->numnv[0]=2;rtk->numnv[1]=3;
    for (i=0;i<nx;i++) P[i+i*nx]=i<3?4.0:2.0;
    P[0+3*nx]=P[3]=.4;
    for (i=0;i<nv;i++) {
        flags[i]=flag(1,i+2);
        for (j=0;j<3;j++) H[j+i*nx]=directions[j+i*3];
        H[3+i*nx]=.2*(i+1);H[4+i*nx]=1.0; /* dormant state must be ignored */
        for (j=0;j<nv;j++) R[i+j*nv]=.1+(i==j?.2:0.0);
    }
    rtkimsaveddr(rtk,H,R,v,flags,nx,nv);
    rtkimglsprior(rtk,x,P,H,R,nx,nv);
    if (!rtk->intg.glsready) return 3;
    for (a=0;a<=6;a++) for (b=a;b<=6;b++) {
        double pos[3],cov[9],xs[5],ps[25],hs[25],rs[25],vs[5];
        if (a==0 && b!=0) continue;
        nr=rtkglsddmat(flags,nv,a,b,rows,ixref);
        if (nr==0) continue;
        matcpy(xs,x,nx,1);matcpy(ps,P,nx,nx);
        for (i=0;i<nr;i++) {
            int r=rows[i],anchor=ixref[i];
            vs[i]=v[r]-(anchor>=0?v[anchor]:0.0);
            for (k=0;k<nx;k++) hs[k+i*nx]=H[k+r*nx]-(anchor>=0?H[k+anchor*nx]:0.0);
            for (j=0;j<nr;j++) {
                int s=rows[j],other=ixref[j];
                rs[i+j*nr]=R[r+s*nv]-(anchor>=0?R[anchor+s*nv]:0.0)
                    -(other>=0?R[r+other*nv]:0.0)+(anchor>=0&&other>=0?R[anchor+other*nv]:0.0);
            }
            if (a==1 && rows[i]+2==a) return 4;
        }
        if (filter(xs,ps,hs,vs,rs,nx,nr,rtk)) return 5;
        if (!rtkglsupdate(rtk->intg.glsQ,rtk->intg.glsF,v,x,rtk->intg.glscov,nv,rows,ixref,nr,pos,cov)) return 6;
        for (i=0;i<3;i++) {
            if (fabs(xs[i]-pos[i])>1e-8) return 7;
            for (j=0;j<3;j++) if (fabs(ps[i+j*nx]-cov[i+j*3])>1e-8) return 8;
        }
        checks++;
    }
    /* Explicit reference cancellation: all retained rows must use anchor C02. */
    nr=rtkglsddmat(flags,nv,1,3,rows,ixref);
    if (nr!=3) return 9;
    for (i=0;i<nr;i++) if (ixref[i]!=0 || rows[i]==1) return 10;
    /* Fixed output deliberately differs by 100 m: actual separation enters PL. */
    {
        double pos[3],cov[9],Pa[9]={0};
        nr=rtkglsddmat(flags,nv,0,0,rows,ixref);
        if (!rtkglsupdate(rtk->intg.glsQ,rtk->intg.glsF,v,x,rtk->intg.glscov,nv,rows,ixref,nr,pos,cov)) return 11;
        rtk->x=pos;rtk->Pa=Pa;rtk->sol.stat=SOLQ_FIX;
        for (i=0;i<3;i++) {rtk->sol.rr[i]=pos[i];Pa[i+i*3]=1e-5;}
        rtk->sol.rr[0]+=100.0;
        rtkimglsassess(rtk);
        if (rtk->intg.glsexp!=43 || !rtk->intg.glsalert || !isinf(rtk->inthpl)) return 12;
        rtk->opt.intamb=0;
        rtkimglsassess(rtk);
        if (!isinf(rtk->inthpl)) return 13;
        rtk->intg.glsready=0;
        rtkimglsassess(rtk);
        if (!isinf(rtk->intvpl)) return 14;
    }
    rtkimfree(rtk);free(rtk);
    printf("PASS: %d GLS/KF subset comparisons, DD reference cancellation, wrong-fixed alarm, missing-model unavailable\n",checks);
    return 0;
}
