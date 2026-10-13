#include "rtklib.h"
#include <assert.h>

int main(void)
{
    double tau,lambda,ep[6]={2025,4,25,7,47,0};
    /* Independent 2-D normal quadrature reference for noncentral chi-square. */
    assert(rtkslopethres(2,0.05,0.2,&tau,&lambda));
    assert(fabs(tau-5.991464547107982)<1E-10);
    assert(fabs(lambda-9.63468897578623)<1E-6);
    assert(!rtkslopethres(0,0.05,0.2,&tau,&lambda));
    assert(!rtkslopethres(2,0.0,0.2,&tau,&lambda));
#if ENABLE_RTK_RANSAC_SLOPE
    {
        rtk_t *rtk=(rtk_t *)calloc(1,sizeof(rtk_t));
        prcopt_t opt=prcopt_default;
        double x[4]={1,2,3,4},P[16]={4,0,0,1.2, 0,3,0,0, 0,0,2,0, 1.2,0,0,5};
        double H[8]={1,0,0,1, 0,1,0,1},R[4]={0.5,0.2,0.2,0.7},v[2]={0,0};
        int flags[2]={(1<<16)|(3<<8),(1<<16)|(4<<8)};
        double det=62.57;
        assert(rtk);
        opt.intmon=1;
        opt.modear=ARMODE_OFF; opt.niter=1;
        rtkiminit(rtk,&opt); rtk->opt=opt; rtk->nx=4;
        rtkimsaveddr(rtk,H,R,v,flags,4,2);
        rtkimslprior(rtk,x,P,H,R,4,2);
        assert(rtk->intg.slready);
        assert(fabs(rtk->intg.slQi[0]-8.7/det)<1E-12);
        assert(fabs(rtk->intg.slQi[1]+6.4/det)<1E-12);
        assert(fabs(rtk->intg.slK[0]-37.56/det)<1E-12);
        assert(fabs(rtk->intg.slK[3]+19.0/det)<1E-12);
        /* Actual estimator response to a reference SD fault c=[1,1]. */
        assert(fabs(rtk->intg.slK[0]+rtk->intg.slK[3]-18.56/det)<1E-12);
        x[3]=0.0; /* filter() excludes zero-valued states */
        rtkimslprior(rtk,x,P,H,R,4,2);
        assert(fabs(rtk->intg.slQi[0]-3.7/(4.5*3.7-0.04))<1E-12);
        /* Detected innovation / missing snapshot must produce unavailable PL. */
        rtk->sol.time=epoch2time(ep); rtk->sol.stat=SOLQ_FLOAT;
        rtk->x=x; rtk->P=P; rtk->intg.init=1;
        rtk->intg.bv[0]=1000;
        rtkimassess(rtk);
        assert(isinf(rtk->inthpl) && isinf(rtk->intvpl));
        rtk->intg.slready=0;
        rtkimassess(rtk);
        assert(isinf(rtk->inthpl) && isinf(rtk->intvpl));
        /* Inject physical SD faults through the actual filter, including the
         * shared reference. Compare their position shifts with reported PL. */
        {
            double origin[4]={-2396539.8,5350742.8,2502969.9,4};
            double nominal[4],faulted[4],post[16],noise[4],vin[2],pos[3],E[9],dx[3],de[3];
            double directions[6]={1,1,-1,0,0,-1},d,b,pmode;
            int a,i,j,sign;
            rtk->opt.intpfa=1E-4;
            rtk->opt.intpmd=1E-3;
            rtk->numnv[0]=2;
            rtkimslprior(rtk,origin,P,H,R,4,2);
            pmode=1E-3/4;
            assert(rtkslopethres(2,1E-4,pmode,&tau,&lambda));
            memcpy(nominal,origin,sizeof(nominal)); memcpy(post,P,sizeof(post));
            memcpy(noise,R,sizeof(noise)); vin[0]=vin[1]=0;
            assert(!filter(nominal,post,H,vin,noise,4,2,rtk));
            ecef2pos(nominal,pos); xyz2enu(pos,E);
            for (a=0;a<3;a++) for (sign=-1;sign<=1;sign+=2) {
                d=0;
                for (i=0;i<2;i++) for (j=0;j<2;j++)
                    d+=directions[i+a*2]*rtk->intg.slQi[i+j*2]*directions[j+a*2];
                b=sign*0.5*sqrt(lambda/d);
                for (i=0;i<2;i++) vin[i]=b*directions[i+a*2];
                memcpy(faulted,origin,sizeof(faulted)); memcpy(post,P,sizeof(post));
                memcpy(noise,R,sizeof(noise));
                assert(!filter(faulted,post,H,vin,noise,4,2,rtk));
                for (i=0;i<2;i++) rtk->intg.bv[i]=vin[i];
                rtk->x=faulted; rtk->P=post;
                rtkimassess(rtk);
                assert(isfinite(rtk->inthpl) && isfinite(rtk->intvpl));
                for (i=0;i<3;i++) dx[i]=faulted[i]-nominal[i];
                matmul("NN",3,1,3,1.0,E,dx,0.0,de);
                assert(hypot(de[0],de[1])<=rtk->inthpl);
                assert(fabs(de[2])<=rtk->intvpl);
            }
        }
        rtkimfree(rtk); free(rtk);
    }
#endif
    puts("PASS: slope thresholds, full-state gain, active states, unavailable PL, reference/target fault envelopes");
    return 0;
}
