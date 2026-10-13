/*------------------------------------------------------------------------------
* rtk_int_common.c : common RTK integrity utilities
*
* notes  : Probability thresholds and ECEF-to-ENU standard deviations.
*-----------------------------------------------------------------------------*/
#include "rtk_int_internal.h"

#ifdef ENABLE_RTK_INTEGRITY
#define SQRT_(x) ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

double rtkimprob(double p, double def)
{
    return p>0.0&&p<1.0?p:def;
}

double rtkimnorminv(double p)
{
    static const double a[]={
        -3.969683028665376e+01, 2.209460984245205e+02,
        -2.759285104469687e+02, 1.383577518672690e+02,
        -3.066479806614716e+01, 2.506628277459239e+00
    };
    static const double b[]={
        -5.447609879822406e+01, 1.615858368580409e+02,
        -1.556989798598866e+02, 6.680131188771972e+01,
        -1.328068155288572e+01
    };
    static const double c[]={
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
        4.374664141464968e+00, 2.938163982698783e+00
    };
    static const double d[]={
        7.784695709041462e-03, 3.224671290700398e-01,
        2.445134137142996e+00, 3.754408661907416e+00
    };
    double q,r;

    if (p<=0.0) return -1E9;
    if (p>=1.0) return 1E9;
    if (p<0.02425) {
        q=sqrt(-2.0*log(p));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
    if (p>1.0-0.02425) {
        q=sqrt(-2.0*log(1.0-p));
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
    q=p-0.5;
    r=q*q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/
           ((((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r)+1.0);
}

void rtkimthres(const prcopt_t *opt, int nset, double *kmain, double *ksub,
                double *kss)
{
    double pfa=rtkimprob(opt?opt->intpfa:0.0,RTKINT_DEFAULT_PFA);
    double pmd=rtkimprob(opt?opt->intpmd:0.0,RTKINT_DEFAULT_PMD);
    double kfa=rtkimnorminv(1.0-pfa*0.5);
    double kmd=rtkimnorminv(1.0-pmd*0.5);
    double pss;

    if (nset<=0) nset=1;
    pss=rtkimprob(pfa/(4.0*nset),RTKINT_DEFAULT_PFA/4.0);
    if (kmain) *kmain=kfa;
    if (ksub)  *ksub =kmd;
    if (kss)   *kss  =rtkimnorminv(1.0-pss);
}

void rtkimsigenu(const double *rr, const double *qr, double *sig)
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

int rtkimsat(int vflg)
{
    int su=(vflg>>8)&0xFF;
    int sf=(vflg>>16)&0xFF;

    return su>0?su:sf;
}


#endif
