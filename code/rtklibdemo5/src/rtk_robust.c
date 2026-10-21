/*------------------------------------------------------------------------------
* rtk_robust.c : iterated IGGIII observation covariance reweighting
*
* notes  : The state prior and linearized innovations are held fixed. Reweight
*          the original R, never assimilate an observation more than once.
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* IGGIII equivalent observation covariance ----------------------------------
* args   : double *x,P,v,H I   fixed prior, innovations and nx x nv design H
*          double *R      IO  original/equivalent observation covariance
*          int nx,nv      I   state and observation dimensions
*          int maxiter    I   weight iteration limit (1: single-pass IGGIII)
*          double tol     I   maximum absolute log-weight change tolerance
*          double relaxation I log-weight update fraction after first pass
*          int *iterations O  completed weight updates
*          double *change O   last maximum absolute log-weight change
* return : 0:converged/single pass, 1:iteration limit with valid last weights,
*          <0:numerical/input failure, R remains unchanged
* notes  : Column-major matrices. Standardization uses the initial unweighted
*          residual covariance throughout, so inflating R cannot inflate the
*          reference noise scale. No persistent state/covariance is updated.
*-----------------------------------------------------------------------------*/
EXPORT int rtkrobust(const double *x, const double *P, const double *v,
                     const double *H, double *R, int nx, int nv, int maxiter,
                     double tol, double relaxation, int *iterations,
                     double *change)
{
    double *p=NULL,*h=NULL,*hp=NULL,*c=NULL,*s=NULL,*r0=NULL,*rw=NULL;
    double *rq=NULL,*q=NULL,*res=NULL,*scale=NULL,*weight=NULL;
    double delta=0.0,k0=2.0,k1=8.5,minweight=1E-6;
    int *ix=NULL,i,j,n=0,iter,status=-1;
    if (iterations) *iterations=0;
    if (change) *change=0.0;
    if (!x||!P||!v||!H||!R||nx<=0||nv<=0||maxiter<1||
        !isfinite(tol)||tol<=0.0||!isfinite(relaxation)||
        relaxation<=0.0||relaxation>1.0) return -1;
    for (i=0;i<nv;i++) {
        if (!isfinite(v[i])) return -1;
        for (j=0;j<nv;j++) if (!isfinite(R[i+j*nv])) return -1;
    }
    ix=imat(nx,1);
    if (!ix) goto done;
    for (i=0;i<nx;i++) if (x[i]!=0.0&&P[i+i*nx]>0.0) ix[n++]=i;
    if (!n) goto done;
    p=mat(n,n); h=mat(n,nv); hp=mat(nv,n); c=mat(nv,nv);
    s=mat(nv,nv); r0=mat(nv,nv); rw=mat(nv,nv); rq=mat(nv,nv);
    q=mat(nv,nv); res=mat(nv,1); scale=mat(nv,1); weight=mat(nv,1);
    if (!p||!h||!hp||!c||!s||!r0||!rw||!rq||!q||!res||!scale||!weight)
        goto done;
    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) {
            p[i+j*n]=P[ix[i]+ix[j]*nx];
            if (!isfinite(p[i+j*n])) goto done;
        }
        for (j=0;j<nv;j++) {
            h[i+j*n]=H[ix[i]+j*nx];
            if (!isfinite(h[i+j*n])) goto done;
        }
    }
    matcpy(r0,R,nv,nv); matcpy(rw,R,nv,nv);
    matmul("TN",nv,n,n,1.0,h,p,0.0,hp);
    matmul("NN",nv,nv,n,1.0,hp,h,0.0,c);
    for (i=0;i<nv*nv;i++) s[i]=c[i]+r0[i];
    if (matinv(s,nv)) goto done;
    matmul("NN",nv,nv,nv,1.0,r0,s,0.0,rq);
    matmul("NT",nv,nv,nv,1.0,rq,r0,0.0,q);
    for (i=0;i<nv;i++) {
        if (!isfinite(q[i+i*nv])||q[i+i*nv]<=0.0) goto done;
        scale[i]=sqrt(q[i+i*nv]); weight[i]=1.0;
    }
    for (iter=0;iter<maxiter;iter++) {
        if (iter) {
            for (i=0;i<nv*nv;i++) s[i]=c[i]+rw[i];
            if (matinv(s,nv)) goto done;
            matmul("NN",nv,nv,nv,1.0,rw,s,0.0,rq);
        }
        matmul("NN",nv,1,nv,1.0,rq,v,0.0,res);
        delta=0.0;
        for (i=0;i<nv;i++) {
            double t=fabs(res[i])/scale[i],target,next,d;
            if (!isfinite(t)) goto done;
            if (t<=k0) target=1.0;
            else if (t>=k1) target=minweight;
            else target=k0/t*pow((k1-t)/(k1-k0),2.0);
            if (target<minweight) target=minweight;
            next=iter?exp((1.0-relaxation)*log(weight[i])+
                         relaxation*log(target)):target;
            d=fabs(log(next)-log(weight[i]));
            if (d>delta) delta=d;
            weight[i]=next;
        }
        for (j=0;j<nv;j++) for (i=0;i<nv;i++) {
            rw[i+j*nv]=r0[i+j*nv]/sqrt(weight[i]*weight[j]);
            if (!isfinite(rw[i+j*nv])) goto done;
        }
        if (iterations) *iterations=iter+1;
        if (change) *change=delta;
        if (delta<=tol||maxiter==1) { status=0; break; }
    }
    if (status<0) status=1;
    matcpy(R,rw,nv,nv);
done:
    free(ix); free(p); free(h); free(hp); free(c); free(s); free(r0);
    free(rw); free(rq); free(q); free(res); free(scale); free(weight);
    return status;
}
