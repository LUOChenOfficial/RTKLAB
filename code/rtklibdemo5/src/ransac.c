/*------------------------------------------------------------------------------
* ransac.c : carrier single-difference consensus selection
*
* notes  : A is column-major np x nv; v is in m and var is in m^2.
*          Sampling, stopping and consensus rules are unchanged.
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include <float.h>

static unsigned int ransacrand(unsigned int *state)
{
    *state^=*state<<13;
    *state^=*state>>17;
    *state^=*state<<5;
    return *state;
}

/* carrier SD RANSAC consensus -----------------------------------------------
* args   : double *A       I   design matrix (np x nv, column-major)
*          double *v       I   SD residuals (nv x 1, m)
*          double *var     I   SD variances (nv x 1, m^2)
*          int    np,nv    I   model dimension and observation count
*          uint   seed     I   deterministic sampling seed
*          uchar  *keep    O   retained observations (nv x 1, 0:reject,1:keep)
* return : consensus size (0: no valid redundant consensus)
* notes  : On failure, a valid keep buffer retains the original observations.
*-----------------------------------------------------------------------------*/
EXPORT int rtkransac(const double *A, const double *v, const double *var,
                     int np, int nv, unsigned int seed, unsigned char *keep)
{
    double As[8*8],vs[8],x[8],Q[8*8],err,emin=DBL_MAX;
    int ix[MAXSAT*NFREQ],i,j,k,trial,valid=0,count,best=0,
        target=RTK_RANSAC_MAX_ITER;
    unsigned char mask[MAXSAT*NFREQ];

    if (!keep||nv<0||nv>MAXSAT*NFREQ) return 0;
    memset(keep,1,nv);
    if (!A||!v||!var||np<4||np>8||nv<=np||
        RTK_RANSAC_MAX_ITER<3||RTK_RANSAC_SIGMA<=0.0||
        RTK_RANSAC_PROBABILITY<=0.0||RTK_RANSAC_PROBABILITY>=1.0)
        return 0;
    for (i=0;i<nv;i++) {
        if (!isfinite(v[i])||!isfinite(var[i])||var[i]<=0.0) return 0;
        for (j=0;j<np;j++)
            if (!isfinite(A[j+i*np])) return 0;
    }
    if (!seed) seed=19938;
    for (trial=0;
         trial<RTK_RANSAC_MAX_ITER&&(trial<target||valid<3);
         trial++) {
        for (i=0;i<nv;i++) ix[i]=i;
        for (i=0;i<np;i++) {
            int sel=i+(int)(ransacrand(&seed)%(nv-i));
            k=ix[i];
            ix[i]=ix[sel];
            ix[sel]=k;
            k=ix[i];
            vs[i]=v[k]/sqrt(var[k]);
            for (j=0;j<np;j++) As[j+i*np]=A[j+k*np]/sqrt(var[k]);
        }
        if (lsq(As,vs,np,np,x,Q)) continue;
        valid++;
        count=0;
        err=0.0;
        for (i=0;i<nv;i++) {
            double r=v[i];
            for (j=0;j<np;j++) r-=A[j+i*np]*x[j];
            mask[i]=isfinite(r)&&fabs(r)<RTK_RANSAC_SIGMA*sqrt(var[i]);
            if (mask[i]) {
                count++;
                err+=r*r/var[i];
            }
        }
        if (count<=np)
            continue; /* no redundant consensus: keep original observations */
        err/=count;
        if (count>best||(count==best&&err<emin)) {
            double fail=1.0-pow(count/(double)nv,np),nreq;
            best=count;
            emin=err;
            memcpy(keep,mask,nv);
            fail=fail<DBL_EPSILON?DBL_EPSILON:fail;
            nreq=ceil(log(1.0-RTK_RANSAC_PROBABILITY)/log(fail));
            target=
                nreq>RTK_RANSAC_MAX_ITER?RTK_RANSAC_MAX_ITER:(int)nreq;
            if (target<3) target=3;
        }
    }
#if RTK_RANSAC_REFIT
    if (best>np) {
        double *Aw=mat(np,nv),*vw=mat(nv,1);
        if (Aw&&vw)
            for (trial=0;trial<5;trial++) {
                int n=0,changed=0;
                for (i=0;i<nv;i++)
                    if (keep[i]) {
                        vw[n]=v[i]/sqrt(var[i]);
                        for (j=0;j<np;j++) Aw[j+n*np]=A[j+i*np]/sqrt(var[i]);
                        n++;
                    }
                if (n<=np||lsq(Aw,vw,np,n,x,Q)) break;
                count=0;
                for (i=0;i<nv;i++) {
                    double r=v[i];
                    for (j=0;j<np;j++) r-=A[j+i*np]*x[j];
                    mask[i]=isfinite(r)&&
                              fabs(r)<RTK_RANSAC_SIGMA*sqrt(var[i]);
                    count+=mask[i];
                    changed|=mask[i]!=keep[i];
                }
                if (count<=np) break; /* preserve last redundant consensus */
                memcpy(keep,mask,nv);
                best=count;
                if (!changed) break;
            }
        free(Aw);
        free(vw);
    }
#endif
    return best;
}
