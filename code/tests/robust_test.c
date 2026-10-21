#include "rtklib.h"
#include <assert.h>

int main(void)
{
    double x[1]={1.0},p[1]={1.0},h[2]={0.0,0.0},v[2]={0.2,0.3};
    double r[4]={1.0,0.2,0.2,1.0},original[4],change;
    int iterations,status;
    matcpy(original,r,2,2);
    status=rtkrobust(x,p,v,h,r,1,2,20,1E-3,.5,&iterations,&change);
    assert(status==0&&iterations==1&&change==0.0);
    assert(memcmp(r,original,sizeof(r))==0);

    /* Exact upper cutoff remains finite; correlated noise keeps its rho. */
    v[0]=8.5;v[1]=.3;
    status=rtkrobust(x,p,v,h,r,1,2,1,1E-3,.5,&iterations,&change);
    assert(status==0&&iterations==1);
    assert(isfinite(r[0])&&fabs(r[0]-1E6)<1E-6);
    assert(fabs(r[1]/sqrt(r[0]*r[3])-.2)<1E-12);
    assert(r[1]==r[2]);

    /* A large innovation with uncertain state needs more than one pass.
     * Iteration must not modify or repeatedly assimilate the fixed prior. */
    h[0]=h[1]=1.0;v[0]=7.0;v[1]=0.0;
    matcpy(r,original,2,2);
    status=rtkrobust(x,p,v,h,r,1,2,20,1E-3,.5,&iterations,&change);
    assert(status>=0&&iterations>1&&iterations<=20);
    assert(status==1||change<=1E-3);
    assert(x[0]==1.0&&p[0]==1.0);
    assert(r[0]>=original[0]&&r[3]>=original[3]);
    assert(fabs(r[1]/sqrt(r[0]*r[3])-.2)<1E-12);
    printf("iteration test: status=%d iterations=%d change=%.9g\n",
           status,iterations,change);
    matcpy(r,original,2,2);
    status=rtkrobust(x,p,v,h,r,1,2,2,1E-3,.5,&iterations,&change);
    assert(status==1&&iterations==2&&change>1E-3&&isfinite(r[0]));

    /* Singular innovation covariance fails atomically. */
    h[0]=h[1]=0.0;
    r[0]=r[1]=r[2]=r[3]=1.0;
    matcpy(original,r,2,2);
    status=rtkrobust(x,p,v,h,r,1,2,20,1E-3,.5,&iterations,&change);
    assert(status<0&&memcmp(r,original,sizeof(r))==0);
    puts("PASS: IGGIII cutoff, covariance correlation, fixed prior and failure handling");
    return 0;
}
