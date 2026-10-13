#include "rtklib.h"
#include <assert.h>

int main(void)
{
    double A[5*16],v[16],var[16],x[5]={0.2,-0.4,0.1,2.0,-3.0};
    unsigned char keep[16],repeat[16];
    int i,j,n;
    for (i=0;i<16;i++) {
        double az=i*2.399963229728653,el=0.2+0.07*i;
        A[0+i*5]=cos(el)*cos(az); A[1+i*5]=cos(el)*sin(az);
        A[2+i*5]=sin(el); A[3+i*5]=i<8?0.19:0.0;
        A[4+i*5]=i>=8?0.24:0.0;
        v[i]=0.0; var[i]=0.0001;
        for (j=0;j<5;j++) v[i]+=A[j+i*5]*x[j];
    }
    assert(rtkransac(A,v,var,5,16,19938,keep)==16);
    v[0]+=3.0; v[12]-=2.0;
    assert(rtkransac(A,v,var,5,16,19938,keep)==14);
    for (i=0;i<16;i++) assert(keep[i]==(i!=0 && i!=12));
    assert(rtkransac(A,v,var,5,16,19938,repeat)==14);
    assert(!memcmp(keep,repeat,sizeof(keep)));
    assert(rtkransac(A,v,var,5,5,19938,keep)==0);
    for (i=0;i<5;i++) assert(keep[i]);
    memset(A,0,sizeof(A));
    assert(rtkransac(A,v,var,5,16,19938,keep)==0);
    for (i=0;i<16;i++) assert(keep[i]);
    var[3]=0.0;
    n=rtkransac(A,v,var,5,16,19938,keep);
    assert(n==0);
    for (i=0;i<16;i++) assert(keep[i]);
    puts("PASS: clean, two outliers, deterministic, insufficient, singular, invalid variance");
    return 0;
}
