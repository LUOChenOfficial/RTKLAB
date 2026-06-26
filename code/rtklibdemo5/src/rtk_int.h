#ifndef RTK_INT_H
#ifndef RTKLIB_IN_RTK_INT
#include "rtklib.h"
#else
#define RTK_INT_H

#define RTKIM_F_NONE   0
#define RTKIM_F_SAT    1
#define RTKIM_F_NFIX   2

#define RTKIM_A_NONE   0
#define RTKIM_A_SAT    1
#define RTKIM_A_FLOAT  2
#define RTKIM_MAX_BIAS_ROWS (MAXOBS*NFREQ*2+1)

struct rtk_tag;

typedef struct {
    int id;
    int mode;
    int sat;
    int fixed;
} rtkim_def_t;

typedef struct {
    int init,act,ready,valid;
    struct rtk_tag *rtk;
    double pos[3];
    double q[6];
    double sig[3];
    int qi;
    double ratio;
    double pdop,hdop,vdop;
    double sep[3];
    double sig0,chisq;
    int dof;
    int refsat[8];
    int ambc;
} rtkim_ss_t;

typedef struct {
    gtime_t time;
    int nu,nr,ns;
    int qi;
    double ratio;
    int fixed;
    int fixed_upd;
} rtkim_ep_t;

typedef struct {
    double hpl,vpl;
    double hpl0,vpl0;
    double hadd,vadd;
    double be,bn,bu;
    double pe,pn,pu;
    double msig[3];
    int bias_loaded,bias_rows;
    int hsrc,hmode,hsat;
    int vsrc,vmode,vsat;
    double hsig[3],vsig[3];
    double hsep[3],vsep[3];
} rtkim_pl_t;

typedef struct {
    int mode;
    int act;
    int id;
    int sat;
    double score;
    double stat[3];
    double thres[3];
} rtkim_fde_t;

typedef struct {
    int ena,init;
    int ndef,nact;
    rtkim_def_t def[MAXSAT+1];
    rtkim_ss_t ss[MAXSAT+1];
    rtkim_ep_t ep;
    rtkim_pl_t pl;
    rtkim_fde_t fde;
    rtkim_fde_t act;
    int child;
    int retry;
    int bnx,bnv;
    int bflg[RTKIM_MAX_BIAS_ROWS];
    double *bH,*bR,*bv;
} rtkim_t;

#endif
#endif

