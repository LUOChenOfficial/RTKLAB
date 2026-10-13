/*------------------------------------------------------------------------------
* rtk_int_internal.h : private RTK integrity interfaces
*-----------------------------------------------------------------------------*/
#ifndef RTK_INT_INTERNAL_H
#define RTK_INT_INTERNAL_H

#include "rtklib.h"

#define RTKINT_DEFAULT_PFA 1e-6
#define RTKINT_DEFAULT_PMD 5e-8

#ifdef ENABLE_RTK_INTEGRITY
double rtkimprob(double p, double def);
double rtkimnorminv(double p);
void rtkimthres(const prcopt_t *opt, int nset, double *kmain, double *ksub,
                double *kss);
void rtkimsigenu(const double *rr, const double *qr, double *sig);
int rtkimsat(int flag);
void rtkimslassess(rtk_t *rtk);
void rtkimdbg(const char *format, ...);
#endif

#endif
