/*------------------------------------------------------------------------------
* rtk_int_config.c : RTK integrity processing options
*
* notes  : Apply the selected build profile after normal positioning options.
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#ifdef ENABLE_RTK_INTEGRITY
/* configure RTK integrity experiment ----------------------------------------
* args   : prcopt_t *opt   IO  processing options
* return : none
* notes  : Apply after ordinary positioning settings. The selected profile
*          can override AR mode, KF iterations, FDE and risk parameters.
*-----------------------------------------------------------------------------*/
EXPORT void rtkimcfg(prcopt_t *opt)
{
    /* Integrity master switches */
    opt->intmon=ENABLE_RTK_RANSAC_ARAIM||ENABLE_RTK_RANSAC_SLOPE||
                  RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS;
    opt->intoutbias=0;

    /* Monitored fault modes */
    opt->intsat=1;
    opt->intamb=1;

    /* Recovery / subset control */
    opt->intfde=1;
    opt->maxintsub=32;

    /* Integrity debug output */
    opt->intoutsub=0;
    opt->intdbgid=0;
    opt->intdbgsat=0;
    opt->intdbgpl=0.0;
    opt->intpfa=1e-4;
    opt->intpmd=1e-3;
#if ENABLE_RTK_RANSAC_ARAIM||ENABLE_RTK_RANSAC_SLOPE||RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS
    /* GLS protects float/fixed output with unconstrained float subsets.
     * The persistent KF stays float; never use fix-and-hold in this mode. */
    opt->modear=
        RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS?ARMODE_CONT:ARMODE_OFF;
    opt->niter=1;
#endif
#if ENABLE_RTK_RANSAC_ARAIM||RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS
    /* GLS streams single/double satellite and wrong-integer hypotheses.
     * FILTER_SS remains available explicitly for the parallel-filter baseline. */
    opt->intamb=RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS;
    opt->intfde=RTK_INT_METHOD==RTK_INT_METHOD_GLS_SS?RTK_GLS_FDE:0;
    opt->maxintsub=MAXSAT+1;
    opt->intpfa=RTK_GLS_PFA;
    opt->intpmd=RTK_GLS_PMD;
#endif
}
#endif
