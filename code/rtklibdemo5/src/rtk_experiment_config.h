/* Experimental build profiles. Override defaults with compiler definitions.
 * These switches retain the existing estimator parameters and output behavior. */
#ifndef RTK_EXPERIMENT_CONFIG_H
#define RTK_EXPERIMENT_CONFIG_H

#ifndef ENABLE_RTK_INTEGRITY
#define ENABLE_RTK_INTEGRITY
#endif
#define RTK_INT_METHOD_FILTER_SS 0
#define RTK_INT_METHOD_LS_SS     1
#define RTK_INT_METHOD_SLOPE     2
#define RTK_INT_METHOD_GLS_SS    3
/* Experimental combinations: change the master switch, then rebuild. */
#ifndef ENABLE_RTK_RANSAC_ARAIM
#define ENABLE_RTK_RANSAC_ARAIM 0 /* screened single KF + current-epoch GLS subsets */
#endif
#ifndef ENABLE_RTK_RANSAC_SLOPE
#define ENABLE_RTK_RANSAC_SLOPE 0 /* screened float KF + innovation slope PL */
#endif
#if ENABLE_RTK_RANSAC_ARAIM && ENABLE_RTK_RANSAC_SLOPE
#error Select only one RANSAC integrity combination: ARAIM or SLOPE
#endif
#ifndef RTK_INT_METHOD
#if ENABLE_RTK_RANSAC_SLOPE
#define RTK_INT_METHOD RTK_INT_METHOD_SLOPE
#elif ENABLE_RTK_RANSAC_ARAIM
#define RTK_INT_METHOD RTK_INT_METHOD_GLS_SS
#else
#define RTK_INT_METHOD RTK_INT_METHOD_GLS_SS /* retain current GLS monitoring without RANSAC */
#endif
#endif
#ifndef ENABLE_RTK_DEBUG_OUTPUT
#define ENABLE_RTK_DEBUG_OUTPUT 1
#endif
#ifndef ENABLE_RTK_ONLINE_MC
#define ENABLE_RTK_ONLINE_MC 1
#endif
#ifndef RTK_ONLINE_MC_SAMPLES
#define RTK_ONLINE_MC_SAMPLES 10000
#endif
#ifndef RTK_ONLINE_MC_RATIO
#define RTK_ONLINE_MC_RATIO 3.0
#endif
#ifndef ENABLE_RTK_SKIP_EPOCH
#define ENABLE_RTK_SKIP_EPOCH 0
#endif
#ifndef ENABLE_RTK_RANSAC
#define ENABLE_RTK_RANSAC (ENABLE_RTK_RANSAC_ARAIM || ENABLE_RTK_RANSAC_SLOPE)
#endif
#if ENABLE_RTK_RANSAC_ARAIM && (!ENABLE_RTK_RANSAC || (RTK_INT_METHOD!=RTK_INT_METHOD_FILTER_SS && RTK_INT_METHOD!=RTK_INT_METHOD_GLS_SS))
#error RANSAC_ARAIM requires RANSAC and FILTER_SS or GLS_SS
#endif
#if ENABLE_RTK_RANSAC_SLOPE && (!ENABLE_RTK_RANSAC || RTK_INT_METHOD!=RTK_INT_METHOD_SLOPE)
#error RANSAC_SLOPE requires RANSAC and the SLOPE integrity method
#endif
#ifndef RTK_RANSAC_MAX_ITER
#define RTK_RANSAC_MAX_ITER 1000
#endif
/* Selected ARAIM profile: split frequency, prior variance, no consensus refit.
 * Explicit compile definitions still allow the historical ablations. */
#ifndef RTK_RANSAC_SPLIT_FREQUENCY
#define RTK_RANSAC_SPLIT_FREQUENCY ENABLE_RTK_RANSAC_ARAIM
#endif
#ifndef RTK_RANSAC_PRIOR_VARIANCE
#define RTK_RANSAC_PRIOR_VARIANCE ENABLE_RTK_RANSAC_ARAIM
#endif
#ifndef RTK_RANSAC_REFIT
#define RTK_RANSAC_REFIT 0 /* selected: no final-consensus refit */
#endif
#ifndef RTK_RANSAC_SIGMA
#define RTK_RANSAC_SIGMA 12.0 /* RANSACRTK::multiThreshold */
#endif
#ifndef RTK_RANSAC_PROBABILITY
#define RTK_RANSAC_PROBABILITY 0.99
#endif
#ifndef ENABLE_BDS_NEW_DATA
#define ENABLE_BDS_NEW_DATA 0 /* 1: BDS-3-reorganized data */
#endif
#if ENABLE_BDS_NEW_DATA
#define BDS_GEO_PRN(prn)              ((prn)>=1&&(prn)<=4)
#define BDS_GEO_FALLBACK_PRN(prn)     BDS_GEO_PRN(prn)
#define BDS_EXCLUDE_OBS_PRN(prn)      ((prn)==9||(prn)==10)
#else
#define BDS_GEO_PRN(prn)              ((prn)>=59&&(prn)<=62)
#define BDS_GEO_FALLBACK_PRN(prn)     BDS_GEO_PRN(prn)
#define BDS_EXCLUDE_OBS_PRN(prn)      ((prn)<19)
#endif
#ifndef ENABLE_RTK_ARAIM_PL_BIAS_TERM
#define ENABLE_RTK_ARAIM_PL_BIAS_TERM 1
#endif
#ifndef ENABLE_RTK_INTEGRITY_CHI2_GATE
#define ENABLE_RTK_INTEGRITY_CHI2_GATE 1
#endif
#ifndef ENABLE_RTK_PARTIAL_AR
#define ENABLE_RTK_PARTIAL_AR 1
#endif
#ifndef ENABLE_RTK_AR_CP_RMS_GATE
#define ENABLE_RTK_AR_CP_RMS_GATE 1
#endif
#ifndef RTK_AR_CP_RMS_GATE_THRESH
#define RTK_AR_CP_RMS_GATE_THRESH 0.05 //m
#endif

#endif
