/* Bounded real-data regression runner. Outputs only to the supplied path.
 * Online MC must be disabled in the linked objects for this smoke test. */
#include "rtklib.h"

int main(int argc, char **argv)
{
    prcopt_t p=prcopt_default;
    solopt_t s=solopt_default;
    filopt_t files={""};
    double ep0[6]={2025,4,25,7,47,0},ep1[6]={2025,4,25,7,48,0};
    char *inputs[]={"F:\\RTKLib-LAB\\data\\Static\\rover.25O",
                    "F:\\RTKLib-LAB\\data\\Static\\rover.25P",
                    "F:\\RTKLib-LAB\\data\\Static\\base.25O"};
    if (argc!=2) return 2;
    p.mode=PMODE_KINEMA; p.nf=2; p.elmin=15*D2R;
    p.weightmode=WEIGHTOPT_SNR; p.ionoopt=IONOOPT_BRDC;
    p.tropopt=TROPOPT_SAAS; p.soltype=0; p.intpref=0;
    p.navsys=SYS_GPS+SYS_CMP;
    p.modear=1; p.bdsmodear=1; p.gpsmodear=1; p.glomodear=0;
    p.minlock=0; p.minfix=10; p.maxout=3; p.chkMaxPLD=120;
    p.thresar[0]=3; p.dynamics=0; p.refpos=3;
    p.eratio[0]=p.eratio[1]=100; p.err[0]=100;
    p.err[1]=p.err[2]=0.003; p.std[0]=30;
    p.intmon=0; p.intoutbias=0;
    p.intfde=0;
#if ENABLE_RTK_RANSAC_ARAIM || ENABLE_RTK_RANSAC_SLOPE
    p.intmon=1; p.modear=ARMODE_OFF; p.niter=1;
    p.intpfa=1E-4; p.intpmd=1E-3;
#endif
    s.timef=1; s.posf=SOLF_XYZ; s.trace=3; s.sstat=0;
    return postpos(epoch2time(ep0),epoch2time(ep1),1.0,0.0,&p,&s,&files,
                   inputs,3,argv[1],"","");
}
