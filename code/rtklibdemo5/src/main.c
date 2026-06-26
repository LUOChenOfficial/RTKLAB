#include "rtklib.h"
#include <stdarg.h>

#define PROGNAME    "rnx2rtkp"          /* program name */
#define MAXFILE     16                  /* max number of input files */

/* help text -----------------------------------------------------------------*/
static const char *help[]={
"",
" usage: rnx2rtkp [option]... file file [...]",
"",
" Read RINEX OBS/NAV/GNAV/HNAV/CLK, SP3, SBAS message log files and ccompute ",
" receiver (rover) positions and output position solutions.",
" The first RINEX OBS file shall contain receiver (rover) observations. For the",
" relative mode, the second RINEX OBS file shall contain reference",
" (base station) receiver observations. At least one RINEX NAV/GNAV/HNAV",
" file shall be included in input files. To use SP3 precise ephemeris, specify",
" the path in the files. The extension of the SP3 file shall be .sp3 or .eph.",
" All of the input file paths can include wild-cards (*). To avoid command",
" line deployment of wild-cards, use \"...\" for paths with wild-cards.",
" Command line options are as follows ([]:default). With -k option, the",
" processing options are input from the configuration file. In this case,",
" command line options precede options in the configuration file.",
"",
" -?        print help",
" -k file   input options from configuration file [off]",
" -o file   set output file [stdout]",
" -ts ds ts start day/time (ds=y/m/d ts=h:m:s) [obs start time]",
" -te de te end day/time   (de=y/m/d te=h:m:s) [obs end time]",
" -ti tint  time interval (sec) [all]",
" -p mode   mode (0:single,1:dgps,2:kinematic,3:static,4:moving-base,",
"                 5:fixed,6:ppp-kinematic,7:ppp-static) [2]",
" -m mask   elevation mask angle (deg) [15]",
" -sys s[,s...] nav system(s) (s=G:GPS,R:GLO,E:GAL,J:QZS,C:BDS,I:IRN) [G|R]",
" -f freq   number of frequencies for relative mode (1:L1,2:L1+L2,3:L1+L2+L5) [2]",
" -v thres  validation threshold for integer ambiguity (0.0:no AR) [3.0]",
" -b        backward solutions [off]",
" -c        forward/backward combined solutions [off]",
" -i        instantaneous integer ambiguity resolution [off]",
" -h        fix and hold for integer ambiguity resolution [off]",
" -e        output x/y/z-ecef position [latitude/longitude/height]",
" -a        output e/n/u-baseline [latitude/longitude/height]",
" -n        output NMEA-0183 GGA sentence [off]",
" -g        output latitude/longitude in the form of ddd mm ss.ss' [ddd.ddd]",
" -t        output time in the form of yyyy/mm/dd hh:mm:ss.ss [sssss.ss]",
" -u        output time in utc [gpst]",
" -d col    number of decimals in time [3]",
" -s sep    field separator [' ']",
" -r x y z  reference (base) receiver ecef pos (m) [average of single pos]",
"           rover receiver ecef pos (m) for fixed or ppp-fixed mode",
" -l lat lon hgt reference (base) receiver latitude/longitude/height (deg/m)",
"           rover latitude/longitude/height for fixed or ppp-fixed mode",
" -y level  output soltion status (0:off,1:states,2:residuals) [0]",
" -x level  debug trace level (0:off) [0]"
};


/* show message --------------------------------------------------------------*/
extern int showmsg(char *format, ...)
{
    va_list arg;
    va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
    fprintf(stderr,"\r");
    return 0;
}
extern void settspan(gtime_t ts, gtime_t te) {};
extern void settime(gtime_t time) {};


/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i=0;i<(int)(sizeof(help)/sizeof(*help));i++) fprintf(stderr,"%s\n",help[i]);
    exit(0);
}


/* 定位测试主函数 -------------------------------------------------------*/
int main(int argc,char *argv) 
{
	gtime_t ts = { 0 }, te = { 0 };
	prcopt_t prcopt = prcopt_default; 
	solopt_t solopt = solopt_default;
	filopt_t filopt={""};
	double tint=1.0; //采样间隔记得改 
	int n=3; //输入文件数
    char basepath[] = "F:\\RTKLib-LAB\\data";
    char resultpath[] = "F:\\RTKLib-LAB\\result";

    //////////////////////////////////////////////////////注意更改文件夹
    char str_envir[] = "\\Build\\";

    char* str[5];
    for (int i = 0; i < 5; i++) {
        str[i] = (char*)malloc(100 * sizeof(char));
        strcpy(str[i], basepath);
        str[i]= strcat(str[i], str_envir);
    }

    //记得改文件名尤其是24/25O
    char* infile[3] = { {strcat(str[0],"rover.21O")},{strcat(str[1],"rover.21P")}};
    char outfile[256] = "";
    infile[2] = strcat(str[3],"base.21O");
    strcpy(outfile, resultpath);
    strcat(outfile, str_envir);
    strcat(outfile, "test.pos");

    //openmodel(strcat(str[4], "model_weight.txt"));


	/* 处理参数设置(与RTKLIB保持一致，仅更改与默认不一致参数) */
	/* Setting1 */
	prcopt.mode = PMODE_KINEMA; /* {PMODE_SINGLE;PMODE_KINEMA; PMODE_DGPS;PMODE_STATIC} */
	prcopt.nf= 2;  //频率
    prcopt.elmin = 10 * D2R;//截止高度角（弧度）
    prcopt.weightmode=WEIGHTOPT_ELEVATION;
    //prcopt.weightmode=WEIGHTOPT_SNR;
    //prcopt.weightmode = WEIGHTOPT_PLD;                   //包含PLD观测值才可以使用
    //prcopt.weightmode = WEIGHTOPT_CV_PV;
    //prcopt.weightmode = WEIGHTOPT_ELE_SNR;
    //prcopt.weightmode = WEIGHTOPT_ELE_PLD;           //包含PLD观测值才可以使用
    //prcopt.weightmode = WEIGHTOPT_SNR_PLD;           //包含PLD观测值才可以使用
    //prcopt.weightmode = WEIGHTOPT_ALL;                     //包含PLD观测值才可以使用

	prcopt.ionoopt = IONOOPT_BRDC;
	prcopt.tropopt = TROPOPT_SAAS;
	prcopt.soltype = 0;
    prcopt.navsys = SYS_GPS +SYS_CMP;

	/* Setting2 */
	prcopt.modear=1;
	prcopt.bdsmodear=1;
    prcopt.gpsmodear = 1;
    prcopt.glomodear = 0;//模糊度固定设置
    prcopt.minlock = 0; //不固定某颗卫星的历元数,新卫星初始化以及maxout
    prcopt.minfix = 10; //holdamb的最小固定历元数
    prcopt.maxout = 3;//连续3个历元卫星不可用则初始化模糊度
    prcopt.chkMaxPLD = 120; //PLD最大阈值
	prcopt.thresar[0]=2.5;
    prcopt.dynamics = 0;
    
#ifdef ENABLE_RTK_INTEGRITY
    /* Integrity master switches */
    prcopt.enable_rtk_integrity_monitor = 0;
    prcopt.enable_rtk_integrity_rbias_export = 0;

    /* Monitored fault modes */
    prcopt.enable_monitor_single_satellite_fault = 1;
    prcopt.enable_monitor_fixed_ambiguity_update_fault = 1;

    /* Recovery / subset control */
    prcopt.enable_rtk_integrity_fde_recovery = 1;
    prcopt.rtk_integrity_max_subset_filters = 32;

    /* Integrity debug output */
    prcopt.enable_rtk_integrity_subset_debug_output = 0;
    prcopt.rtk_integrity_debug_subset_id = 0;
    prcopt.rtk_integrity_debug_satellite = 0;
    prcopt.rtk_integrity_debug_pl_threshold = 0.0;
#endif
#if ENABLE_RTK_SKIP_EPOCH
    prcopt.rtk_skip_epoch_time[0] = 2021;
    prcopt.rtk_skip_epoch_time[1] = 2;
    prcopt.rtk_skip_epoch_time[2] = 4;
    prcopt.rtk_skip_epoch_time[3] = 9;
    prcopt.rtk_skip_epoch_time[4] = 45;
    prcopt.rtk_skip_epoch_time[5] = 2;
    prcopt.rtk_skip_epoch_satellite = 108; /* sat号 */
#endif

	/* Position */
	prcopt.refpos=3;//0
	if (prcopt.refpos==0) {
        //////////////////////////////////////////////////////注意更改基准站坐标
        double rb[3] = { -2191836.0667,5182445.0573,2993201.7774 };
		prcopt.rb[0]= rb[0];
		prcopt.rb[1]= rb[1];
		prcopt.rb[2]= rb[2];
	}
	/* Statistics */
	prcopt.eratio[0]=100.0;
	prcopt.eratio[1]=100.0;
	prcopt.err[0] = 100.0;
	prcopt.err[1]=0.003;
	prcopt.err[2]=0.003;
	prcopt.std[0]=30.0;

    /*Statistical test*/
    prcopt.sizes[0] = 4;
    prcopt.sizes[1] = 4;

    /*Open*/
    static int arr1[] = { 5, 16, 24, 30 };  // GPS dd fixed satellite list -- GPS PRN+2
    static int arr2[] = { 134, 136, 140, 149 };  // BDS dd fixed satellite list -- BDS PRN+104
    /*Build*/
    //static int arr1[] = { 3, 10, 23, 32 };
    //static int arr2[] = { 114, 116, 128, 129 };
    /*Wall*/
    //static int arr1[] = { 5, 16, 19, 30 };
    //static int arr2[] = { 111, 114, 115, 129 };
    /*Tree*/
   // static int arr1[] = { 3, 9, 10, 29 };
    //static int arr2[] = { 114, 116, 128, 129 };

    prcopt.satlist[0] = arr1; prcopt.satlist[1] = arr2;

    
	/* Output */
	solopt.timef= 1;
	solopt.posf = SOLF_XYZ;
	solopt.sstat = 0;
	solopt.trace = 0;

	postpos(ts, te, tint, 0.0, &prcopt, &solopt, &filopt, infile, n ,outfile, "", "");
    closemodel();
    // 释放内存
    for (int i = 0; i < 5; i++) {
        free(str[i]);
    }
}



/* 各种定位功能主要代码（参考rnx2rtkp.c rnx2rtkp main） --------------------*/
/*
int main(int argc, char **argv)
{
    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    gtime_t ts={0},te={0};
    double tint=0.0,es[]={2000,1,1,0,0,0},ee[]={2000,12,31,23,59,59},pos[3];
    int i,j,n,ret;
    char *infile[MAXFILE],*outfile="",*p;
    
    prcopt.mode  =PMODE_KINEMA;
    prcopt.navsys=0;
    prcopt.refpos=1;
    prcopt.glomodear=1;
    solopt.timef=0;
    sprintf(solopt.prog ,"%s ver.%s %s",PROGNAME,VER_RTKLIB,PATCH_LEVEL);
    sprintf(filopt.trace,"%s.trace",PROGNAME);
    
    // load options from configuration file
    for (i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-k")&&i+1<argc) {
            resetsysopts();
            if (!loadopts(argv[++i],sysopts)) return -1;
            getsysopts(&prcopt,&solopt,&filopt);
        }
    }
    for (i=1,n=0;i<argc;i++) {
        if      (!strcmp(argv[i],"-o")&&i+1<argc) outfile=argv[++i];
        else if (!strcmp(argv[i],"-ts")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",es,es+1,es+2);
            sscanf(argv[++i],"%lf:%lf:%lf",es+3,es+4,es+5);
            ts=epoch2time(es);
        }
        else if (!strcmp(argv[i],"-te")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",ee,ee+1,ee+2);
            sscanf(argv[++i],"%lf:%lf:%lf",ee+3,ee+4,ee+5);
            te=epoch2time(ee);
        }
        else if (!strcmp(argv[i],"-ti")&&i+1<argc) tint=atof(argv[++i]);
        else if (!strcmp(argv[i],"-k")&&i+1<argc) {++i; continue;}
        else if (!strcmp(argv[i],"-p")&&i+1<argc) prcopt.mode=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-f")&&i+1<argc) prcopt.nf=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-sys")&&i+1<argc) {
            for (p=argv[++i];*p;p++) {
                switch (*p) {
                    case 'G': prcopt.navsys|=SYS_GPS;
                    case 'R': prcopt.navsys|=SYS_GLO;
                    case 'E': prcopt.navsys|=SYS_GAL;
                    case 'J': prcopt.navsys|=SYS_QZS;
                    case 'C': prcopt.navsys|=SYS_CMP;
                    case 'I': prcopt.navsys|=SYS_IRN;
                }
                if (!(p=strchr(p,','))) break;
            }
        }
        else if (!strcmp(argv[i],"-m")&&i+1<argc) prcopt.elmin=atof(argv[++i])*D2R;
        else if (!strcmp(argv[i],"-v")&&i+1<argc) prcopt.thresar[0]=atof(argv[++i]);
        else if (!strcmp(argv[i],"-s")&&i+1<argc) strcpy(solopt.sep,argv[++i]);
        else if (!strcmp(argv[i],"-d")&&i+1<argc) solopt.timeu=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-b")) prcopt.soltype=1;
        else if (!strcmp(argv[i],"-c")) prcopt.soltype=2;
        else if (!strcmp(argv[i],"-i")) prcopt.modear=2;
        else if (!strcmp(argv[i],"-h")) prcopt.modear=3;
        else if (!strcmp(argv[i],"-t")) solopt.timef=1;
        else if (!strcmp(argv[i],"-u")) solopt.times=TIMES_UTC;
        else if (!strcmp(argv[i],"-e")) solopt.posf=SOLF_XYZ;
        else if (!strcmp(argv[i],"-a")) solopt.posf=SOLF_ENU;
        else if (!strcmp(argv[i],"-n")) solopt.posf=SOLF_NMEA;
        else if (!strcmp(argv[i],"-g")) solopt.degf=1;
        else if (!strcmp(argv[i],"-r")&&i+3<argc) {
            prcopt.refpos=prcopt.rovpos=0;
            for (j=0;j<3;j++) prcopt.rb[j]=atof(argv[++i]);
            matcpy(prcopt.ru,prcopt.rb,3,1);
        }
        else if (!strcmp(argv[i],"-l")&&i+3<argc) {
            prcopt.refpos=prcopt.rovpos=0;
            for (j=0;j<3;j++) pos[j]=atof(argv[++i]);
            for (j=0;j<2;j++) pos[j]*=D2R;
            pos2ecef(pos,prcopt.rb);
            matcpy(prcopt.ru,prcopt.rb,3,1);
        }
        else if (!strcmp(argv[i],"-y")&&i+1<argc) solopt.sstat=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-x")&&i+1<argc) solopt.trace=atoi(argv[++i]);
        else if (*argv[i]=='-') printhelp();
        else if (n<MAXFILE) infile[n++]=argv[i];
    }
    if (!prcopt.navsys) {
        prcopt.navsys=SYS_GPS|SYS_GLO;
    }
    if (n<=0) {
        showmsg("error : no input file");
        return -2;
    }
    ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,"","");
    
    if (!ret) fprintf(stderr,"%40s\r","");
    return ret;
}*/
