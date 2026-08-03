/* -------------------------------------------------------------------------
* SUPREME_Options.h : Process Configuration Options in Multi-GNSS Solutions
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2018.12.07
* ------------------------------------------------------------------------- */
#ifndef SUPREME_OPTIONS_H_HH
#define SUPREME_OPTIONS_H_HH

#include "SUPREME_GnssTime.h"
#include "SUPREME_Defvar.h"
#include <vector>


/// Porcess System: PRO_SYS_???
#define PRO_SYS_GPS                1                       // Process system G: GPS      (America)
#define PRO_SYS_GLO                2                       // Process system R: GLONASS  (Russa)
#define PRO_SYS_GAL                4                       // Process system E: GALILEO  (Europ)
#define PRO_SYS_BDS                8                       // Process system C: BDS      (China)
#define PRO_SYS_QZS               32                       // Process system J: QZSS     (Japan)
#define PRO_SYS_SBS               64                       // Process system S: SBAS     
#define PRO_SYS_IRN              128                       // Process system  : IRNSS    (India)

#define MODE_POST                  0                       // Post processing mode
#define MODE_REALTIME              1                       // Realtime processing mode

#define MODE_SPP                   0                       // Standard point positioning
#define MODE_PPP                   1                       // Precise point positioning
#define MODE_RPP                   2                       // Relative point positioning
#define MODE_NET                   3                       // GNSS net processing

#define MODE_KINEMATIC             0                       // Kinematic process
#define MODE_STATIC                1                       // Static process

#define MODE_FILTER_FORWARD        0                       // Forward filting mode
#define MODE_FILTER_BACKWARD       1                       // Backward filting mode
#define MODE_FILTER_COMBINE        2                       // Combine filting mode

#define MODE_IF_COMBINE            0                       // Uncombine model
#define MODE_UNCOMBINE             1                       // Ionosphere-free combine model
#define MODE_UOFC_COMBINE          2                       // Uofc combine model
#define MODE_UNCOMBINE_ION_CONT    3						//电离层GIM约束
#define GIM_const				   1
#define P1P2_const				   2	

#define MODE_PARAEST_LSQFILTING    0                       // Least square filting
#define	MODE_PARAEST_KALMANFILTING 1                       // Kalman filting
#define MODE_PARAEST_RMSFILTING    2                       // Root mean square filting

#define MODE_PPP_FLOAT	           0                       // PPP float solutions
#define MODE_PPP_FIXED             1                       // PPP ambiguity resolution

#define TROP_EST_NONE              0                       // Not estimate troposphere delay
#define TROP_EST_WET               1                       // Estimate trop wet delay

/// Trop Model
#define TROP_MODEL_HOPFIELD        0
#define TROP_MODEL_SAASTAMOINEN    1
#define TROP_MODEL_UNB3            2
#define TROP_MODEL_GPT2            3
#define TROP_MODEL_GPT2w           4
#define TROP_MODEL_IGGTROP         5

#define TROP_MF_NMF                0
#define TROP_MF_GMF                1
#define TROP_MF_VMF                2

/// Ion Model
#define ION_MODEL_FREE             0                        // Ionosphere-free model
#define ION_MODEL_EST_WHITE        1                        // Estimate ionosphere as white noise
#define ION_MODEL_EST_RANDWALK     2                        // Estimate ionosphere as random walk
#define ION_MODEL_GIM_CORRECT      3                        // Correct ionosphere by GIM
#define ION_MODEL_GIM_DCB_CORRECT  4                        // Corrected by GIM and SatDcb
#define ION_MODEL_GIM_CONSTRAINT   5                        // Estimate ionosphere by GIM as psuedo observation


#define RCLK_MODEL_WHITE_NOISE     0                        // Estimate receiver clock as white noise
#define RCLK_MODEL_RANDWALK        1                        // Estimate receiver clock as random walk

#define THREADMAXNUM               6              // Max Thread Number

#define OUTPATH_LENGTH             400            // Out name
#define PATH_LENGTH                400            // File path length

/// General File Path
struct GenFilePath
{
	char trp_fl[PATH_LENGTH];      /* TROP file path               */
	char ion_fl[PATH_LENGTH];      /* ION file path                */

	char atx_fl[PATH_LENGTH];      /* ATX file path                */
	char pco_fl[PATH_LENGTH];      /* PCO information file         */
	char pcv_fl[PATH_LENGTH];      /* PCV information file         */

	//char pla_fl[PATH_LENGTH];      /* DE405 planet ephemeris       */
	//char sec_fl[PATH_LENGTH];      /* UTC-TAI information file     */
	//char iau_fl[PATH_LENGTH];      /* IAU2000 file                 */

	char dcb_p1c1_fl[PATH_LENGTH]; /* DCB-P1C1                     */
	char dcb_p1p2_fl[PATH_LENGTH]; /* DCB-P1P2                     */
	char dcb_p2c2_fl[PATH_LENGTH]; /* DCB-P2C2                     */
	char dcb_mgex[PATH_LENGTH];    /* DCB-mgex  for BDS                     */
	char ifcb_f[PATH_LENGTH];

	char erp_fl[PATH_LENGTH];      /* ERP file path                */
	char blq_fl[PATH_LENGTH];      /* BLQ file path                */
	char crd_fl[PATH_LENGTH];      /* station coordinate file      */
};

/// Rinex File Path
struct RinexFilePath
{
	char  obsf[PATH_LENGTH];       /* observation data file path   */
	char  navf[PATH_LENGTH];       /* navigation data file path    */
	char  prephf[PATH_LENGTH];     /* precise ephemeris file path  */
	char  preclkf[PATH_LENGTH];    /* precise clock file path      */
};

/// Gnss Output File Path
struct GnssOutputFilePath
{
	char time_str[20];                   /* time in string                  */

	char OutName[OUTPATH_LENGTH];        /* out name                        */

	char out_fp[PATH_LENGTH];            /* output positon etc parameter    */
	char out_ion[PATH_LENGTH];           /* output ionosphere slant delay   */
	char out_res[PATH_LENGTH];           /* output residuals                */
	char out_vbl[PATH_LENGTH];           /* output variables                */

	char out_sat[PATH_LENGTH];           /* output satellite use and delete */
	char out_log[PATH_LENGTH];           /* output log information          */
	char out_trace[PATH_LENGTH];         /* output trace                    */
	char out_cycle_slip[PATH_LENGTH];    /* satellite cycle slip            */

	/// out full path
	char out_fp_fpath[PATH_LENGTH];      /* output positon etc parameter full path    */
	char out_ion_fpath[PATH_LENGTH];     /* output ionosphere delay                   */
	char out_res_fpath[PATH_LENGTH];     /* output residuals full path                */
	char out_vbl_fpath[PATH_LENGTH];     /* output variables full path                */

	char out_sat_fpath[PATH_LENGTH];     /* output satellite use and delete full path */
	char out_log_fpath[PATH_LENGTH];     /* output log information full path          */
	char out_trace_fpath[PATH_LENGTH];   /* output trace full path                    */
	char out_IFB[PATH_LENGTH];
};

/// Gnss File Path
struct FilePath
{
	RinexFilePath      rinex_f;
	GenFilePath        gen_f;
	GnssOutputFilePath out_f;
};

struct estpara_t
{
	unsigned int xyz;
	unsigned int rclk;
	unsigned int trop;
	unsigned int ion;
	unsigned int amb;
};

struct ppp_option_t
{
	char ProjName[300];              // Project name
	char StaName[50];                // Station name
	unsigned int StaNo;              // Station ID No.

	unsigned int ion_const;			 //电离层约束
	unsigned int trop_model;         // troposphere model: 0: Estimate (random walk) 1: Gradient 2: Grid or ZPD corrected
	unsigned int ion_model;          // ionosphere model: 0:IF model 1: white noise 2:random walk 3:GIM corrected 4: GIM and SDCB corrected 5: GIM as psuedo obs
	unsigned int clk_model;          // receiver clock model: 0: white noise  1: random-walk noise (only for external maser)
	unsigned int math_model;         // math model: 0: IF model 1: uncombined model

	unsigned int GIM_P1P2;
	unsigned int IS_ISB;
	unsigned int IS_DCB;
	unsigned int ISB_MODE;

	unsigned int mode_ppp;           // process mode: 0: static  1: kinematic
	unsigned int mode_filter;        // filter mode： 0: forward 1: backward 2: combined
	unsigned int process_time;       // 0: post process 1: realtime process
	unsigned int process_mode;       // process mode: 0: spp 1: ppp 2: relative positioning 3: net process
	unsigned int System;             // system need to process

	unsigned int Gnuplot;            // graph plot by gnuplot
	unsigned int trace;              // trace mode

	unsigned int ssr_id;
	unsigned int nav_id;
	unsigned int use_orbclk;
	char ssr_mountpoint[50];
	unsigned int ssr_com;            // SSR Orb: 0->APC   1->COM
	char nav_mountpoint[50];
	unsigned int Sta_id;             // GNSS data id
	char IGGNtripPath[500];          // IGGNtrip path in linux
	double time_delay;

	unsigned int interval;           // sampling interval
	unsigned int freqn;              // frequency number
	unsigned int Freqn;

	double ele_mask;                 // satellite elevation mask
	double sigmaP;                   // psuedo-range prior variance
	double sigmaL;                   // carrier phase prior variance

	double antenna_offset[3];        // antenna offset:NEU
	char antenna_type[21];           // receiver antenna type
	char receiver_type[21];          // receiver type

	double trop_noise, trop_gradientNnoise, trop_gradientEnoise; // troposphere process noise
	double ion_noise;                // ionosphere process noise
	double rclk_noise;               // receiver clock noise

	estpara_t estpara;
	
	FilePath IOfile;                 // Input and Output file path

	unsigned int threadn;            // thread number
	//unsigned int batchmode;          // batch process mode
	unsigned int ContinueMode;       // Continue process mode
	unsigned int para_est_mode;      // parameter estimation mode
	char stationlist[PATH_LENGTH];   // station list path
	char CentreName[10];             // analysis centre
	char ObsPath[PATH_LENGTH];       // Observation path
	char PrePath[PATH_LENGTH];       // Precise path
	char GenPath[PATH_LENGTH];       // General file path
	gnsstime starttime, endtime;     // batch process start and end time

	gnsstime CurrentEpoch;           // current process epoch time
	int dt;                          // current time - previous time

	double nmea, lat, lon;
	double GDOP;

	unsigned int SatN_SPP, SatN_PPP;
	double StaPos_CRD[3];
	double StaPos[3];
	double StaPos_SPP[3];
	double StaPos_PPP[3];
	double ProXYZ[3];
	double PreXYZ[3];
	double Rclk;
	double ISB[4];
	double ISB_pre[4];
	char ISBsystem[4];
	char ISBsystem_pre[4];
	unsigned int isbsit[4];
	unsigned int isbsit_pre[4];

	double SunPos[3], MoonPos[3];

	double phw1[500], phw2[500];

	char Sysm[255];
	int SysNum;

	unsigned int epochnum;

	int clkJump;                       // receiver clock jump
	double obs0[GNSS_SATNO_NUM][4];    //GPS_SATNUM

	int QC_Prn[GNSS_SATNO_NUM];            // Cycle Slip flag
	int Time_dif[GNSS_SATNO_NUM];           //gf 时间距离
	double gf_pre[GNSS_SATNO_NUM];

	double repairvalue[repairstanum][3];
	double repairsimvalue[repairstanum][3];
	int recordrepair[repairstanum];
	int recordsim[repairstanum];

	int recordfre[10];
	int Recordfre;

	//模拟周跳序列（每个历元只添加一个周跳）
	int simulcycletime[sinulstanum];
	int simulcyclesit[sinulstanum];
	int simulcyclevalue[sinulstanum][3];
	int simnum;
	int simnumc;
	int simnums;

	//猜Z检验量数值
	double Zavemin;
	double Zavemax;
	double Zsigma;
	int Zcount;

	int rateorppp;
	int KLOorGIM;
	int intervset;
	int schemeset;
	int cccount;
	int frennn;
	double pppsrate;

	ppp_option_t* next;

	int IFB_num;
	int CUR_mum;
	int PRE_num;
	int ifb_mode;
};

struct rel_option_t
{
	//
};

struct net_option_t
{
	//
};

/// Multi-GNSS PPP Configuration
class ppp_config
{
public:
	ppp_config();
	~ppp_config();

	unsigned int StaNum;                 // Station number
	unsigned int batchmode;              // batch process mode

	ppp_option_t* popt;                  // PPP option

	int ReadConfig(const char* cfg_path);                          // Read PPP Configuration
	int WriteConfig(ppp_option_t* popt, char* cfg_path);     // Write PPP Configuration

private:
	int add_opt(ppp_option_t* p, ppp_option_t** tail);
};
#endif
