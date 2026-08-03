/* -------------------------------------------------------------------------
* SUPREME_ParaEst.cpp : GNSS Parameter Estimation
*
*           Copyright (C) by C. Zhao  rights reserved
*           Contact caszcb@163.com
*
* Create : 2019.07.30
* ------------------------------------------------------------------------- */
#ifndef SUPREME_PARAMETER_ESTIMATE_H_HH
#define SUPREME_PARAMETER_ESTIMATE_H_HH
#include <vector>
#include <deque>
#include "SUPREME_Options.h"
#include "SUPREME_Matrix.h"
#include "SUPREME_ObsRNX.h"
#include "SUPREME_OrbClk.h"
#include "SUPREME_GenData.h"
#include "SUPREME_DataProcess.h"
#include <algorithm>
#include <set>
#include <iomanip>

using namespace std;

#define VAR_POS     SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_CLK     SQR(100000.0)       /* init variance receiver clock (m^2) */
#define VAR_ISB		SQR(100000.0)				/* init variance  ISB (m^2) */
#define VAR_ZTD     SQR( 0.3)         /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define ThresMW     5.0
#define THresGF		0.15	

typedef struct ION_x
{
	vector <int>time;

	double a;
	double b;
	double c;


	int num;

	vector<double>dgf_value;
	vector<double>gf_value;
	vector<double>ele;
};

class LeastSquarePPP
{
public:
	LeastSquarePPP();
	LeastSquarePPP(ppp_option_t& popt);
	~LeastSquarePPP();

	unsigned int _traceMode;
	unsigned int _OutputFile;
	vector <int>num_zhy; vector <int>num_zhy_pre; vector <int>num_zhy_accu; vector <int>num_zhy_accu_pre; vector<vector<string>> type_pre; vector<vector<string>> type_;
	vector<vector<int>>type_int; vector<vector<int>>type_int_pre;
	unsigned int m_SatValidN;                  // Valid satellite number
	double m_Sta[3];                           // Station position
	double m_TropZenith;                       // Zenith troposphere delay
	double m_TropGradientN, m_TropGradientE;   // Troposphere gradient

	int gf_count[MAX_GNSS_PRN_NUM];
	double gf_cur[MAX_GNSS_PRN_NUM];
	double gf_pre[MAX_GNSS_PRN_NUM];

	double gf_cur13[MAX_GNSS_PRN_NUM];
	double gf_pre13[MAX_GNSS_PRN_NUM];

	vector<unsigned int> prnnameL, prnnametotal;
	vector<vector<double>> L_Cur;


	void ION_Initialize(ION_x& ion);
	void ION_del_one_first(ION_x& ion);

	int mw_count[MAX_GNSS_PRN_NUM];
	int mw_count13[MAX_GNSS_PRN_NUM];
	double mw_cur[MAX_GNSS_PRN_NUM];
	double mw_pre[MAX_GNSS_PRN_NUM];
	double mw_avg[MAX_GNSS_PRN_NUM];

	double mw_cur13[MAX_GNSS_PRN_NUM];
	double mw_pre13[MAX_GNSS_PRN_NUM];
	double mw_avg13[MAX_GNSS_PRN_NUM];

	double mw_var[MAX_GNSS_PRN_NUM];

	int iCycle[MAX_GNSS_PRN_NUM];
	gnsstime tLast[MAX_GNSS_PRN_NUM];
	vector<unsigned int> m_Prn;                // Satellite prn
	vector<int> m_Slot;                        // Glonass slot number
	SatInfo m_SatInfo;                         // Satellite information

	unsigned int flag_epoch_count;         // Process epoch count

	void calculate_dion(unsigned int satno);
	double get_ion_xishu(int type, unsigned int prn, NavData* m_NavData);
	double get_wth(unsigned int prn, int flag, NavData* m_NavData);

	// ★ 修改1：Cycle_slip 新增两个参数 satclks_prn / satclks_prn_prev
	// 原来：
	//   int Cycle_slip(..., vector<double>& satclks, vector<double>& satclks_prev);
	// 改为：
	int Cycle_slip(ppp_option_t& popt, gnsstime& t,
		vector<zhydata>& prepppdata, vector<zhyinfo>& presatinfo,
		PPPObsData& pppdata, PreciseData& predata,
		vector<double>& ion_gimsub, vector<int>& commprn, vector<double>& isb,
		vector<double>& trop, vector<double>& trop_prev,
		vector<double>& satclks, vector<double>& satclks_prev,
		vector<unsigned int>& satclks_prn,       // ★ 新增
		vector<unsigned int>& satclks_prn_prev); // ★ 新增

	int Cycle_slipTDCP(ppp_option_t& popt, vector <zhydata>& prepppdata, vector <zhyinfo>& presatinfo, PPPObsData& pppdata, PreciseData& predata, vector<double>& ion_gimsub, vector<int>& commprn);

	void QC_Scdia(vector<double> L_P, vector<double> L_C, vector<double> P_P, vector<double> P_C, vector<int> Obstype, unsigned int prn, PreciseData& predata, zhydata data);

	int LSQ_qc(MatrixT& b, MatrixT& l, MatrixT& p, MatrixT& v, MatrixT& nn, MatrixT& x, int ind, unsigned int prn);

	int QC_zhy_CycleSlip(vector<zhydata>& data, PPPObsData& pppdata);

	int QC_detect(ppp_option_t& popt, vector <zhydata>& prepppdata, PPPObsData& pppdata, PreciseData& predata, GenData& gendata);

	int detslp_MW(ppp_option_t& popt, PPPObsData& pppdata, PreciseData& predata, GenData& gendata, vector<vector<double>> L_Cur, vector<vector<double>>& mwNvalue, int frenummwgf);

	int detslp_GF(ppp_option_t& popt, PPPObsData& pppdata, PreciseData& predata, GenData& gendata, vector<vector<double>> L_Cur, vector<vector<double>>& gfNvalue, int frenummwgf);
	int detslp_GF13(ppp_option_t& popt, PPPObsData& pppdata, PreciseData& predata, GenData& gendata);
	int Lsq_Option_Init(ppp_option_t& popt);
	double setvaluecal(int eopchtime, int satnum);

	int SetISBindex(PPPObsData& inputdata);

	int Lsq_ParaInit(ppp_option_t& popt, unsigned int SatN);

	int KF_ParaInt(ppp_option_t& popt, unsigned int SatN);

	int Lsq_CalMatrix(ppp_option_t& popt, GenData& gendata, PPPObsData& pppdata, PreciseData& predata, ClockData& clkdata, int flag_zhy);

	int Lsq_LSMatrix();

	int Lsq_FiltingMatrix();

	int Lsq_StateUpdate(ppp_option_t& popt);
	int Lsq_StateUpdateSPP(ppp_option_t& popt);

	int Lsq_Filter(ppp_option_t& popt, gnsstime& t, GenData& gendata, ObsEpochData& obsdata, PreciseData* predata, ClockData* clkdata, vector<zhydata>& prepppdata, vector<zhyinfo>& presatinfo, int zhy_flag, deque<vector<int>>& compareprn, deque<vector<double>>& iongim123);

	int Lsq_Quality_Check(ppp_option_t& popt, PPPObsData& inputdata);

	// Results Output
	int PPP_CRD_Trop_Output(ppp_option_t& popt);

	int PPP_Res_Output(ppp_option_t& popt);

	int PPP_STEC_Sat_Output(ppp_option_t& popt);

	int PPP_STEC_Output(ppp_option_t& popt, PPPObsData& inputdata);

	int LsqSPP_Output(ppp_option_t& popt);
	int LsqPPP_Output(ppp_option_t& popt, PPPObsData& pppdata);
	//*********************************************

	int UdState2KF(ppp_option_t* popt, unsigned int SatN, PPPObsData& pppdata, vector<double>ion_GIM, PreciseData& predata);

	int Udpos(ppp_option_t* popt);

	int udclk(ppp_option_t* popt);

	int udion(ppp_option_t* popt, unsigned int satN, vector<int>& prn_);

	int udtrop(ppp_option_t* popt);

	int udbias(ppp_option_t* popt, unsigned int satN_, vector<int>& prn_);

	int udisb(ppp_option_t* popt);

	int udifb(ppp_option_t* popt, vector<int>& prn_);

	int Residual_editing(MatrixT& v, MatrixT& r, PPPObsData& pppdata);

	int DelSat_ElevThresold(ppp_option_t& popt, SatInfo& satinfo, PPPObsData& pppdata);

	void DeleteSatData_num(PPPObsData& input, unsigned int index);

	void dele_cur(zhydata& cur, unsigned int index);

	void cal_IFB_num(ppp_option_t& popt, PPPObsData& pppdata);

	//############################################
	int Initialize();
	int Clear();
private:

	double gf_ion[MAX_GNSS_PRN_NUM];
	double gf_sigemaIon[MAX_GNSS_PRN_NUM];

	std::vector<double> ion;       // 当前历元的电离层延迟
	std::vector<double> ion_prev;  // 前一历元的电离层延迟

	std::vector<double> trop;      // 当前历元对流层延迟（备用，实际由 Lsq_Filter 局部变量传入）
	std::vector<double> trop_prev; // 前一历元对流层延迟（同上）

	// ★ 修改2：删除原有 satclk / satclk_prev / prn_prev，
	//    替换为命名一致的 satclks / satclks_prev，
	//    并新增 satclks_prn / satclks_prn_prev 用于 PRN 索引查找
	//
	// 原来（删除）：
	//   std::vector<double> satclk_prev;
	//   std::vector<double> satclk;
	//   std::vector<unsigned int> prn_prev;
	//
	// 改为：
	std::vector<double>       satclks;              // 当前历元卫星钟差（顺序与 satclks_prn 对应）
	std::vector<double>       satclks_prev;         // 前一历元卫星钟差（顺序与 satclks_prn_prev 对应）
	std::vector<unsigned int> satclks_prn;          // ★ 新增：当前历元 satclks 对应的 PRN 顺序
	std::vector<unsigned int> satclks_prn_prev;     // ★ 新增：前一历元 satclks 对应的 PRN 顺序

	ION_x Ion_X[MAX_GNSS_PRN_NUM];

	double sigema30[MAX_GNSS_PRN_NUM];

	double sigema60[MAX_GNSS_PRN_NUM];

	double PSD[MAX_GNSS_PRN_NUM];


	char _StaName[10];
	// PPP option setting
	int _SystemIndex[255];
	int _IFB_Index[25];
	unsigned int _MathModel, _Mode, _Pmode, _FreqNum, _System, _SystemN, _FilterMode, _TimeMode, _GIM_P1P2, _IS_ISB, _IS_DCB, _ISBMODE;
	unsigned int _ISBNum, _GPS_ISB, _GLO_ISB, _GAL_ISB, _BDS_ISB, _BDS2_ISB, _BDS3_ISB;
	unsigned int _SystemN_pre, _ParaN_const_pre;
	unsigned int _IonModel, _TropModel, _RclkModel;
	unsigned int est_xyz, est_rclk, est_trop, est_ion, est_cbias, est_pbias;
	double m_SigemaP, m_SigemaL;
	double _TropProcessNoise, _TropGNnoise, _TropGEnoise, _IonProcessNoise, _RclkNoise;


	vector<double> _ResAmb;

	gnsstime _previousEpoch;

	// Sat num used in ppp
	unsigned int _GpsSatNum, _GloSatNum, _GalSatNum, _BdsSatNum;


	vector<unsigned int> _Prn_previous, _Prn_common;    // Prn in 1st epoch; Common Prn in both 1st and 2nd epoch
	unsigned int _PsuedoRclkNum, _AmbNum, _IonNum, _TropNum, _ObsEqNum;
	unsigned int XYZ_Rclk_Trop_Num = 0, XYZ_Rclk_Num = 0, XYZ_Num = 0; // Parameter number in filter
	unsigned int _CommSatN, _PsuedoObsNum;
	unsigned int _ParaN_const, _ParaN_sat, _ParaN_EstTotal, _ParaN_iondcb, _PsuedoObsN;
	unsigned int _CbiasNum;

	unsigned int _Index_Cbias, _Index_Pbias;
	unsigned int _RefSatPrn, _RefSatValidFlag;
	int _RefSatIndex;

	double Antdiff[3];
	double CoarseXYZ[3];
	double CoarseXYZ_pre[3];
	double m_Sta_tmp[3];
	double m_Rclk;

	double m_dt;

	double _TideCorrection[3];

	MatrixT B, P, Q, L;
	MatrixT B_LS, P_LS, L_LS, Nbb, V, X;
	MatrixT B_x, P_x, L_x, _pp, _xx;
	MatrixT X_zhy, H_zhy, V_zhy, R_zhy, P_zhy;
	MatrixT X_Trans, P_Trans;

	vector<double> BIAS_;
	vector<double> BIAS_1;
	vector<double> BIAS_2;
	vector<double> BIAS_3;
	vector<double> BIAS_4;
	vector<double> ION_;
	double XIgema0;

	ObsFileHeader _ObsHeader;
	vector<unsigned int> _Cycle_Slip;
	vector<unsigned int> _Cycle_Sliptest;
	vector<unsigned int> _Cycle_Sliptest13;
	vector<unsigned int> _Cycle_Sliptestgf;
	vector<unsigned int> _Cycle_Sliptestgf13;

	int LsqEpochClear();
};
class SUB
{
public:
	SUB();
	~SUB();

	deque<vector<int>> compareprn;
	deque<vector<double>> iongim123;
	deque<vector<double>> isb123;
};

void Rinexouthead(ppp_option_t& popt, ofstream& rinexo);

void Rinexoutbody(ofstream& rinexo, CodeDCB& dcbdata, ppp_option_t& popt, vector <zhydata>& prepppdata);

void Rinexoutbodypre(ofstream& rinexo, CodeDCB& dcbdata, ppp_option_t& popt, vector <zhydata>& prepppdata);

void Rinexoutbodypsim(ofstream& rinexo, CodeDCB& dcbdata, ppp_option_t& popt, vector <zhydata>& prepppdata);

int LeastSquare(MatrixT& b, MatrixT& l, MatrixT& p, MatrixT& v, MatrixT& nn, MatrixT& x);

int LeastSquareQcc(MatrixT& b, MatrixT& l, MatrixT& p, MatrixT& v, MatrixT& nn, MatrixT& cc, MatrixT& x);

int LeastSquareconstraint(MatrixT& b, MatrixT& l, MatrixT& p, MatrixT& c, MatrixT& w, MatrixT& v, MatrixT& vz, MatrixT& nn, MatrixT& cc, MatrixT& x, MatrixT& Pxx);

double Kfilter(MatrixT& X, MatrixT& H, MatrixT& Pp, MatrixT& R, MatrixT& V, MatrixT& X_Trans, MatrixT& P_Trans, unsigned int flag_epoch_count);

void deletecurrentobsdata(vector<zhydata>& data, int index);


#endif