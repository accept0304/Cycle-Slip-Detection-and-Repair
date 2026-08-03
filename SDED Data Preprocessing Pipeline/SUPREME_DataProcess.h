#ifndef SUPREME_DATA_PROCESSING_H_HH
#define SUPREME_DATA_PROCESSING_H_HH

#include "SUPREME_Options.h"
#include "SUPREME_GnssTime.h"
#include "SUPREME_GenData.h"
#include "SUPREME_NavRNX.h"
#include "SUPREME_ObsRNX.h"
#include "SUPREME_ClkRNX.h"
#include <vector>
using namespace std;
/// SPP output parameters
class SppOutputPara
{
public:
	SppOutputPara();
	SppOutputPara(unsigned int satnum);
	~SppOutputPara();

	unsigned int m_Size;
	double *Spos, *Sclk, *Trop, *Elev, *Azim, *Resp;
	double StaPos_spp[3], StaPosCov_spp[3];
	double Rclk_spp, RclkCov_spp;

	void DeleSat(unsigned int i);
};

/// Precise Point Positioning Process Data
class PPPObsData:public CodeDCB
{
public:
	PPPObsData();
	PPPObsData(unsigned int satn);
	~PPPObsData();

	double ThresGF, ThresMW;                                    // Circle slip detect threshold

	double Ion_GPSA[4], Ion_GPSB[4];                            // Klobuchar parameter

	unsigned int m_SatCount;                                    // Satellite number
	unsigned int *m_Valid_flag;                                 // Satellite valid flag: 1->useful 0-> unuseful
	unsigned int *m_Prn;                                        // Satellite prn
	double* m_Pr1, * m_Pr2, * m_Pr3, * m_Pr4, * m_Pr5,* m_L1, * m_L2, * m_L3, * m_L4, * m_L5,* m_D1;         // Pusedorange and carrier phase

	int *_num; int cur_num; int *_num_accu;

	vector<vector<string>>Prn_Otype;
	vector<vector<int>>Prn_Otype_int;
	vector<vector<string>>Prn_Otype_channel;
	vector<string> epoch_otype;
	
	int Fcb_num;

	gnsstime m_obst;

	int Initialize(unsigned int satn);                          // Initialize by satllite number
	int CheckObsData(ppp_option_t& popt, GenData& gendata, unsigned int satno, gnsstime& t, PreciseData *predata, ClockData *clkdata);

	int GetBlockData(ppp_option_t& popt, GenData *gendata, ObsEpochData& obsdata, PreciseData *predata, ClockData *clkdata, int flag);

	int Clear();
	int P_detec(GenData* gendata, ppp_option_t& popt, unsigned int satno);
private:
	int fl, pflag,f2,pflag2, f3,f4,f5,pflag3, pflag4, pflag5;
	int obstype_num;
	double c1, p1, p2, p3,p4,p5, l1, l2, l3,l4,l5,d1,d2,d3,d4,d5;
	char obs_type[25];
	int GetGPSData(CodeDCB &dcbdata, ObsEpochData& obsdata, SatelliteData &epochdata, unsigned int freqn, int flag,int mode);
	int GetGLOData(CodeDCB &dcbdata, ObsEpochData& obsdata, SatelliteData &epochdata, unsigned int freqn, int flag, int mode);
	int GetGALData(CodeDCB &dcbdata, ObsEpochData& obsdata, SatelliteData &epochdata, unsigned int freqn, int flag, int mode);
	int GetBDSData(CodeDCB &dcbdata, ObsEpochData& obsdata, SatelliteData &epochdata, unsigned int freqn, int flag,int prn, int mode);
};

void SetCodePri(int sys, int freq, const char *pri);

int RankDetect(ppp_option_t& popt, unsigned int satn);

void DeleteSatData(PPPObsData& input, unsigned int index);

typedef struct zhydata
{
	double  m_L1[200];
	double  m_L2[200];
	double  m_L3[200];
	double  m_L4[200];
	double  m_L5[200];
	double  m_P1[200];
	double  m_P2[200];
	double  m_P3[200];
	double  m_P4[200];
	double  m_P5[200];
	unsigned int m_SatCount;
	unsigned int m_Prn[200];
}zhydata;
#endif