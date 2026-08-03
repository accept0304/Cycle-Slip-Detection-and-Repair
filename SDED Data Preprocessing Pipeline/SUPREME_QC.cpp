/* -------------------------------------------------------------------------
* SUPREME_QC.cpp : Quality Check and Control
*
*           Copyright (C) by C. Zhao 2016-2018 All rights reserved
*           Contact zhaoC.@whigg.ac.cn
*
* Create : 2018.03.16
* ------------------------------------------------------------------------- */
#include <math.h>

#include "SUPREME_QC.h"
#include "SUPREME_CommonFunction.h"
#include "SUPREME_Options.h"
#include "SUPREME_ParaEst.h"

#define MAX_SLIP_GAP            120
#define CYCLE_SLIP_VAR          10
#define THRES_MW_JUMP           10

/// psuedorange quality check
#define PSEUDORANGE_C1P1_K1         10
#define PSEUDORANGE_P1P2_K2         30

//static int gf_count[MAX_GNSS_PRN_NUM] = { 0 };
//static double gf_tow[MAX_GNSS_PRN_NUM] = { 0.0 };
//static double gf_old[MAX_GNSS_PRN_NUM] = { 0.0 };
//
//static int mw_count[MAX_GNSS_PRN_NUM] = { 0 };
//static double mw_tow[MAX_GNSS_PRN_NUM] = { 0.0 };
//static double mw_avg[MAX_GNSS_PRN_NUM] = { 0.0 };
//static double mw_var[MAX_GNSS_PRN_NUM] = { 0.0 };
//
//
//void CycleSlipInit()
//{
//	memset(gf_count, 0, sizeof(int)*MAX_GNSS_PRN_NUM);
//	memset(gf_tow, 0.0, sizeof(double)*MAX_GNSS_PRN_NUM);
//	memset(gf_old, 0.0, sizeof(double)*MAX_GNSS_PRN_NUM);
//
//	memset(mw_count, 0, sizeof(int)*MAX_GNSS_PRN_NUM);
//	memset(mw_tow, 0.0, sizeof(double)*MAX_GNSS_PRN_NUM);
//	memset(mw_avg, 0.0, sizeof(double)*MAX_GNSS_PRN_NUM);
//	memset(mw_var, 0.0, sizeof(double)*MAX_GNSS_PRN_NUM);
//}

/// Receiver Clock Slip Detect and Repair
int RecClkSlipRepair(ppp_option_t &popt, unsigned int Satn, unsigned int *Prn, double *Pr1, double *Pr2, double *L1, double *L2)
{
	unsigned int FreqNum = 0;
	int i, sat, orb_n = 0, validGps, cjGps;
	int bObserved[GNSS_SATNO_NUM];
	double delta0 = 0.0, delta1 = 0.0, d1, d2, d3, d4, ddd1, ddd2;
	double lam[2] = { 0 };
	double CJ_F1, CJ_F2;

	FreqNum = popt.freqn;

	for (i = 0; i<GNSS_SATNO_NUM; i++) bObserved[i] = 0;

	validGps = cjGps = 0;

	for (i = 0; i<Satn; i++)
	{
		sat = Prn[i];
		if (popt.math_model == 0 || popt.math_model == 4) 
		{
			if (GetSystem(sat) == 'G' || 'E' ) 
			{
				lam[0] = Get_WaveLength(sat, 1, orb_n);//L1
				lam[1] = Get_WaveLength(sat, 2, orb_n);//L2
			}
			else if (GetSystem(sat) == 'C') 
			{
				lam[0] = Get_WaveLength(sat, 1, orb_n);//B1
				lam[1] = Get_WaveLength(sat, 3, orb_n);//B3
			}
			else if (GetSystem(sat) == 'R')
			{
				lam[0] = Get_WaveLength(sat, 1, orb_n);//G1    //但是频分多址需要传递 进来广播星历参参数
				lam[1] = Get_WaveLength(sat, 3, orb_n);//G2
			}

		}
		else
		{
			lam[0] = Get_WaveLength(sat, 1, orb_n);
			lam[1] = Get_WaveLength(sat, 2, orb_n);

		}
		
		//if (sat>GPS_SATNUM) continue;

		if (Pr1[i] * L1[i] == 0.0) continue;
		if ((popt.math_model == 0||popt.math_model==4) && (Pr2[i] * L2[i]) == 0.0) continue;

		if (popt.obs0[sat - 1][0] * popt.obs0[sat - 1][2] == 0.0)
			continue;
		if ((popt.obs0[sat - 1][1] * popt.obs0[sat - 1][3] == 0.0) && FreqNum == 2)
			continue;

		validGps++;

		d1 = Pr1[i] - popt.obs0[sat - 1][0];            // P1
		d2 = Pr2[i] - popt.obs0[sat - 1][1];            // P2
		d3 = (L1[i] - popt.obs0[sat - 1][2])*lam[0];    // L1
		d4 = (L2[i] - popt.obs0[sat - 1][3])*lam[1];    // L2

		if (fabs(d1 - d3)>290000)   //ms clock jump
		{
			delta0 += d1 - d3;
			delta1 += d2 - d4;
			cjGps++;
		}
	}

	// Receiver clock slip detect
	if (cjGps != 0 && cjGps == validGps)
	{
		d1 = delta0 / cjGps;
		d2 = delta1 / cjGps;

		CJ_F1 = 0.0;   //flag for clock jump
		CJ_F2 = 0.0;
		CJ_F1 = d1 / LIGHTSPEED * 1000.0;
		CJ_F2 = round(CJ_F1);

		if ((fabs(CJ_F1 - CJ_F2)<2.5E-2) && (FreqNum == 2))
		{
			popt.clkJump += (int)CJ_F2;
			printf("*** WARNING: clock jump=%d(ms)\n", popt.clkJump);

			FILE *trace_fp = fopen(popt.IOfile.out_f.out_trace_fpath, "a+");
			if (trace_fp)
			{
				fprintf(trace_fp, "  WARNING: Receiver Clock Jump=%d(ms)\n", popt.clkJump);
			}
			if (trace_fp) fclose(trace_fp);
		}
		else if (FreqNum == 1)
		{
			popt.clkJump += (int)CJ_F2;
			printf("*** WARNING: clock jump=%d(ms)\n", popt.clkJump);

			FILE *trace_fp = fopen(popt.IOfile.out_f.out_trace_fpath, "a+");
			if (trace_fp)
			{
				fprintf(trace_fp, "  WARNING: Receiver Clock Jump=%d(ms)\n", popt.clkJump);
			}
			if (trace_fp) fclose(trace_fp);
		}
	}

	// Receiver clock slip repaire
	for (i = 0; i<Satn; i++)
	{
		sat = Prn[i];

		//if (sat>GPS_SATNUM) continue;

		bObserved[sat - 1] = 1;

		popt.obs0[sat - 1][0] = Pr1[i];
		popt.obs0[sat - 1][1] = Pr2[i];
		popt.obs0[sat - 1][2] = L1[i];
		popt.obs0[sat - 1][3] = L2[i];

		ddd1 = popt.clkJump*LIGHTSPEED / 1000.0;
		ddd2 = popt.clkJump*LIGHTSPEED / 1000.0;

		//repair for phase observations
		if (L1[i] != 0.0) L1[i] += ddd1 / lam[0];

		if (L2[i] != 0.0) L2[i] += ddd2 / lam[1];
	}

	for (i = 0; i<GNSS_SATNO_NUM; i++)
	{
		if (bObserved[i] == 0)
			popt.obs0[i][0] = popt.obs0[i][1] = popt.obs0[i][2] = popt.obs0[i][3] = 0.0;
	}

	return 1;
}

int RC_Repair(ppp_option_t& popt, PPPObsData& inputdata, NavData* navdata, double* Ele, unsigned int flag_epoch_count ) 
{
	unsigned int FreqNum = 0;
	int i, sat, orb_n = 0, validGps, cjGps;
	int bObserved[GNSS_SATNO_NUM];
	double delta0 = 0.0, delta1 = 0.0, d1, d2, d3, d4, ddd1, ddd2;
	double lam[5] = { 0 };
	double CJ_F1, CJ_F2;
	int Satn = inputdata.m_SatCount;
	FreqNum = popt.freqn;
	double g1 = 0.0, g0 = 0.0; double ele = 0.0; double KEY = 0.0,ThresGF=0.15;//gamp的阈值
	for (i = 0; i < GNSS_SATNO_NUM; i++) bObserved[i] = 0;

	validGps = cjGps = 0;

	for (i = 0; i < Satn; i++)
	{
		sat = inputdata.m_Prn[i];
		if (popt.math_model == 0 || popt.math_model == 4)
		{
			if (GetSystem(sat) == 'G')
			{
				lam[0] = Get_WaveLength(sat, 1, orb_n);//L1
				lam[1] = Get_WaveLength(sat, 2, orb_n);//L2
				lam[2] = Get_WaveLength(sat, 3, orb_n);//L3
			}
			else if (GetSystem(sat) == 'E') 
			{
				lam[0] = Get_WaveLength(sat, 1, orb_n);//L1
				lam[1] = Get_WaveLength(sat, 2, orb_n);//L2
				lam[2] = Get_WaveLength(sat, 3, orb_n);//L1
				lam[3] = Get_WaveLength(sat, 4, orb_n);//L2
				lam[4] = Get_WaveLength(sat, 5, orb_n);//L1
			}
			else if (GetSystem(sat) == 'C')
			{
				lam[0] = Get_WaveLength(sat, 1, orb_n);//L1
				lam[1] = Get_WaveLength(sat, 3, orb_n);//L2
				lam[2] = Get_WaveLength(sat, 2, orb_n);//L1
				lam[3] = Get_WaveLength(sat, 4, orb_n);//L2
				lam[4] = Get_WaveLength(sat, 5, orb_n);//L1
			}
			else if (GetSystem(sat) == 'R')
			{
				orb_n= Get_Nav_Orbitn(sat, navdata);
				lam[0] = Get_WaveLength(sat, 1, orb_n);//L1
				lam[1] = Get_WaveLength(sat, 2, orb_n);//L2
				lam[2] = Get_WaveLength(sat, 3, orb_n);//L1
				lam[3] = Get_WaveLength(sat, 4, orb_n);//L2
				lam[4] = Get_WaveLength(sat, 5, orb_n);//L1
			}

		}
		else
		{
			lam[0] = Get_WaveLength(sat, 1, orb_n);
			lam[1] = Get_WaveLength(sat, 2, orb_n);

		}

		//if (sat>GPS_SATNUM) continue;

		if (inputdata.m_Pr1[i] * inputdata.m_L1[i] == 0.0) continue;
		if ((popt.math_model == 0 || popt.math_model == 4) && (inputdata.m_Pr2[i] * inputdata.m_L2[i] == 0.0)) continue;

		if (popt.obs0[sat - 1][0] * popt.obs0[sat - 1][2] == 0.0)
			continue;
		if ((popt.obs0[sat - 1][1] * popt.obs0[sat - 1][3] == 0.0) && (popt.math_model == 0 || popt.math_model == 4))
			continue;
		d1 = inputdata.m_Pr1[i] - popt.obs0[sat - 1][0];            // P1
		d2 = inputdata.m_Pr2[i] - popt.obs0[sat - 1][1];            // P2
		d3 = (inputdata.m_L1[i] - popt.obs0[sat - 1][2]) * lam[0];    // L1
		d4 = (inputdata.m_L2[i] - popt.obs0[sat - 1][3]) * lam[1];    // L2
		if ((popt.math_model == 0 || popt.math_model == 4))   // 历元数量  在加进来一个参数 
		{
			ele = Ele[i]*180/PI;
			g1 = inputdata.m_L1[i] * lam[0] - inputdata.m_L2[i] * lam[1];//gf
			if (ele < popt.ele_mask) { ele = popt.ele_mask; }
			if (ele >= 15)KEY = ThresGF;
			else { KEY = -ThresGF / 15 * ele + 2 * ThresGF; } 
			g0 = popt.gf_pre[sat - 1];
			popt.gf_pre[sat - 1] = g1;
			if (g0 != 0 && fabs(g1 - g0) > min(KEY * (flag_epoch_count-popt.Time_dif[sat-1]), 1.5))continue;
		}
		validGps++;
		
		if (fabs(d1 - d3) > 290000)   //ms clock jump   dP-dL//
		{
			delta0 += d1 - d3;
			delta1 += d2 - d4;
			cjGps++;
		}
	}

	// Receiver clock slip detect
	if (cjGps != 0 && cjGps == validGps)
	{
		d1 = delta0 / cjGps;
		d2 = delta1 / cjGps;

		CJ_F1 = 0.0;   //flag for clock jump
		CJ_F2 = 0.0;
		CJ_F1 = d1 / LIGHTSPEED * 1000.0;  //M  传化成毫秒  且钟跳是 有整数特性
		CJ_F2 = round(CJ_F1);
		//CJ_F2 = d2 / LIGHTSPEED * 1000.0;
		if ((fabs(CJ_F1 - CJ_F2) < 2.5E-2) && (popt.math_model == 0 || popt.math_model == 4))
		{
			popt.clkJump += (int)CJ_F2;
			printf("*** WARNING: clock jump=%d(ms)\n", popt.clkJump);

			FILE* trace_fp = fopen(popt.IOfile.out_f.out_trace_fpath, "a+");
			if (trace_fp)
			{
				fprintf(trace_fp, "  WARNING: Receiver Clock Jump=%d(ms)\n", popt.clkJump);
			}
			if (trace_fp) fclose(trace_fp);
		}
		else  
		{
			popt.clkJump += (int)CJ_F2;
			printf("*** WARNING: clock jump=%d(ms)\n", popt.clkJump);

			FILE* trace_fp = fopen(popt.IOfile.out_f.out_trace_fpath, "a+");
			if (trace_fp)
			{
				fprintf(trace_fp, "  WARNING: Receiver Clock Jump=%d(ms)\n", popt.clkJump);
			}
			if (trace_fp) fclose(trace_fp);
		}
	}

	// Receiver clock slip repaire
	for (i = 0; i < Satn; i++)
	{
		sat = inputdata.m_Prn[i];
		
		if (popt.math_model == 0 || popt.math_model == 4)
		{
			if (GetSystem(sat) == 'G')
			{
				lam[0] = Get_WaveLength(sat, 1, orb_n);//L1
				lam[1] = Get_WaveLength(sat, 2, orb_n);//L2
				lam[2] = Get_WaveLength(sat, 3, orb_n);//L3
			}
			else if (GetSystem(sat) == 'E')
			{
				lam[0] = Get_WaveLength(sat, 1, orb_n);//L1
				lam[1] = Get_WaveLength(sat, 2, orb_n);//L2
				lam[2] = Get_WaveLength(sat, 3, orb_n);//L1
				lam[3] = Get_WaveLength(sat, 4, orb_n);//L2
				lam[4] = Get_WaveLength(sat, 5, orb_n);//L1
			}
			else if (GetSystem(sat) == 'C')
			{
				lam[0] = Get_WaveLength(sat, 1, orb_n);//L1
				lam[1] = Get_WaveLength(sat, 2, orb_n);//L2
				lam[2] = Get_WaveLength(sat, 3, orb_n);//L1
				lam[3] = Get_WaveLength(sat, 4, orb_n);//L2
				lam[4] = Get_WaveLength(sat, 5, orb_n);//L1
			}
			else if (GetSystem(sat) == 'R')
			{
				orb_n = Get_Nav_Orbitn(sat, navdata);
				lam[0] = Get_WaveLength(sat, 1, orb_n);//L1
				lam[1] = Get_WaveLength(sat, 2, orb_n);//L2
				lam[2] = Get_WaveLength(sat, 3, orb_n);//L1
				lam[3] = Get_WaveLength(sat, 4, orb_n);//L2
				lam[4] = Get_WaveLength(sat, 5, orb_n);//L1
			}
		}
		else
		{
			lam[0] = Get_WaveLength(sat, 1, orb_n);
			lam[1] = Get_WaveLength(sat, 2, orb_n);

		}
		//if (sat>GPS_SATNUM) continue;

		bObserved[sat - 1] = 1;

		popt.obs0[sat - 1][0] = inputdata.m_Pr1[i];
		popt.obs0[sat - 1][1] = inputdata.m_Pr2[i];
		popt.obs0[sat - 1][2] = inputdata.m_L1[i];
		popt.obs0[sat - 1][3] = inputdata.m_L2[i];
		popt.Time_dif[sat - 1] = flag_epoch_count;
		

		ddd1 = popt.clkJump * LIGHTSPEED / 1000.0;
		ddd2 = popt.clkJump * LIGHTSPEED / 1000.0;

		//repair for phase observations
		if (inputdata.m_L1[i] != 0.0) inputdata.m_L1[i] += ddd1 / lam[0];

		if (inputdata.m_L2[i] != 0.0) inputdata.m_L2[i] += ddd2 / lam[1];

		if (inputdata.m_L3[i] != 0.0) inputdata.m_L3[i] += ddd2 / lam[2];

		if (inputdata.m_L4[i] != 0.0) inputdata.m_L4[i] += ddd2 / lam[3];

		if (inputdata.m_L5[i] != 0.0) inputdata.m_L5[i] += ddd2 / lam[4];
	}

	for (i = 0; i < GNSS_SATNO_NUM; i++)
	{
		if (bObserved[i] == 0)
			popt.obs0[i][0] = popt.obs0[i][1] = popt.obs0[i][2] = popt.obs0[i][3] = 0.0;
	}

	return 1;

	
}



/// Quality check of cycle slip and reciver clock slip
int QC_CycleSlip_ClockSlip(ppp_option_t &popt, PPPObsData &inputdata, NavData *navdata)
{
	int i = 0, j = 0, satn = inputdata.m_SatCount;
	int valid_satn = 0;
	double tow = popt.CurrentEpoch.m_SecondofWeek;

	unsigned int FreqN = popt.freqn;

	// receiver clock slip detect and repair.
	RecClkSlipRepair(popt, inputdata.m_SatCount, inputdata.m_Prn, inputdata.m_Pr1, inputdata.m_Pr2, inputdata.m_L1, inputdata.m_L2);

	
	// Users can add Cycle slip detection method

	return 1;
}

int PseudoRange_C1P1_Detect(const double &c1, const double &p1)
{
	if (fabs(c1 - p1) > PSEUDORANGE_C1P1_K1)
		return 1;

	return 0;
}

// Quality chech of O-C
int QC_L(MatrixT L)
{
	double Threshold = 100;
	int n = L.GetRow();
	double *data = (double*)(malloc(sizeof(double)*n));
	for (unsigned int i = 0; i < n; i++)
	{
		data[i] = L(i + 1, 1);
	}

	double median = 0;
	median = GetMedian(data, n);

	for (unsigned int i = 0; i < n; i++)
	{
		double tmp = abs(L(i + 1, 1));
		if (tmp == 0) continue;
		if (abs(tmp - abs(median))>Threshold)
		{
			return (i + 1);
		}
	}

	if (data) free(data);

	return 0;
}

// Quality chech based on residuals. (Users can add their own QC methods)
int LeastSquarePPP::Lsq_Quality_Check(ppp_option_t &popt, PPPObsData &inputdata)
{
	int re = 0;

	// get common satellite
	vector<unsigned int> Prn_index_previous, Prn_index_current;
	_Prn_common.clear();
	for (int i = 0; i < m_SatValidN; ++i)
	{
		for (int j = 0; j < _Prn_previous.size(); ++j)
		{
			if (m_Prn[i] == _Prn_previous[j])
			{
				int k = 0;
				for (k = 0; k < _Cycle_Slip.size(); ++k)
				{
					if (m_Prn[i] == _Cycle_Slip[k])
						break;
				}
				if (k == _Cycle_Slip.size())
				{
					_Prn_common.push_back(m_Prn[i]);
					Prn_index_previous.push_back(j);
					Prn_index_current.push_back(i);
				}
			}
		}
	}

	vector<double> res_p, res_c, res, res_amb;

	int index_amb = -1;
	if (_FreqNum == 1)
	{
		_ResAmb.clear();
		for (int i = 0; i < _Prn_common.size(); ++i)
		{
			_ResAmb.push_back(V(est_xyz * 3 + est_trop + i + 1, 1));
			res_amb.push_back(fabs(V(est_xyz * 3 + est_trop + i + 1, 1)));
		}

		index_amb = MaxVectorIndex(res_amb);
	}

	for (int i = 0; i < m_SatValidN; ++i)
	{
		res_c.push_back(fabs(V(_PsuedoObsN + _PsuedoObsNum + _ObsEqNum*i + 1, 1)));
		res_p.push_back(fabs(V(_PsuedoObsN + _PsuedoObsNum + _ObsEqNum*i + 2, 1)));
	}

	int index_c = MaxVectorIndex(res_c);
	int index_p = MaxVectorIndex(res_p);

	int index_c_sat = 0, index_p_sat = 0;
	if (_FreqNum == 1)
	{
		index_c_sat = index_c;
		index_p_sat = index_p;
	}

	double thres = 0;

	if (_FreqNum == 1)
		thres = 5;

	if (res_c[index_c] > thres)
	{
		DeleteSatData(inputdata, index_c_sat);

		re = 2;

		return re;
	}

	//thres = VarErr(popt, m_Prn[index_p_sat], popt.System, m_SatInfo.m_Ele[index_p_sat], 1, 1);
	//thres = sqrt(thres) * 3;

	if (_FreqNum == 1)
		thres = 0.15;

	if ((res_p[index_p] > thres))
	{
		DeleteSatData(inputdata, index_p_sat);
		re = 2;
	}

	if (re == 2)
	{
		return re;
	}

	return re;
}

