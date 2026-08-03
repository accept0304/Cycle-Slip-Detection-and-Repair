#include"lambda.h"
#include "SUPREME_CommonFunction.h"
#include "SUPREME_Coordinate.h"
#include "SUPREME_Corrections.h"
#include "SUPREME_DataProcess.h"
#include "SUPREME_GraphPlot.h"
#include "SUPREME_Ionosphere.h"
#include "SUPREME_OrbClk.h"
#include "SUPREME_ParaEst.h"
#include "SUPREME_QC.h"
#include "SUPREME_SunMoon.h"
#include "SUPREME_Tide.h"
#include "SUPREME_Troposphere.h"
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric>
#include <set>
#include <stdlib.h>
#include <map>
#include <deque>

/// Function: Delete satellite less elevation threshold and in earth shadow
/* Parameter:
* SatelliteAntenna        antsat             satellite information for get satellite type
* GnssKalmanPara          kalmpara           kalman parameter
* double                  sun[3]             sun position
* PPPDataPara             input              PPP data and parameter
* ObsEpochData            obsdata            observation data
* PreciseData             prenav             precise ephemeris and navigation
* ----------------------------------------------------------------------------------- */

struct EpochSolution {
	double xyz[3]; 
	double rclk;   
	gnsstime t;    
};
std::ostream& operator<<(std::ostream& os, const gnsstime& t) {
	os << t.m_Year << "-"
		<< std::setw(2) << std::setfill('0') << t.m_Month << "-"
		<< std::setw(2) << std::setfill('0') << t.m_Day << " "
		<< std::setw(2) << std::setfill('0') << t.m_Hour << ":"
		<< std::setw(2) << std::setfill('0') << t.m_Min << ":"
		<< std::fixed << std::setprecision(3) << t.m_Sec;
	return os;
}

static std::map<gnsstime, EpochSolution> epoch_solutions;
static EpochSolution current_solution;  
static EpochSolution prev_solution;    
static bool has_prev_solution = false;  
double chisqr_arr[100] =
{
	/* chi-sqr(n) (alpha=0.001) */
	10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
	31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
	46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
	61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
	74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
	88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
	101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
	113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
	126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
	138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};

void IGGIII(MatrixT& v, MatrixT& p, int satnum);

void LeastSquarePPP::DeleteSatData_num(PPPObsData& input, unsigned int index)
{
	if (input.m_SatCount == 0) return; int a = index;
	vector<int>num = num_zhy;
	vector<int>num_ = num_zhy_accu;
	vector<vector<string>>type = input.Prn_Otype;
	vector<vector<int>>typeint = input.Prn_Otype_int;
	vector<vector<string>>type_chan = input.Prn_Otype_channel;
	int endindex = input.m_SatCount - 1;
	for (int i = 0; i < num_zhy_accu.size(); i++)
	{
		if (i > index)
		{
			num_[i] = num_[i] + (num[endindex] - num[index]);
		}
	}
	type_chan[index] = type_chan[endindex];
	num[index] = num[endindex];
	type[index] = type[endindex];
	typeint[index] = typeint[endindex];
	num.pop_back(); num_.pop_back(); type.pop_back(); typeint.pop_back(); type_chan.pop_back();
	input.Prn_Otype_int = typeint;
	input.Prn_Otype = type;
	num_zhy = num;
	num_zhy_accu = num_;
	input.Prn_Otype_channel = type_chan;
	type_chan.clear();
	num.clear();
	num_.clear();
	type.clear();
	input.m_Valid_flag[index] = input.m_Valid_flag[input.m_SatCount - 1];
	input.m_Prn[index] = input.m_Prn[input.m_SatCount - 1];
	input.m_Pr1[index] = input.m_Pr1[input.m_SatCount - 1];
	input.m_Pr2[index] = input.m_Pr2[input.m_SatCount - 1];
	input.m_Pr3[index] = input.m_Pr3[input.m_SatCount - 1];
	input.m_Pr4[index] = input.m_Pr4[input.m_SatCount - 1];
	input.m_Pr5[index] = input.m_Pr5[input.m_SatCount - 1];

	input.m_L1[index] = input.m_L1[input.m_SatCount - 1];
	input.m_L2[index] = input.m_L2[input.m_SatCount - 1];
	input.m_L3[index] = input.m_L3[input.m_SatCount - 1];
	input.m_L4[index] = input.m_L4[input.m_SatCount - 1];
	input.m_L5[index] = input.m_L5[input.m_SatCount - 1];

	input.cur_num -= input._num[index];
	input._num[index] = input._num[input.m_SatCount - 1];

	input.m_Prn[input.m_SatCount - 1] = 0;
	--input.m_SatCount;
}

void LeastSquarePPP::dele_cur(zhydata& cur, unsigned int index)
{
	cur.m_Prn[index] = cur.m_Prn[cur.m_SatCount - 1];
	cur.m_P1[index] = cur.m_P1[cur.m_SatCount - 1];
	cur.m_P2[index] = cur.m_P2[cur.m_SatCount - 1];
	cur.m_P3[index] = cur.m_P3[cur.m_SatCount - 1];
	cur.m_P4[index] = cur.m_P4[cur.m_SatCount - 1];
	cur.m_P5[index] = cur.m_P5[cur.m_SatCount - 1];

	cur.m_L1[index] = cur.m_L1[cur.m_SatCount - 1];
	cur.m_L2[index] = cur.m_L2[cur.m_SatCount - 1];
	cur.m_L3[index] = cur.m_L3[cur.m_SatCount - 1];
	cur.m_L4[index] = cur.m_L4[cur.m_SatCount - 1];
	cur.m_L5[index] = cur.m_L5[cur.m_SatCount - 1];
	--cur.m_SatCount;
}

int LeastSquarePPP::DelSat_ElevThresold(ppp_option_t& popt, SatInfo& satinfo, PPPObsData& pppdata)
{

	char sys = 0;
	unsigned int m = 0, prn = 0;
	double satpos[3] = { 0.0 }, satvel[3] = { 0.0 }, satpco[3] = { 0 }; // Satellite position, velocity, phase centre offset
	double ele = 0.0, azimuth = 0.0;                                    // Satellite elevation, azimuth

	for (unsigned int i = 0; i < pppdata.m_SatCount; ++i)
	{
		if (pppdata.m_Prn[i] == 0) continue;

		sys = GetSysPrn(pppdata.m_Prn[i], prn);

		for (m = 0; m < satinfo.m_Size; ++m)
		{
			if (pppdata.m_Prn[i] == satinfo.m_Prn[m])
				break;
		}
		if (m == satinfo.m_Size)
		{
			if (popt.math_model == 4)
			{
				popt.CUR_mum = popt.CUR_mum - num_zhy[i];
				DeleteSatData_num(pppdata, i--);

			}
			else
			{
				DeleteSatData(pppdata, i--);
			}
			continue;
		}
		else
		{
			if (satinfo.flag[m] == 0)
			{
				if (popt.math_model == 4)
				{
					popt.CUR_mum = popt.CUR_mum - num_zhy[i];
					DeleteSatData_num(pppdata, i--);
				}
				else
				{
					DeleteSatData(pppdata, i--);
				}
				continue;
			}

			ele = satinfo.m_Ele[m];
			azimuth = satinfo.m_Azi[m];
			satpos[0] = satinfo.m_Xs[m * 3 + 0];
			satpos[1] = satinfo.m_Xs[m * 3 + 1];
			satpos[2] = satinfo.m_Xs[m * 3 + 2];
		}

		if (ele < ELE_THRES(popt.ele_mask))
		{
			// Log in trace file
			FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
			if (logfp)
				fprintf(logfp, "  Delete Sat for lower %d° or in Earth Shadow in PPP: %c%02d(Ele:%4.1f°) \n",
					(int)popt.ele_mask, sys, prn, ele * R2D);
			if (logfp) fclose(logfp);

			if (popt.math_model == 4)
			{
				popt.CUR_mum = popt.CUR_mum - num_zhy[i];
				DeleteSatData_num(pppdata, i--);
			}
			else
			{
				DeleteSatData(pppdata, i--);
			}
		}

	}

	return 1;
}

void Rinexouthead(ppp_option_t& popt, ofstream& rinexo)
{
	string list;
	list = "C2IC6IC7IC7QC7XC7DC7PC7ZC1DC1PC1XC5DC5PC5XL2IL6IL7IL7QL7XL7DL7PL7ZL1PL1DL1XL5PL5DL5X";

	if (0)
	{
		rinexo << popt.Recordfre << ' ';
		for (int i = 0; i < popt.Recordfre; i++)
		{
			rinexo << list[(popt.recordfre[i] - 1) * 3 + 0] << list[(popt.recordfre[i] - 1) * 3 + 1] << list[(popt.recordfre[i] - 1) * 3 + 2] << ' ';
			if (i == popt.Recordfre - 1)
			{
				rinexo << endl;
			}
		}
	}
	if (1)
	{
		rinexo << "     3.03           OBSERVATION DATA    M: Mixed            RINEX VERSION / TYPE" << endl;
		rinexo << "Geo++ RINEX Logger  Geo++               20210108 043033 UTC PGM / RUN BY / DATE" << endl;
		rinexo << "************************************************************COMMENT" << endl;
		rinexo << "This file was generated by the Geo++ RINEX Logger App       COMMENT" << endl;
		rinexo << "for Android devices (Version 2.1.6). If you encounter       COMMENT" << endl;
		rinexo << "any issues, please send an email to android@geopp.de        COMMENT" << endl;
		rinexo << "Filtering Mode: BEST                                        COMMENT" << endl;
		rinexo << "************************************************************COMMENT" << endl;
		rinexo << "Geo++                                                       MARKER NAME" << endl;
		rinexo << "SPACE BORNE                                                 MARKER TYPE" << endl;
		rinexo << "Geo++               Geo++                                   OBSERVER / AGENCY" << endl;
		rinexo << "unknown             Xiaomi              MI 8                REC # / TYPE / VERS" << endl;
		rinexo << "unknown             MI 8                                    ANT # / TYPE" << endl;
		rinexo << " -2489600.6580  4037671.5665  4249721.6504                  APPROX POSITION XYZ" << endl;
		rinexo << "        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N" << endl;
		rinexo << "G    2 C1C C5Q L1C L5Q                                      SYS / # / OBS TYPES" << endl;
		rinexo << "R    2 C1C L1C                                              SYS / # / OBS TYPES" << endl;
		rinexo << "E    6 C1C C5Q L1C L5Q                                      SYS / # / OBS TYPES" << endl;
		rinexo << "C    2 C2I L2I                                              SYS / # / OBS TYPES" << endl;
		rinexo << "J    2 C1C L1C                                              SYS / # / OBS TYPES" << endl;
		rinexo << "  2021     1     8     4    30   52     GPS                 TIME OF FIRST OBS" << endl;
		rinexo << " 24 R01  1 R02 -4 R03  5 R04  6 R05  1 R06 -4 R07  5 R08  6 GLONASS SLOT / FRQ #" << endl;
		rinexo << "    R09 -2 R10 -5 R11  0 R12 -1 R13 -2 R14 -7 R15  0 R16 -1 GLONASS SLOT / FRQ #" << endl;
		rinexo << "    R17  4 R18 -3 R19  3 R20  2 R21  4 R22 -3 R23  3 R24  2 GLONASS SLOT / FRQ #" << endl;
		rinexo << "G L1C                                                       SYS / PHASE SHIFT" << endl;
		rinexo << "G L5Q -0.25000                                              SYS / PHASE SHIFT" << endl;
		rinexo << "R L1C                                                       SYS / PHASE SHIFT" << endl;
		rinexo << "E L1B                                                       SYS / PHASE SHIFT" << endl;
		rinexo << "E L1C +0.50000                                              SYS / PHASE SHIFT" << endl;
		rinexo << "E L5Q -0.25000                                              SYS / PHASE SHIFT" << endl;
		rinexo << "C L2I                                                       SYS / PHASE SHIFT" << endl;
		rinexo << "J L1C                                                       SYS / PHASE SHIFT" << endl;
		rinexo << "J L5Q -0.25000                                              SYS / PHASE SHIFT" << endl;
		rinexo << " C1C    0.000 C1P    0.000 C2C    0.000 C2P    0.000        GLONASS COD/PHS/BIS" << endl;
		rinexo << "                                                            END OF HEADER" << endl;
	}

}

void Rinexoutbody(ofstream& rinexo, CodeDCB& dcbdata, ppp_option_t& popt, vector <zhydata>& prepppdata)
{
	vector<int>tempprn; int S0prn = 0; double DCBdis = 0;

	for (int i = 0; i < prepppdata[1].m_SatCount; i++)
	{
		for (int j = 0; j < prepppdata[0].m_SatCount; j++)
		{
			if (prepppdata[1].m_Prn[i] == prepppdata[0].m_Prn[j]) { tempprn.push_back(prepppdata[1].m_Prn[i]); break; }
		}
	}
	rinexo << "> " << popt.CurrentEpoch.m_Year << " " << setw(2) << setfill('0') << popt.CurrentEpoch.m_Month << " "
		<< setw(2) << setfill('0') << popt.CurrentEpoch.m_Day << " " << setw(2) << setfill('0') << popt.CurrentEpoch.m_Hour << " "
		<< setw(2) << setfill('0') << popt.CurrentEpoch.m_Min << " " << setw(10) << setfill('0') << fixed << setprecision(7) << round(popt.CurrentEpoch.m_Sec) << "  0" << " " << tempprn.size() << endl;
	int tempsim = 0;
	for (int i = 0; i < prepppdata[1].m_SatCount; i++)
	{
		tempsim = 0;
		for (int j = 0; j < tempprn.size(); j++)
		{
			if (prepppdata[1].m_Prn[i] == tempprn[j])
			{
				for (int k = 0; k < repairstanum; k++)
				{
					if (popt.recordrepair[k] == tempprn[j] && popt.recordrepair[k] != 0)
					{
						for (int m = 0; m < popt.simnum; m++)
						{
							if (popt.simulcyclesit[m] == tempprn[j])
							{
								switch (GetSystem(prepppdata[1].m_Prn[i]))
								{
								case 'G': { S0prn = prepppdata[1].m_Prn[i] - MIN_GPS_SATNO + 1; DCBdis = (dcbdata.Get_DCB_P1C1(prepppdata[1].m_Prn[i]) * LIGHTSPEED * 1.0e-9); break; }
								case 'R': { S0prn = prepppdata[1].m_Prn[i] - MIN_GLO_SATNO + 1; break; }
								case 'C': { S0prn = prepppdata[1].m_Prn[i] - MIN_BDS_SATNO + 1; break; }
								case 'B': { S0prn = prepppdata[1].m_Prn[i] - MIN_BDS_SATNO + 1; break; }
								case 'E': { S0prn = prepppdata[1].m_Prn[i] - MIN_GAL_SATNO + 1; break; }
								default:break;
								}
								rinexo << GetSystem(prepppdata[1].m_Prn[i]) << setw(2) << setfill('0') << S0prn << " ";
								if (prepppdata[1].m_P1[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P1[i] - DCBdis << "   "; }
								if (prepppdata[1].m_P2[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P2[i] << "   "; }
								if (prepppdata[1].m_P3[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P3[i] << "   "; }
								if (prepppdata[1].m_L1[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L1[i] + popt.simulcyclevalue[m][0] - popt.repairvalue[k][0] << "   "; }
								if (prepppdata[1].m_L2[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L2[i] + popt.simulcyclevalue[m][1] - popt.repairvalue[k][1] << "   "; }
								if (prepppdata[1].m_L3[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L3[i] + popt.simulcyclevalue[m][2] - popt.repairvalue[k][2]; }
								rinexo << endl;
								tempsim = 1; S0prn = 0; DCBdis = 0;
								break;
							}
						}
						if (tempsim == 0)
						{
							switch (GetSystem(prepppdata[1].m_Prn[i]))
							{
							case 'G': { S0prn = prepppdata[1].m_Prn[i] - MIN_GPS_SATNO + 1; DCBdis = (dcbdata.Get_DCB_P1C1(prepppdata[1].m_Prn[i]) * LIGHTSPEED * 1.0e-9); break; }
							case 'R': { S0prn = prepppdata[1].m_Prn[i] - MIN_GLO_SATNO + 1; break; }
							case 'C': { S0prn = prepppdata[1].m_Prn[i] - MIN_BDS_SATNO + 1; break; }
							case 'B': { S0prn = prepppdata[1].m_Prn[i] - MIN_BDS_SATNO + 1; break; }
							case 'E': { S0prn = prepppdata[1].m_Prn[i] - MIN_GAL_SATNO + 1; break; }
							default:break;
							}
							rinexo << GetSystem(prepppdata[1].m_Prn[i]) << setw(2) << setfill('0') << S0prn << " ";
							if (prepppdata[1].m_P1[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P1[i] - DCBdis << "   "; }
							if (prepppdata[1].m_P2[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P2[i] << "   "; }
							if (prepppdata[1].m_P3[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P3[i] << "   "; }
							if (prepppdata[1].m_L1[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L1[i] << "   "; }
							if (prepppdata[1].m_L2[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L2[i] << "   "; }
							if (prepppdata[1].m_L3[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L3[i]; }
							rinexo << endl;
							S0prn = 0; DCBdis = 0;
						}
					}
				}
				break;
			}
		}
	}
}

void Rinexoutbodypre(ofstream& rinexo, CodeDCB& dcbdata, ppp_option_t& popt, vector <zhydata>& prepppdata)
{
	vector<int>tempprn; int S0prn = 0; double DCBdis = 0;

	for (int i = 0; i < prepppdata[1].m_SatCount; i++)
	{
		for (int j = 0; j < prepppdata[0].m_SatCount; j++)
		{
			if (prepppdata[1].m_Prn[i] == prepppdata[0].m_Prn[j]) { tempprn.push_back(prepppdata[1].m_Prn[i]); break; }
		}
	}
	rinexo << "> " << popt.CurrentEpoch.m_Year << " " << setw(2) << setfill('0') << popt.CurrentEpoch.m_Month << " "
		<< setw(2) << setfill('0') << popt.CurrentEpoch.m_Day << " " << setw(2) << setfill('0') << popt.CurrentEpoch.m_Hour << " "
		<< setw(2) << setfill('0') << popt.CurrentEpoch.m_Min << " " << setw(10) << setfill('0') << fixed << setprecision(7) << round(popt.CurrentEpoch.m_Sec) << "  0" << " " << tempprn.size() << endl;
	for (int i = 0; i < prepppdata[1].m_SatCount; i++)
	{
		for (int j = 0; j < tempprn.size(); j++)
		{
			if (prepppdata[1].m_Prn[i] == tempprn[j])
			{
				switch (GetSystem(prepppdata[1].m_Prn[i]))
				{
				case 'G': { S0prn = prepppdata[1].m_Prn[i] - MIN_GPS_SATNO + 1; DCBdis = (dcbdata.Get_DCB_P1C1(prepppdata[1].m_Prn[i]) * LIGHTSPEED * 1.0e-9); break; }
				case 'R': { S0prn = prepppdata[1].m_Prn[i] - MIN_GLO_SATNO + 1; break; }
				case 'C': { S0prn = prepppdata[1].m_Prn[i] - MIN_BDS_SATNO + 1; break; }
				case 'B': { S0prn = prepppdata[1].m_Prn[i] - MIN_BDS_SATNO + 1; break; }
				case 'E': { S0prn = prepppdata[1].m_Prn[i] - MIN_GAL_SATNO + 1; break; }
				default:break;
				}
				rinexo << GetSystem(prepppdata[1].m_Prn[i]) << setw(2) << setfill('0') << S0prn << " ";
				if (prepppdata[1].m_P1[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P1[i] - DCBdis << "   "; }
				if (prepppdata[1].m_P2[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P2[i] << "   "; }
				if (prepppdata[1].m_P3[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P3[i] << "   "; }
				if (prepppdata[1].m_L1[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L1[i] << "   "; }
				if (prepppdata[1].m_L2[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L2[i] << "   "; }
				if (prepppdata[1].m_L3[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L3[i]; }
				rinexo << endl;
				S0prn = 0; DCBdis = 0;
				break;
			}
		}
	}
}

void Rinexoutbodypsim(ofstream& rinexo, CodeDCB& dcbdata, ppp_option_t& popt, vector <zhydata>& prepppdata)
{
	vector<int>tempprn; int S0prn = 0; double DCBdis = 0;

	for (int i = 0; i < prepppdata[1].m_SatCount; i++)
	{
		for (int j = 0; j < prepppdata[0].m_SatCount; j++)
		{
			if (prepppdata[1].m_Prn[i] == prepppdata[0].m_Prn[j]) { tempprn.push_back(prepppdata[1].m_Prn[i]); break; }
		}
	}
	rinexo << "> " << popt.CurrentEpoch.m_Year << " " << setw(2) << setfill('0') << popt.CurrentEpoch.m_Month << " "
		<< setw(2) << setfill('0') << popt.CurrentEpoch.m_Day << " " << setw(2) << setfill('0') << popt.CurrentEpoch.m_Hour << " "
		<< setw(2) << setfill('0') << popt.CurrentEpoch.m_Min << " " << setw(10) << setfill('0') << fixed << setprecision(7) << round(popt.CurrentEpoch.m_Sec) << "  0" << " " << tempprn.size() << endl;
	int tempsim = 0;
	for (int i = 0; i < prepppdata[1].m_SatCount; i++)
	{
		tempsim = 0;
		for (int j = 0; j < tempprn.size(); j++)
		{
			if (prepppdata[1].m_Prn[i] == tempprn[j])
			{
				switch (GetSystem(prepppdata[1].m_Prn[i]))
				{
				case 'G': { S0prn = prepppdata[1].m_Prn[i] - MIN_GPS_SATNO + 1; DCBdis = (dcbdata.Get_DCB_P1C1(prepppdata[1].m_Prn[i]) * LIGHTSPEED * 1.0e-9); break; }
				case 'R': { S0prn = prepppdata[1].m_Prn[i] - MIN_GLO_SATNO + 1; break; }
				case 'C': { S0prn = prepppdata[1].m_Prn[i] - MIN_BDS_SATNO + 1; break; }
				case 'B': { S0prn = prepppdata[1].m_Prn[i] - MIN_BDS_SATNO + 1; break; }
				case 'E': { S0prn = prepppdata[1].m_Prn[i] - MIN_GAL_SATNO + 1; break; }
				default:break;
				}
				for (int m = 0; m < popt.simnum; m++)
				{
					if (popt.simulcyclesit[m] == tempprn[j])
					{
						rinexo << GetSystem(prepppdata[1].m_Prn[i]) << setw(2) << setfill('0') << S0prn << " ";
						if (prepppdata[1].m_P1[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P1[i] - DCBdis << "   "; }
						if (prepppdata[1].m_P2[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P2[i] << "   "; }
						if (prepppdata[1].m_P3[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P3[i] << "   "; }
						if (prepppdata[1].m_L1[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L1[i] + popt.simulcyclevalue[m][0] << "   "; }
						if (prepppdata[1].m_L2[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L2[i] + popt.simulcyclevalue[m][1] << "   "; }
						if (prepppdata[1].m_L3[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L3[i] + popt.simulcyclevalue[m][3]; }
						rinexo << endl;
						tempsim = 1; S0prn = 0; DCBdis = 0;
						break;
					}
				}
				if (tempsim == 0)
				{
					rinexo << GetSystem(prepppdata[1].m_Prn[i]) << setw(2) << setfill('0') << S0prn << " ";
					if (prepppdata[1].m_P1[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P1[i] - DCBdis << "   "; }
					if (prepppdata[1].m_P2[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P2[i] << "   "; }
					if (prepppdata[1].m_P3[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_P3[i] << "   "; }
					if (prepppdata[1].m_L1[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L1[i] << "   "; }
					if (prepppdata[1].m_L2[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L2[i] << "   "; }
					if (prepppdata[1].m_L3[i] != 0.0) { rinexo << setw(13) << setfill(' ') << fixed << setprecision(3) << prepppdata[1].m_L3[i]; }
					rinexo << endl;
					S0prn = 0; DCBdis = 0;
					break;
				}
			}
		}
	}

	ofstream ofsim;
	ofsim.open("C:\\Users\\15149\\Desktop\\shiy\\sim", ios::app);
	ofsim << "> " << popt.CurrentEpoch.m_Year << " " << setw(2) << setfill('0') << popt.CurrentEpoch.m_Month << " "
		<< setw(2) << setfill('0') << popt.CurrentEpoch.m_Day << " " << setw(2) << setfill('0') << popt.CurrentEpoch.m_Hour << " "
		<< setw(2) << setfill('0') << popt.CurrentEpoch.m_Min << " " << setw(10) << setfill('0') << fixed << setprecision(7) << popt.CurrentEpoch.m_Sec << "  0" << " " << tempprn.size() << ' ' << popt.simnum << ' ' << popt.simnums << endl;
	for (int m = 0; m < popt.simnum; m++)
	{
		ofsim << popt.simulcyclesit[m] << ' ';
		for (int n = 0; n < 3; n++)
		{
			ofsim << "!" << int(popt.simulcyclevalue[m][n]) << ' ';
		}
		if (m == popt.simnum - 1) { ofsim << endl; }
	}
}

/// Spp MatrixT for calculate
int SPP_CalMatrixT(ppp_option_t& popt, PPPObsData& sppdata, double* Sta, gnsstime& t,
	PreciseData* predata, ClockData* clkdata, pcv_t* AntSat,
	MatrixT& B, MatrixT& L, SppOutputPara& sppout, SatInfo& satinfo, CodeDCB& m_dcb,
	GIM_t iondata,
	vector<double>& ion, vector<double>& stec, vector<double>& trop,
	vector<double>& satclks,      
	vector<double>& satclks_prev,
	vector<unsigned int>& satclks_prn,     
	vector<unsigned int>& satclks_prn_prev  
)
{
	unsigned int FreqNum = 0, MathModel = 0, Pmode = 0;
	unsigned int i = 0, k = 0, SatCount = 0;
	unsigned int sysm[255] = { 0 }, sysnum = 1;
	int sittemp = 0;

	// 初始化



	std::cout << "SPP_CalMatrixT 开始: ion.size=" << ion.size() << ", trop.size=" << trop.size() << std::endl;

	FreqNum = popt.freqn;                    // Frequence number
	MathModel = popt.math_model;             // Math model
	Pmode = popt.process_mode;               // Process mode

	SatCount = B.GetRow();
	for (int k = 0; k < 4; k++)//
	{
		if (popt.ISBsystem[k] == 'C' || popt.ISBsystem[k] == 'B' || popt.ISBsystem[k] == 'E' || popt.ISBsystem[k] == 'R')
		{
			sittemp++;
			popt.isbsit[k] = sittemp;
		}
	}
	/*for (int k = 0; k < 4; k++)
	{
		if (popt.ISBsystem[k] != ' ')
		{
			sittemp++;
			switch (popt.ISBsystem[k])
			{
			case 'R': popt.isbsit[0] = sittemp; break;
			case 'C': popt.isbsit[1] = sittemp; break;
			case 'B': popt.isbsit[2] = sittemp; break;
			case 'E': popt.isbsit[3] = sittemp; break;
			default: break;
			}
		}
	}*/

	for (i = 0; i < SatCount; ++i)
	{
		if (sppdata.m_Prn[i] == 0)	continue;

		gnsstime obst_b = t;
		gnsstime obst = t;
		unsigned int satno = sppdata.m_Prn[i];
		int SlotNo = 0;

		char sys = 0;
		unsigned int prn = 0;
		sys = GetSysPrnGREC3C2(satno, prn);

		// Get navigation ephemeris when realtime ppp,realtime spp and post spp
		void* eph = NULL;
		int b = (popt.System & PRO_SYS_GLO);
		if ((Pmode == MODE_SPP && !(popt.use_orbclk)))//
		{
			eph = Find_NavEph(satno, obst_b, predata->m_NavData);          // Only get navigation

			if (eph == NULL)
			{
				// Log in file
				FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
				if (logfp)
					fprintf(logfp, "  Find navigation of SAT (%c%02d) is NULL. \n", sys, prn);
				if (logfp) fclose(logfp);
				// log end

				continue;
			}
		}
		else if (popt.process_time == MODE_POST && (popt.System & PRO_SYS_GLO))  // 0: post process 1: realtime process 格洛纳斯需要广播星历
		{
			if (GetSystem(satno) == 'R') {
				eph = Find_NavEph(satno, obst_b, predata->m_NavData);
				if (eph == NULL)
				{
					// Log in file
					FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
					if (logfp)
						fprintf(logfp, "  Find GLONASS navigation of SAT (%c%02d) is NULL. \n", sys, prn);
					if (logfp) fclose(logfp);
					// log end

					continue;
				}
			}

		}
		// find navigation end

		double satpos_B[3] = { 0 }, satpos_P[3] = { 0 }, satpos[3] = { 0 }, satvel[3] = { 0 }, satpco[3] = { 0 }, satclk_B = 0.0, satclk_P = 0.0, satclk = 0.0, satdcb = 0.0;
		double Ele = 0.0, Azimuth = 0.0, satvel_B[3] = { 0 }, satvel_P[3] = { 0 }; double tgd1 = 0.0, tgd2 = 0.0;
		double Ele_P = 0.0, Azimuth_P = 0.0;
		if (popt.process_time == MODE_REALTIME)
		{
			if (Cal_Sat_PosVel_Real(sppdata.m_Prn[i], obst_b, sppdata.m_Pr1[i], Sta, eph, satpos_B, satvel_B, satclk_B, popt.use_orbclk, tgd1, tgd2) == 0)
			{
				// Log in file
				FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
				if (logfp)
					fprintf(logfp, "  Cal pos/vel and clock of SAT (%c%02d) failed in realtime SPP. \n", sys, prn);
				if (logfp) fclose(logfp);
				// log end

				continue;
			}
		}
		if (popt.process_time == MODE_POST)
		{
			if (popt.use_orbclk)
			{
				if (predata == NULL || clkdata == NULL)
				{
					// Log in file
					FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
					if (logfp)
						fprintf(logfp, "  Precise Ephemeris and Clock Data are NULL in SPP. \n");
					if (logfp) fclose(logfp);
					// log end

					printf("Error:Precise Ephemeris and Clock Data are null.\n");
					return 0;
				}

				double SatPCO[3] = { 0.0 };

				int clk_len = 0, pre_len = 0;
				int ml = 10, ml_clk = 5;

				Clk_EpochData* ptr_clk = clkdata->Get_Clock_Interp_Data(sppdata.m_Prn[i], obst, clk_len, ml_clk);
				PreEpochData* ptr_pre = predata->Get_Precise_InterpData(sppdata.m_Prn[i], obst, pre_len, ml);


				gnsstime sun_t = obst;
				double* gmst = NULL;
				GetSunMoonPos(sun_t, popt.SunPos, popt.MoonPos, gmst);

				if (!Cal_Sat_PosVel_Post(satno, obst, sppdata.m_Pr1[i], Sta, ptr_pre, pre_len, ptr_clk, clk_len, satpos_P, satvel_P, satclk_P))
					continue;

				if (ptr_clk) free(ptr_clk);
				if (ptr_pre) free(ptr_pre);
			}
			else
			{
				if (eph == NULL) return 0;
				if (Cal_Sat_PosVel_Real(satno, obst, sppdata.m_Pr1[i], Sta, eph, satpos, satvel, satclk, 0, tgd1, tgd2) == 0)
				{
					// Log in file
					FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
					if (logfp)
						fprintf(logfp, "  Cal pos/vel and clock of SAT (%c%02d) failed by navigation in SPP. \n", sys, prn);
					if (logfp) fclose(logfp);
					// log end

					continue;
				}
			}
		}

#if 0
		for (int qq = 0; qq < 3; qq++)
		{
			satpos[qq] = satpos_B[qq];
			satvel[qq] = satvel_B[qq];
		}

		satclk = satclk_B;
#endif // 0
#if 1
		for (int qq = 0; qq < 3; qq++)
		{
			satpos[qq] = satpos_P[qq];
			satvel[qq] = satvel_P[qq];
		}

		satclk = satclk_P;
#endif // 0







		// get satellite elevation and azimuth
		Cal_Sat_EleAzimuth(Sta, satpos, Ele, Azimuth);
		
		double azel[2];
		azel[0] = Azimuth; azel[1] = Ele;
		const double azimuth = Azimuth;

		double P1 = 0.0;       // Raw observations
		double P2 = 0.0;
		double IonDelay = 0, var_ion = 0, S_TEC = 0;                       // Ionosphere delay for single frequency

		if (FreqNum == 1)
		{
			double w1, w2;
			P1 = sppdata.m_Pr1[i];

			if (1) 
			{
				if (GetSystem(satno) == 'C')
				{
					double ptr_dcb = m_dcb.Get_DCB_B1B3(sppdata.m_Prn[i]);
					satdcb = ptr_dcb * 1E-9 * LIGHTSPEED;
					double gamma = BDS_B1_B3_gamma;
					P1 = P1 + satdcb * BDS_B1_B3_gamma;
				}
				if (GetSystem(satno) == 'G')
				{
					double ptr_dcb = m_dcb.Get_DCB_P1P2(sppdata.m_Prn[i]);
					satdcb = ptr_dcb * 1E-9 * LIGHTSPEED;

					P1 = P1 + satdcb * BDS_L1_L2_gamma;
				}
			}
			
		}
		double Lat = 0.0, Lon = 0.0, H = 0.0;
		double TropDelay = 0.0, Trop_dry = 0.0, Trop_wet = 0.0;
		double Trop_dry_s = 0, Trop_wet_s = 0;
		double Trop_dry_func = 0, Trop_wet_func = 0;

		Coordinate blh;
		blh.XYZ.setX(Sta[0]); blh.XYZ.setY(Sta[1]); blh.XYZ.setZ(Sta[2]);
		blh._xyz2blh();
		Lat = blh.BLH._B; Lon = blh.BLH._L; H = blh.BLH._H;
		double blh3[3] = { Lat, Lon, H };
	
		double freq;

		freq = Get_Frequency(satno, 1, 0);


		//double iondalay_p1p2 = IonDelay;
		IonModel_Grid(obst.m_gtime, freq, &iondata, blh3, azel, 3, &IonDelay, &var_ion, &S_TEC);
		ion.push_back(IonDelay);

	


		stec.push_back(S_TEC);
		//double chazhi = IonDelay - iondalay_p1p2;
		TropDelay = Trop_NMF_Saastamoinen(obst.m_Doy, Lat, H, Ele, Trop_dry, Trop_wet, 0);
		
		trop.push_back(TropDelay);  
	
		satclks.push_back(satclk);  
		satclks_prn.push_back(satno);


		double RelativeDelay = 0.0;
		RelativeDelay = RelativeEffect(satpos, satvel);

		double SagnacDelay = 0.0;
		SagnacDelay = SagnacEffect(satpos, Sta);

		double rs = 0.0;
		rs = sqrt((Sta[0] - satpos[0]) * (Sta[0] - satpos[0]) + (Sta[1] - satpos[1]) * (Sta[1] - satpos[1]) + (Sta[2] - satpos[2]) * (Sta[2] - satpos[2]));

		// coefficient of position
		int est_xyz = 1;
		if (est_xyz)
		{
			B(k + 1, 1) = (Sta[0] - satpos[0]) / rs;
			B(k + 1, 2) = (Sta[1] - satpos[1]) / rs;
			B(k + 1, 3) = (Sta[2] - satpos[2]) / rs;
		}

		// coefficient of receiver clock
		B(k + 1, est_xyz * 3 + popt.estpara.rclk) = 1.0;
		printf("ISB parameters: R=%d, C=%d, B=%d, E=%d\n",
			popt.isbsit[0], popt.isbsit[1], popt.isbsit[2], popt.isbsit[3]);


		switch (sys) {
		case 'G':
			// GPS - 无ISB
			break;
		case 'R':
			if (popt.isbsit[0] > 0) {
				int col = est_xyz * 3 + popt.estpara.rclk + popt.isbsit[0];
				if (col <= B.GetCol()) {
					B(k + 1, col) = 1.0;
				}
				else {
					//	printf("ERROR: Column index %d out of bounds for GLONASS PRN %d\n", col, satno);
				}
			}
			break;
		case 'C':
			if (popt.isbsit[1] > 0) {
				int col = est_xyz * 3 + popt.estpara.rclk + popt.isbsit[1];
				if (col <= B.GetCol()) {
					B(k + 1, col) = 1.0;
				}
				else {
					//	printf("ERROR: Column index %d out of bounds for BDS PRN %d\n", col, satno);
				}
			}
			break;
		case 'E':
			if (popt.isbsit[3] > 0) {
				int col = est_xyz * 3 + popt.estpara.rclk + popt.isbsit[3];
				if (col <= B.GetCol()) {
					B(k + 1, col) = 1.0;
				}
				else {
					//printf("ERROR: Column index %d out of bounds for Galileo PRN %d\n", col, satno);
				}
			}
			break;
		case 'B':
			if (popt.isbsit[2] > 0) {
				int col = est_xyz * 3 + popt.estpara.rclk + popt.isbsit[2];
				if (col <= B.GetCol()) {
					B(k + 1, col) = 1.0;
				}
				else {
					//printf("ERROR: Column index %d out of bounds for BDS2 PRN %d\n", col, satno);
				}
			}
			break;
		default:
			//printf("WARNING: Unknown system '%c' for PRN %d\n", sys, satno);
			break;
		}
		// Observation - Computation (O-C)
		int a = sppdata.m_Prn[i];
		rs = rs + SagnacDelay;
		double DTS = satclk * LIGHTSPEED;
		L(k + 1, 1) = P1 - (rs - DTS + TropDelay + IonDelay - RelativeDelay);

		sppout.Elev[k] = Ele;                                 // Satellite Elevation
		sppout.Azim[k] = Azimuth;                             // Satellite Azimuth
		sppout.Trop[k] = Trop_dry + Trop_wet;                 // Troposphere zenith Delay
		sppout.Sclk[k] = satclk;                              // Satellite Clock offset
		sppout.Spos[k * 3 + 0] = satpos[0];                   // Satellite Position
		sppout.Spos[k * 3 + 1] = satpos[1];
		sppout.Spos[k * 3 + 2] = satpos[2];

		//% Satellite information storage
		satinfo.flag[k] = 1;
		satinfo.m_Prn[k] = satno;
		satinfo.m_Ele[k] = Ele;
		satinfo.m_Azi[k] = Azimuth;
		satinfo.m_Clk[k] = satclk;
		satinfo.m_Xs[k * 3 + 0] = satpos[0];
		satinfo.m_Xs[k * 3 + 1] = satpos[1];
		satinfo.m_Xs[k * 3 + 2] = satpos[2];
		satinfo.m_Vs[k * 3 + 0] = satvel[0];
		satinfo.m_Vs[k * 3 + 1] = satvel[1];
		satinfo.m_Vs[k * 3 + 2] = satvel[2];
		//%--------------------------------

		++k;
	}
	/*std::cout << "SPP_CalMatrixT 结束: ion.size=" << ion.size() << ", 内容=[";
	for (size_t i = 0; i < ion.size(); ++i) {
		std::cout << ion[i] << (i < ion.size() - 1 ? ", " : "");
	}
	std::cout << "]" << std::endl;*/
	std::cout << "satclks: [";
	for (size_t i = 0; i < satclks.size(); ++i) {
		std::cout << satclks[i] << " ";
	}
	std::cout << "]" << std::endl;
	
	sppout.m_Size = k;
	satinfo.m_Size = k;
	return k;
}

double CalDOP(MatrixT Q)
{
	double gdop = 0;
	double var_lat = 0, var_lon = 0, var_h = 0;
	double var_x = sqrt(Q(1, 1)), var_y = sqrt(Q(2, 2)), var_z = sqrt(Q(3, 3));
	double var_rclk = sqrt(Q(4, 4));

	Cart_Crd var_xyz(var_x, var_y, var_z);
	Coordinate var_cord(var_xyz);
	var_cord._xyz2blh();

	var_lat = var_cord.BLH._B; var_lon = var_cord.BLH._L; var_h = var_cord.BLH._H;

	//gdop = sqrt(SQR(var_lat) + SQR(var_lon) + SQR(var_h) + SQR(var_rclk));
	gdop = sqrt(SQR(var_x) + SQR(var_y) + SQR(var_z) + SQR(var_rclk));

	return gdop;
}

/// Estimate Parameter by LST
int SPPEst(ppp_option_t& popt, gnsstime& t, GenData& gendata,
	PreciseData* predata, ClockData* clkdata,
	unsigned int& satnum, PPPObsData& sppdata, SatInfo& satinfo,
	vector<double>& ion, vector<double>& STEC, vector<double>& trop,
	vector<double>& satclks, vector<double>& satclks_prev,
	vector<unsigned int>& satclks_prn,      
	vector<unsigned int>& satclks_prn_prev  
)

{
	int est_xyz = 1;

	int gd = 0, itera = 0, bkf = 1;
	int SatCount = 0, satntmp = 0;

	double Sta[3] = { 0.0 };
	double xp[3] = { 0.0 }, dx[3] = { 0.0 }, maxv = 0.0;
	double drclk = 0.0, rclkp = 0.0;
	pcv_t* ansat = gendata.m_pcvs;
	GIM_t Iondata = gendata.GIMdata;

	if (Norm(popt.StaPos_CRD, 3) > 0) {
		Sta[0] = popt.StaPos_CRD[0];
		Sta[1] = popt.StaPos_CRD[1];
		Sta[2] = popt.StaPos_CRD[2];
	}
	else if (Norm(popt.ProXYZ, 3) > 0) {
		Sta[0] = popt.ProXYZ[0];
		Sta[1] = popt.ProXYZ[1];
		Sta[2] = popt.ProXYZ[2];
	}

	// Check SSR ---------------------------------------
	if ((popt.process_mode == MODE_SPP) && (!popt.use_orbclk)) {
		for (int i = 0; i < sppdata.m_SatCount; ++i) {
			unsigned int satno = sppdata.m_Prn[i];
			char sys = 0;
			unsigned int prn = 0;
			sys = GetSysPrn(satno, prn);

			void* eph = NULL;
			eph = Find_NavEph(satno, t, predata->m_NavData);

			if (eph == NULL) {
				FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
				if (logfp)
					fprintf(logfp, "  Find navigation or SSR of SAT (%c%02d) is NULL. \n", sys, prn);
				if (logfp) fclose(logfp);
				DeleteSatData(sppdata, i);
				satnum--;
				i--;
				continue;
			}
		}
	}
	// Check SSR end --------------------------------------
	MatrixT V_1(satnum, 1);
	int zhy_ = 0;
	int kk = 0;
	do {
		kk++;
		int i = 0, j = 0, k = 0, gps = 0, glo = 0, gal = 0, bds3 = 0, bds2 = 0;
		int B_2 = 0, B_3 = 0, is_issb = 0;
		int flag = 0, ParaNum = est_xyz * 3 + popt.estpara.rclk;

		gd = 0;
		popt.ISBsystem[4] = ' ';
		for (i = 0; i < satnum; ++i) {
			unsigned int satno = sppdata.m_Prn[i];
			char sys = GetSystem_GREC2C3(satno);
			printf("Satellite %d: PRN=%u, System=%c\n", i, satno, sys);
			switch (GetSystem_GREC2C3(sppdata.m_Prn[i])) {
			case 'G': { gps = 1; break; }
			case 'R': { glo = 1; popt.ISBsystem[0] = 'R'; break; }
			case 'C': { bds3 = 1; popt.ISBsystem[1] = 'C'; break; }
			case 'B': { bds2 = 1; popt.ISBsystem[2] = 'B'; break; }
			case 'E': { gal = 1; popt.ISBsystem[3] = 'E'; break; }
			default: break;
			}
		}
		is_issb = gps + glo + bds3 + bds2 + gal - 1;
		ParaNum = ParaNum + is_issb;
		std::cout << "ParaNum " << ParaNum << std::endl;

		if ((popt.System & PRO_SYS_GPS) && (gps == 0)) {
			gd = 0;
			break;
		}

		if (satnum < ParaNum) return 0;

		MatrixT B;                      // Coefficient MatrixT
		MatrixT L;                      // Constant MatrixT
		MatrixT P;                      // Weight MatrixT
		B.Resize(satnum, ParaNum); L.Resize(satnum, 1); P.Resize(satnum, satnum);

		SppOutputPara sppout(satnum);

		
		SatCount = SPP_CalMatrixT(popt, sppdata, Sta, t, predata, clkdata, ansat,
			B, L, sppout, satinfo, gendata.m_Dcb, Iondata,
			ion, STEC, trop,
			satclks, satclks_prev,
			satclks_prn, satclks_prn_prev);  
		for (size_t i = 0; i < satclks.size(); ++i) {
			std::cout << satclks[i] << (i < satclks.size() - 1 ? ", " : "");
		}
		std::cout << "]" << std::endl;

		if (SatCount < ParaNum) return 0;

		satntmp = satnum;
		for (i = 0; i < satnum; ++i) {
			if (satnum < ParaNum) break;
			for (j = 0; j < satinfo.m_Size; ++j) {
				if (sppdata.m_Prn[i] == satinfo.m_Prn[j]) break;
			}
			if (j == satinfo.m_Size) DeleteSatData(sppdata, i--);
			satnum = satinfo.m_Size;
		}

		if (satnum < ParaNum) return 0;

		MatrixT Q(ParaNum, ParaNum);
		MatrixT X(ParaNum, 1);
		MatrixT V(satnum, 1);
		double xigema0 = 0, vtpv, sig = 0.0;
		MatrixT BT(ParaNum, satnum);
		MatrixT BN(satnum, ParaNum);
		MatrixT BNB(satnum, satnum);
		MatrixT VTPV;
		MatrixT Qvv(satnum, satnum);
		for (i = 0; i < satnum; ++i) {
			char sys = 0;
			double sinel = 0.0, sys_weight = 0.0, power = 0.0;

			sys = GetSystem(sppdata.m_Prn[i]);
			sinel = sin(sppout.Elev[i]);
			sys_weight = (sys == 'G') ? 3 : 1.0;
			power = sinel * sinel * sys_weight;

			P(i + 1, i + 1) = power;
		}

		if (LeastSquare(B, L, P, V, Q, X) < 0) {
			printf("MatrixT invert failed!\n");
			itera = 20;
			break;
		}

		if (isnan(V(1, 1))) {
			printf("SPP calculate NAN!\n");
			gd = 0;
			break;
		}

		for (int m = 0; m < satnum; ++m) {
			sppout.Resp[m] = V(m + 1, 1);
		}

		if (est_xyz) {
			Sta[0] += X(1, 1); Sta[1] += X(2, 1); Sta[2] += X(3, 1);
		}

		popt.Rclk = X(est_xyz * 3 + popt.estpara.rclk, 1);
		for (int i = 0; i < is_issb; i++) {
			popt.ISB[i] = X(5 + i, 1);
		}
		drclk = (popt.Rclk - rclkp) / LIGHTSPEED;
		rclkp = popt.Rclk;

		if (est_xyz) {
			dx[0] = fabs(Sta[0] - xp[0]);
			dx[1] = fabs(Sta[1] - xp[1]);
			dx[2] = fabs(Sta[2] - xp[2]);

			xp[0] = Sta[0]; xp[1] = Sta[1]; xp[2] = Sta[2];
		}

		satntmp = satnum;
		sppdata.m_SatCount = satntmp;
		double resm = 0.0;
		if (flag || (dx[0] < 200 && dx[1] < 200 && dx[2] < 200)) {
			maxv = V.GetMax(k);
			sppout.Elev[k] = maxv > 15 ? 0.0 : sppout.Elev[k];

			for (i = 0; i < sppout.m_Size; ++i) {
				if (abs(sppout.Elev[i]) < ELE_THRES(popt.ele_mask)) {
					if (popt.trace) {
						FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
						if (sppout.Elev[i] == 0) {
							if (logfp)
								fprintf(logfp, "  Delete Sat for too large residuals in SPP: %c%3d\n",
									GetSystem(sppdata.m_Prn[i]), sppdata.m_Prn[i]);
						}
						else {
							if (logfp)
								fprintf(logfp, "  Delete Sat for lower %d in SPP: %c%3d(Ele:%4.1f) \n",
									(int)popt.ele_mask, GetSystem(sppdata.m_Prn[i]), sppdata.m_Prn[i], (sppout.Elev[i] * R2D));
						}
						if (logfp) fclose(logfp);
					}
					DeleteSatData(sppdata, i);
					sppout.DeleSat(i);
					i--;

					if (RankDetect(popt, sppout.m_Size)) {
						FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
						if (logfp)
							fprintf(logfp, "  Equation Number < Parameter Number in SPP for ElevMask. \n");
						if (logfp) fclose(logfp);
						return 0;
					}
					continue;
				}
				resm += (sppout.Resp[i] * sppout.Resp[i]);
			}
			flag = 1;
		}
		else {
			flag = 0;
		}
		satnum = satntmp;

		double GDOP = 0.0;
		if (est_xyz)
			GDOP = CalDOP(Q);

		gd = (GDOP < 30) ? 1 : 0;
		if (GDOP < 30) {
			gd = 1;
		}
		else {
			FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
			if (logfp)
				fprintf(logfp, "  GDOP is over 30 in SPP (GDOP=%5.1f). \n", GDOP);
			if (logfp) fclose(logfp);
			gd = 0;
		}

		if (gd)
			gd = sqrt(resm / satnum) < 100 ? 1 : 0;

		popt.GDOP = GDOP;
		popt.StaPos_SPP[0] = Sta[0];
		popt.StaPos_SPP[1] = Sta[1];
		popt.StaPos_SPP[2] = Sta[2];

		if (est_xyz) {
			if (dx[0] < 1.0e-4 && dx[1] < 1.0e-4 && dx[2] < 1.0e-4) {
				bkf = 0;
			}
		}
		else {
			if (drclk < 1e-4)
				bkf = 0;
		}

		if ((++itera) > 100) {
			FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
			if (logfp)
				fprintf(logfp, "  Iteration time is over 20 times in SPP (iteration=%3d). \n", itera);
			if (logfp) fclose(logfp);
			gd = 0;
			break;
		}
		if (gd) {
			EpochSolution sol;
			sol.xyz[0] = Sta[0];
			sol.xyz[1] = Sta[1];
			sol.xyz[2] = Sta[2];
			sol.rclk = popt.Rclk;
			sol.t = t;
			epoch_solutions[t] = sol;
			std::cout << "存储历元解: 时间=" << t
				<< ", X=" << Sta[0] << ", Y=" << Sta[1] << ", Z=" << Sta[2]
				<< ", Rclk=" << popt.Rclk << std::endl;
			std::cout << "存储时间: " << t.m_Year << "-" << t.m_Month << "-" << t.m_Day
				<< " " << t.m_Hour << ":" << t.m_Min << ":" << t.m_Sec << std::endl;
			while (epoch_solutions.size() > 2) {
				auto oldest_it = epoch_solutions.begin(); // 获取时间最早的历元
				std::cout << "删除最旧历元: 时间=" << oldest_it->first << std::endl;
				epoch_solutions.erase(oldest_it);
			}
		}
		std::cout << "epoch_solutions 大小: " << epoch_solutions.size() << std::endl;
		for (const auto& sol_entry : epoch_solutions) {
			std::cout << "当前存储的历元: 时间=" << sol_entry.first
				<< ", X=" << sol_entry.second.xyz[0] << ", Y=" << sol_entry.second.xyz[1]
				<< ", Z=" << sol_entry.second.xyz[2] << ", Rclk=" << sol_entry.second.rclk << std::endl;
		}
	} while (bkf);

	return gd;
}
void IGGIII(MatrixT& v, MatrixT& p, int satnum)
{

	double k0 = 1.5, k1 = 5;


	for (int i = 1; i < satnum + 1; i++)
	{
		if (v(i, 1) > k0 && v(i, 1) <= k1)
			p(i, i) = (p(i, i) * (k0 * pow((k1 - v(i, 1)) / (k1 - k0), 2))) / v(i, 1);
		if (v(i, 1) > k1)
			p(i, i) = p(i, i) / 100;
	}

};
/// Code Standard Point Positioning
int CodeSPP(ppp_option_t& popt, ObsEpochData& obsdata, GenData& gendata,
	PreciseData* predata, ClockData* clkdata, SatInfo& satinfo,
	vector<double>& ion, vector<double>& stec, vector<double>& trop,
	vector<double>& satclks, vector<double>& satclks_prev,
	vector<unsigned int>& satclks_prn,      
	vector<unsigned int>& satclks_prn_prev  
)

{
	unsigned int SatCount = 0;
	int status = 0;

	SatCount = obsdata.m_SatCount;

	PPPObsData sppdata(SatCount);	printf("经过 GetBlockData 有效卫星数：%d\n", SatCount);
	SatCount = sppdata.GetBlockData(popt, &gendata, obsdata, predata, clkdata, -1);
	printf("原始观测卫星数：%d\n", obsdata.m_SatCount);
	printf("经过 GetBlockData 后有效卫星数：%d\n", SatCount);


	if (RankDetect(popt, SatCount))
	{
		// Log in file
		FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
		if (logfp)
			fprintf(logfp, "  Satellite Number < Parameter Number in SPP. \n");
		if (logfp) fclose(logfp);
		// log end

		return 0;
	}

	// SPP Estimation
	satinfo.Initialize(SatCount);



	status = SPPEst(popt, obsdata.m_EpochTime, gendata, predata, clkdata,
		SatCount, sppdata, satinfo,
		ion, stec, trop,
		satclks, satclks_prev,
		satclks_prn, satclks_prn_prev);  


	popt.SatN_SPP = SatCount;
	for (int r = 0; r < satinfo.m_Size; r++)
	{
		satinfo.m_Prn[r];
	}
	// SPP Estimation end

	if (status == 0)
	{
		// Log in file
		FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
		if (logfp)
			fprintf(logfp, "  SPP failed...\n");
		if (logfp) fclose(logfp);
		// log end

		return 0;
	}

	return 1;
}

LeastSquarePPP::LeastSquarePPP()
{
	Initialize();
}

LeastSquarePPP::LeastSquarePPP(ppp_option_t& popt)
{
	Initialize();

	_traceMode = popt.trace;
}

LeastSquarePPP::~LeastSquarePPP()
{
	Clear();
}

void LeastSquarePPP::ION_Initialize(ION_x& ion)
{
	ion.dgf_value.clear();
	ion.gf_value.clear();
	ion.ele.clear();
	ion.num = 0;
	ion.time.clear();
	ion.a = 0;
	ion.b = 0;
	ion.c = 0;
}
void LeastSquarePPP::ION_del_one_first(ION_x& ion)
{
	if (ion.num < 10)return;
	vector<double>a = ion.dgf_value; vector<int>b = ion.time;
	vector<double>a1 = ion.gf_value;
	a.erase(a.begin()); b.erase(b.begin());
	a1.erase(a1.begin());
	ion.dgf_value = a;
	ion.gf_value = a1;
	ion.time = b;
	ion.num -= 1;
	/*if (ion.num1 < 50)return;
	vector<double>a1 = ion.gf_value; vector<int>b1 = ion.time1;
	a1.erase(a1.begin()); b1.erase(b1.begin());
	ion.gf_value = a1;
	ion.num1 -= 1;
	ion.time1 = b1;*/
}
int LeastSquarePPP::Initialize()
{
	memset(_SystemIndex, 0, sizeof(char) * 255);

	flag_epoch_count = 0;
	_traceMode = 0;
	_MathModel = 0; _Mode = 0; _FreqNum = 0; _System = 0;
	est_xyz = 0; est_rclk = 0; est_trop = 0; est_ion = 0;
	_GIM_P1P2 = 0; _IS_DCB = 0; _IS_ISB = 0; _ISBMODE = 0;
	_GpsSatNum = _GloSatNum = _GalSatNum = _BdsSatNum = 0;


	_RefSatValidFlag = 0;    // no Reference Satellite
	_RefSatPrn = 0; _RefSatIndex = -1;

	CoarseXYZ[0] = CoarseXYZ[1] = CoarseXYZ[2] = 0.0;
	m_Rclk = m_dt = 0;

	m_SatValidN = 0;
	m_Sta[0] = m_Sta[1] = m_Sta[2] = 0.0;
	m_TropZenith = 0;
	for (int i = 0; i < MAX_GNSS_PRN_NUM; i++)
	{
		Ion_X->a = 0; Ion_X->b = 0; Ion_X->c = 0; Ion_X->num = 0;
		gf_count[i] = mw_count[i] = iCycle[i] = mw_count13[i] = 0;
		gf_cur[i] = gf_pre[i] = mw_cur[i] = mw_avg[i] = mw_var[i] = mw_pre[i] = 0.0;
		gf_cur13[i] = gf_pre13[i] = mw_cur13[i] = mw_avg13[i] = mw_pre13[i] = 0.0;
		//mw_avg13[i] = 0; mw_cur13[i] = 0; mw_pre13[i] = 0;
	}

	memset(this, 0, sizeof(LeastSquarePPP));

	return 1;
}



double LeastSquarePPP::get_ion_xishu(int type, unsigned int prn, NavData* m_NavData)
{
	unsigned int satno = prn;
	int orbn = 0; double fre1 = 0, fre_need = 0;
	if (GetSystem(satno) == 'R')
	{
		orbn = Get_Nav_Orbitn(satno, m_NavData);
		fre1 = Get_Frequency(satno, 1, orbn);
		fre_need = Get_Frequency(satno, type, orbn);
	}
	else if (GetSystem(satno) == 'G')
	{
		fre1 = Get_Frequency(satno, 1, orbn);
		if (type == 3)
		{
			fre_need = Get_Frequency(satno, 5, orbn);
		}
		else fre_need = Get_Frequency(satno, type, orbn);
	}
	else if (GetSystem(satno) == 'C')
	{
		fre1 = Get_Frequency(satno, 1, orbn);
		if (type == 3)
		{
			fre_need = Get_Frequency(satno, 2, orbn);
		}
		else if (type == 2)
		{
			fre_need = Get_Frequency(satno, 3, orbn);
		}
		else
		{
			fre_need = Get_Frequency(satno, type, orbn);
		}
	}
	else
	{
		fre1 = Get_Frequency(satno, 1, orbn);
		fre_need = Get_Frequency(satno, type, orbn);

	}
	return SQR(fre1 / fre_need);
}

double LeastSquarePPP::get_wth(unsigned int prn, int flag, NavData* m_NavData)
{
	unsigned int satno = prn;
	int orbn = 0;
	if (GetSystem(satno) == 'R')
	{
		orbn = Get_Nav_Orbitn(satno, m_NavData);
		return Get_WaveLength(satno, flag, orbn);
	}
	else if (GetSystem(satno) == 'C')
	{
		if (flag == 2)
		{
			return Get_WaveLength(satno, 3, orbn);
		}
		else if (flag == 3)
		{
			return Get_WaveLength(satno, 2, orbn);
		}
		else
		{
			return Get_WaveLength(satno, flag, orbn);
		}
	}
	else if (GetSystem(satno) == 'G')
	{
		if (flag == 3)
		{
			double a = Get_WaveLength(satno, 5, orbn);
			return Get_WaveLength(satno, 5, orbn);
		}
		else
		{
			return Get_WaveLength(satno, flag, orbn);
		}

	}
	else return Get_WaveLength(satno, flag, orbn);
}


//flag_epoch_count
int LeastSquarePPP::Cycle_slip(
	ppp_option_t& popt, gnsstime& t,
	vector<zhydata>& prepppdata, vector<zhyinfo>& presatinfo,
	PPPObsData& pppdata, PreciseData& predata,
	vector<double>& ion_gimsub, vector<int>& commprn,
	vector<double>& isb,
	vector<double>& trop, vector<double>& trop_prev,
	vector<double>& satclks, vector<double>& satclks_prev,
	vector<unsigned int>& satclks_prn,      // 新增：当前历元PRN顺序
	vector<unsigned int>& satclks_prn_prev  // 新增：前一历元PRN顺序
)
{
	vector<double> L_P, L_C; vector<vector<double>>L_Pre, L_Cur, L_Pre1, L_Cur1;
	vector<double> P_P, P_C; vector<vector<double>>P_Pre, P_Cur;
	vector<vector<string>>fre_comm, fre_comm1; vector<string> fre;
	vector<int>obstype; vector< vector<int>>Obstype;
	vector<unsigned int> prn_commsat, prn_commsat1;
	vector<double> ion_gimsub1;
	vector<double> rs_1;
	vector<double> rs_2;
	vector<double> repairvalue;
	//vector<double> trop_prev;  // 前一历元对流层延迟
	vector<unsigned int>prn_all; vector<int>typeint_all; vector<string>fre_all;
	int YQ = 0; int comm_num = 0;
	typedef struct zhy
	{
		double e_x;
		double e_y;
		double e_z;
	};
	zhy e1, e2;
	vector<zhy>e_1;
	vector<zhy>e_2;
	vector<double>ele;
	int Gpsscount[5] = {};// GLO - GAL - BDS3 - BDS2
	int ISBco = 0;
	double L_pre = 0, L_cur = 0; double P_pre = 0, P_cur = 0; int m = 0;
	double rs1 = 0.0, rs2 = 0.0;
	double xyz1[3] = { 0 }, xyz2[3] = { 0 };
	double satpos1[3] = { 0 }, satpos2[3] = { 0 };
	double RKpos_cur[3] = { 0 };
	double RKpos_pre[3] = { 0 };
	int num = 0;
	num = 11;

	int epoch01 = 1500;//周跳历元
	int starsitu = 5;

	auto get_satclk = [](
		const vector<double>& clks,
		const vector<unsigned int>& prns,
		unsigned int target_prn) -> double
		{
			for (size_t idx = 0; idx < prns.size(); ++idx) {
				if (prns[idx] == target_prn) {
					return clks[idx];
				}
			}
			return 0.0;  
		};

	for (int i = 0; i < 3; i++)
	{
		RKpos_cur[i] = CoarseXYZ[i];
		RKpos_pre[i] = CoarseXYZ_pre[i];
	}
	vector<double>rand_num;
	vector<double> getRlenth1, getRfrequ1;
	vector<vector <double>> GetRlenth, GetRfrequ, GetRlenth1, GetRfrequ1;
	double getRlenth = 0.0, getRfrequ = 0.0;
	int systemset[4] = {};//GLO-GAL-BDS3-BDS2
	int glosit = 0, galsit = 0, bds3sit = 0, bds2sit = 0, maxsit = 0;



	int frenumcontrol1 = popt.frennn;		
	int frenumcontrol[5] = { 1,1,1,1,1 };	
	if (frenumcontrol1 == 2) { frenumcontrol[1] = 2; }
	if (frenumcontrol1 == 3) { frenumcontrol[1] = 2; frenumcontrol[2] = 3; }
	int ITDE_OTDE = popt.schemeset;

	double kkkcal = 3.0;

	double Zxave = 0.0, Zavemin = 0.0;
	popt.Zavemin = 0.0;
	

	for (int i = 0; i < ion_gimsub.size(); i++)
	{
		Zxave += (ion_gimsub[i] / ion_gimsub.size());
	}
	for (int i = 0; i < ion_gimsub.size(); i++)
	{
		Zavemin += ((ion_gimsub[i] - Zxave) * (ion_gimsub[i] - Zxave) / ion_gimsub.size());
	}
	if (abs(Zavemin) <= 0.01) { Zavemin = kkkcal * sqrt(abs(Zavemin)); }
	else if ((abs(Zavemin) > 0.01) && (abs(Zavemin) <= 0.029)) { Zavemin = (kkkcal + 0.2) * sqrt(abs(Zavemin)); }
	else { Zavemin = (kkkcal + 0.5) * sqrt(abs(Zavemin)); }

	popt.Zavemin = abs(Zavemin - (popt.Zsigma * 2.58));


	//for (int i = 0; i < 5; i++)
	//{
	//	if (i < frenumcontrol1) { frenumcontrol[i] = i + 1; }
	//	else { frenumcontrol[i] = 1; }
	//}
	// 在循环开始前添加


	for (int i = 0; i < prepppdata[0].m_SatCount; i++)
	{
		for (int j = 0; j < prepppdata[1].m_SatCount; j++)
		{
			if (prepppdata[0].m_Prn[i] == prepppdata[1].m_Prn[j] && prepppdata[1].m_Prn[j] < 400)
			{
				unsigned int numb = prepppdata[1].m_Prn[j];
				prn_commsat1.push_back(numb);
				for (int k = 0; k < type_int[j].size(); k++)//当前历元卫星频率数量
				{
					for (int l = 0; l < type_int_pre[i].size(); l++)
					{
						if (type_int[j][k] == frenumcontrol[0] && type_int_pre[i][l] == frenumcontrol[0])//第一个频率 
						{
							L_cur = prepppdata[1].m_L1[j] * get_wth(numb, 1, predata.m_NavData); L_pre = prepppdata[0].m_L1[i] * get_wth(numb, 1, predata.m_NavData);
							getRlenth = get_wth(numb, 1, predata.m_NavData); getRfrequ = 299792458.0 / getRlenth;
						}
						else if (type_int[j][k] == frenumcontrol[1] && type_int_pre[i][l] == frenumcontrol[1])//第二个频率 
						{
							L_cur = prepppdata[1].m_L2[j] * get_wth(numb, 2, predata.m_NavData); L_pre = prepppdata[0].m_L2[i] * get_wth(numb, 2, predata.m_NavData);
							getRlenth = get_wth(numb, 2, predata.m_NavData); getRfrequ = 299792458.0 / getRlenth;
						}
						else if (type_int[j][k] == frenumcontrol[2] && type_int_pre[i][l] == frenumcontrol[2])//第三个频率
						{
							L_cur = prepppdata[1].m_L3[j] * get_wth(numb, 3, predata.m_NavData); L_pre = prepppdata[0].m_L3[i] * get_wth(numb, 3, predata.m_NavData);
							getRlenth = get_wth(numb, 3, predata.m_NavData); getRfrequ = 299792458.0 / getRlenth;
						}
						else if (type_int[j][k] == frenumcontrol[3] && type_int_pre[i][l] == frenumcontrol[3])//第四个频率 
						{
							L_cur = prepppdata[1].m_L4[j] * get_wth(numb, 4, predata.m_NavData); L_pre = prepppdata[0].m_L4[i] * get_wth(numb, 4, predata.m_NavData);
						}
						else if (type_int[j][k] == frenumcontrol[4] && type_int_pre[i][l] == frenumcontrol[4])//第五个频率 
						{
							//if (flag_epoch_count == epoch01 && prn_commsat.size() > starsitu)//flag_epoch_count  代表历元数 
							L_cur = prepppdata[1].m_L5[j] * get_wth(numb, 5, predata.m_NavData); L_pre = prepppdata[0].m_L5[i] * get_wth(numb, 5, predata.m_NavData);
						}
						if (L_pre != 0) {
							prn_all.push_back(numb); typeint_all.push_back(type_int[j][k]); fre_all.push_back(type_[j][k + 1]);
							L_P.push_back(L_pre); L_C.push_back(L_cur); L_pre = 0; fre.push_back(type_[j][k + 1]); obstype.push_back(type_int[j][k]);
							getRlenth1.push_back(getRlenth); getRfrequ1.push_back(getRfrequ);
							break;
						}
					}
				}
				fre_comm1.push_back(fre); L_Pre1.push_back(L_P); L_Cur1.push_back(L_C);
				Obstype.push_back(obstype); GetRlenth1.push_back(getRlenth1); GetRfrequ1.push_back(getRfrequ1);
				fre.clear(); L_P.clear(); L_C.clear(); obstype.clear(); getRlenth1.clear(); getRfrequ1.clear();
				break;
			}
		}
	}

	
	for (int i = 0; i < prn_commsat1.size(); i++)
	{
		
		if (i == prn_commsat1.size() - 1)
		{
			ISBco = Gpsscount[0] + Gpsscount[1] + Gpsscount[2] + Gpsscount[3] + Gpsscount[4] - 1;
			int temp001 = 1;
			for (int j = 0; j < 4; j++)
			{
				if (Gpsscount[j + 1] != 0) { systemset[j] = temp001; temp001++; }
				maxsit = systemset[j] > maxsit ? systemset[j] : maxsit;
			}
		}
	}

	ofstream ofsatnumgg;
	//ofsatnumgg.open("C:\\Users\\15149\\Desktop\\shiy\\ofsatnumgg.out", ios::app);
	//ofsatnumgg << setw(6) << setfill(' ') << fixed << gpsnnum << setw(6) << setfill(' ') << fixed << galinnum << endl;

	ofsatnumgg.open("D:\\data10.0\\out\\satnum.out", ios::app);
	ofsatnumgg << setw(6) << setfill(' ') << fixed << prn_commsat1.size() << endl;

	FILE* fp = NULL;
	ofstream ofsitu, ofx, ofx1, ofx2, ofx3, ofxion, ofxslip, ofxorigincount, ofsigmaxpre, ofrecord, ofXX, oftime;
	string base, Zcountbase, ofsitub, ofxb, ofx1b, ofx2b, ofx3b, ofxionb, ofxslipb, ofxorigincountb, ofsigmaxpreb, ofrecordb, ofXXb, oftimeb, X11b;
	base = "D:\\data10.0\\out\\shiy\\"; Zcountbase = to_string(popt.cccount);
	base += Zcountbase;
	ofsitub = "\\ofminsitu.out";
	ofxb = "\\sigmafinal.out";
	ofx2b = "\\sigmafinal1.out";
	ofx3b = "\\sigmafinalK.out";
	ofxionb = "\\ion.out";
	ofxslipb = "\\cslip.out";
	ofxorigincountb = "\\ofxorigincount.out";
	ofsigmaxpreb = "\\ofsigmapre.out";
	ofrecordb = "\\ofrecord.out";
	ofXXb = "\\ofXX.out";
	oftimeb = "\\oftime.out";
	X11b = "\\X11.out";

	if (1)//细节输出
	{
		ofxslip.open(base + ofxslipb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\cslip.out"
		ofxion.open(base + ofxionb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ion.out"
		ofsitu.open(base + ofsitub, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ofminsitu.out"
		ofrecord.open(base + ofrecordb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ofrecord.out"
		ofsigmaxpre.open(base + ofsigmaxpreb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ofsigmapre.out"
		oftime.open(base + oftimeb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\oftime.out"
		ofxorigincount.open(base + ofxorigincountb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ofxorigincount.out"
		ofXX.open(base + ofXXb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ofXX.out"
		ofx.open(base + ofxb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\sigmafinal.out"
		ofx3.open(base + ofx3b, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\sigmafinalK.out"
		ofx2.open(base + ofx2b, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\sigmafinal1.out"
		ofx1b = "\\repair.txt";
		ofx1.open(base + ofx1b, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\repair.txt"
	}
	else
	{
		string baserepir1, baserepir2, baserepir3, baserepir4, baserepir5, baserepir6, base1;
		base1 = "C:\\Users\\15149\\Desktop\\shiy\\";
		baserepir1 = "\\"; baserepir2 = popt.StaName; baserepir3 = "_repair"; baserepir4 = to_string(popt.schemeset); baserepir5 = ".txt";
		baserepir6 = to_string(popt.cccount);
		ofx1b += baserepir1; ofx1b += baserepir2; ofx1b += baserepir3; ofx1b += baserepir6; ofx1b += baserepir5;
		base1 += baserepir4;
		ofx1.open(base1 + ofx1b, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\repair.txt"
	}

	for (int i = 0; i < ISBco; i++)
	{
		//aveofd[4 + i] = isb[i];
	}

	int cstanum = 0, cstacyc[3] = {};
	//生成模拟周跳位置和数值
	int ccstasit = 0, ccstasit1 = 0, ccsimulvalue = 0, plusorminus = 0;
	int cycle_source = 100;
	for (int i = 0; i < popt.simnumc; i++)
	{
		if (flag_epoch_count == popt.simulcycletime[i] && popt.simulcycletime[i] != 0)
		{
			int tempcount = 0; popt.simnums++;
			ccstasit = rand() % prn_commsat.size(); ccstasit1 = rand() % frenumcontrol1;
			plusorminus = rand() % 2;
			if (plusorminus > 0) { ccsimulvalue = (rand() % cycle_source + 1); }
			else { ccsimulvalue = -(rand() % cycle_source + 1); }

			for (int j = 0; j < popt.simnumc; j++)
			{
				if (popt.simulcyclesit[j] == prn_commsat[ccstasit])
				{
					cstanum = prn_commsat[ccstasit]; cstacyc[ccstasit1] += ccsimulvalue;
					if (GetSystem(cstanum) == 'G')
					{
						if (ccstasit1 == 1)
						{
							popt.simulcyclevalue[j][2] += ccsimulvalue;
							break;
						}
						else
						{
							popt.simulcyclevalue[j][ccstasit1] += ccsimulvalue;
							break;
						}
					}
					else
					{
						popt.simulcyclevalue[j][ccstasit1] += ccsimulvalue;
						break;
					}
					popt.simulcyclevalue[j][ccstasit1] += ccsimulvalue;
					break;
				}
				else
				{
					tempcount++;
				}
				if (tempcount == popt.simnumc)
				{
					cstanum = prn_commsat[ccstasit]; cstacyc[ccstasit1] = ccsimulvalue;
					popt.simnum += 1;
					if (GetSystem(cstanum) == 'G')
					{
						if (ccstasit1 == 1)
						{
							popt.simulcyclesit[popt.simnum - 1] = prn_commsat[ccstasit]; popt.simulcyclevalue[popt.simnum - 1][2] = ccsimulvalue;
							break;
						}
						else
						{
							popt.simulcyclesit[popt.simnum - 1] = prn_commsat[ccstasit]; popt.simulcyclevalue[popt.simnum - 1][ccstasit1] = ccsimulvalue;
							break;
						}
					}
					else
					{
						popt.simulcyclesit[popt.simnum - 1] = prn_commsat[ccstasit]; popt.simulcyclevalue[popt.simnum - 1][ccstasit1] = ccsimulvalue;
						break;
					}
				}
			}
			break;
		}
	}

	vector<vector<double>> L_Curtt;
	L_Curtt = L_Cur;
	for (int i = 0; i < prn_commsat.size(); i++)
	{
		int temprpairsta = 0;
		for (int j = 0; j < repairstanum; j++)
		{
			if (popt.recordrepair[j] != prn_commsat[i])
			{
				temprpairsta++;
			}
			else
			{
				break;
			}
		}
		if (temprpairsta == repairstanum)
		{
			for (int j = 0; j < repairstanum; j++)
			{
				if (popt.recordrepair[j] == 0)
				{
					popt.recordrepair[j] = prn_commsat[i];
					break;
				}
			}
		}
		if (prn_commsat[i] == cstanum)
		{
			for (int j = 0; j < frenumcontrol1; j++)
			{
				L_Curtt[i][j] += (cstacyc[j] * GetRlenth[i][j]);
			}
		}
	}

	for (int i = 0; i < prn_commsat.size(); i++)
	{
		for (int j = 0; j < presatinfo[0].m_Size + 1; j++) {
			if (prn_commsat[i] == presatinfo[0].m_Prn[j] && presatinfo[0].m_Prn[j] < 400) {

				satpos1[0] = presatinfo[0].m_Xs[3 * j + 0];
				satpos1[1] = presatinfo[0].m_Xs[3 * j + 1];
				satpos1[2] = presatinfo[0].m_Xs[3 * j + 2];
				for (int k = 0; k < 3; ++k)
				{
					rs1 += SQR(RKpos_pre[k] - satpos1[k]);// (Xr-Xs)^2 +(Yr-Ys)^2+(Zr-Zs)^2
				}
				rs1 = sqrt(rs1);

				rs_1.push_back(rs1);
				rs1 = 0.0;
			}
		}
	}

	comm_num = frenumcontrol1 * prn_commsat.size();
	double sigema0 = 0.0;
	MatrixT B1, P1, L1, V1, VZ, Nbb1, X1, Pp1, Sigema0, Sigema1, C1, W1, Ncc1, Qvv, Pxx;
	MatrixT v__1(comm_num, 1);
	MatrixT p__1(comm_num, comm_num); MatrixT Sigema0__1(1, 1);

	int constnu, all_const, weizhishu;
	int mk_zhy = 0;

	all_const = prn_commsat.size();
	constnu = comm_num;

	weizhishu = 3 + _SystemN + all_const;   
	int necess_para = 4;
	B1.Resize(constnu, weizhishu);
	Nbb1.Resize(weizhishu, weizhishu);
	X1.Resize(weizhishu, 1);

	P1.Resize(constnu, constnu);
	Qvv.Resize(constnu, constnu);
	L1.Resize(constnu, 1);
	V1.Resize(constnu, 1);
	Pp1.Resize(constnu, constnu);

	Pxx.Resize(weizhishu, weizhishu);
	C1.Resize(all_const, weizhishu);
	W1.Resize(weizhishu, 1);
	Ncc1.Resize(all_const, all_const);
	bool tdcp_output_done = false;
	bool fslip_written = false;

	for (int i = 0; i < all_const; i++)
	{
		for (int j = 0; j < presatinfo[1].m_Size + 1; j++)
		{
			if (prn_commsat[i] == presatinfo[1].m_Prn[j])
			{
				ele.push_back(presatinfo[1].Ele[j]);
				break;
			}
		}
	}
	for (int i = 0; i < all_const; i++)
	{
		if (i == 0)mk_zhy = 0;
		else mk_zhy += fre_comm[i - 1].size();
		for (int k = 0; k < fre_comm[i].size(); k++)
		{
			P1(mk_zhy + 1 + k, mk_zhy + 1 + k) = SQR((0.1) / (sin(ele[i])));// 2 * (SQR(0.003) / SQR(sin(ele[i])));SQR(gf_sigemaIon[prn_commsat[k]])
		}
	}
	int shuliang = comm_num;

	typedef struct xy
	{
		int x;
		int y;
	};

	vector<double> insigma;
	multiset<double> insigma_mset;
	vector<int> count_sit_min;
	vector<int> storesitu;

	int count_num_a = 0;
	int count_insigma = 0;

	int count_dowhile = 0;
	int controlB = 0;
	double listminvalue = 0.0;
	vector<double> dioncoe2;
	vector<vector<double>> dioncoe1;
	for (int i = 0; i < GetRlenth.size(); i++)
	{
		for (int j = 0; j < GetRlenth[i].size(); j++)
		{
			if (j > 0)
			{
				dioncoe2.push_back(((GetRfrequ[i][0] * GetRfrequ[i][0]) / (GetRfrequ[i][j] * GetRfrequ[i][j])));
			}
			else
			{
				dioncoe2.push_back(1.0);
			}
		}
		if (dioncoe2.size() > 0)
		{
			dioncoe1.push_back(dioncoe2);
			dioncoe2.clear();
		}
	}

	//估计ISB
	if (1)
	{
		double num_dtr = _SystemN + 1;
		double sigemais = 0.0;
		MatrixT Bis, Lis, Vis, Xis, Nbbis, Sigemais;
		Bis.Resize(constnu, num_dtr);
		Lis.Resize(constnu, 1);
		Vis.Resize(constnu, 1);
		Xis.Resize(num_dtr, 1);
		Nbbis.Resize(num_dtr, num_dtr);

		for (int i = 0; i < all_const; i++)
		{
			for (int j = 0; j < presatinfo[1].m_Size + 1; j++) {

				if (prn_commsat[i] == presatinfo[1].m_Prn[j])
				{
					satpos2[0] = presatinfo[1].m_Xs[3 * j + 0];
					satpos2[1] = presatinfo[1].m_Xs[3 * j + 1];
					satpos2[2] = presatinfo[1].m_Xs[3 * j + 2];

					for (int k = 0; k < 3; ++k)
					{
						rs2 += (RKpos_cur[k] - satpos2[k]) * (RKpos_cur[k] - satpos2[k]);
						xyz2[k] = RKpos_cur[k] - satpos2[k];																// (Xr-Xs)^2 +(Yr-Ys)^2+(Zr-Zs)^2
					}
					rs2 = sqrt(rs2);
					e2.e_x = xyz2[0] / rs2;
					e2.e_y = xyz2[1] / rs2;
					e2.e_z = xyz2[2] / rs2;
					e_2.push_back(e2);
					rs_2.push_back(rs2);
					rs2 = 0.0;
				}
			}
		}
		for (int m = 0; m < all_const; m++)
		{
			if (m == 0)mk_zhy = 0;
			else { mk_zhy += fre_comm[m - 1].size(); }
			for (int i = 0; i < fre_comm[m].size(); i++)
			{
				Lis(mk_zhy + i + 1, 1) = (L_Cur[m][i] - rs_2[m]) + (rs_1[m] - L_Pre[m][i]);
				if (fabs(Lis(mk_zhy + i + 1, 1)) > 20)
				{
					double a = Lis(mk_zhy + i + 1, 1);
					a = fabs(a);
				}
				Bis(mk_zhy + i + 1, 1) = 1.0;

				if (ISBco > 0)
				{
					switch (GetSystem(prn_commsat[m]))
					{
					case 'G':break;
						//case 'R':Bis(mk_zhy + i + 1, 1 + systemset[0]) = 1.0; break;
					case 'C':Bis(mk_zhy + i + 1, 1 + systemset[1]) = 1.0; break;
					case 'B':Bis(mk_zhy + i + 1, 1 + systemset[2]) = 1.0; break;
					case 'E':Bis(mk_zhy + i + 1, 1 + systemset[3]) = 1.0; break;
					default:break;
					}
				}
			}
		}
		
	
		LeastSquare(Bis, Lis, P1, Vis, Nbbis, Xis);

		Sigemais = Vis.Trans() * P1 * Vis;
		sigemais = Sigemais.GetValue(1, 1);
		sigemais = sqrt(sigemais / (constnu - _SystemN));
		rs_2.clear();

		ofstream ofisbs;
		//ofisbs.open("D:\\data10.0\\out\\isbs.out", ios::app);

		double xisb[5] = {};
		for (int i = 0; i < _SystemN - 1; i++)
		{
			xisb[i] = Xis.GetValue(i + 2, 1);
		}
		ofisbs << abs(xisb[0]) << '\t' << xisb[1] << '\t' << xisb[2] << '\t' << xisb[3] << '\t' << xisb[4] << endl;
	}


	vector<double> sV;
	vector<double> time1;
	double time11 = 0.0;

	int clip_map = 0;			
	int clip_range = 100;		
	vector<int> circ;			
	int porm = 2;				
	int cclip_value = 0;		
	int crow = 0, ccol = 0;
	vector<int> scclip_sit;
	vector<int> scclip_value;
	vector<vector<double>> recordsucc;
	vector<int> origisit; int countorigin = 0.0;
	circ.clear();
	vector<double> standV; multiset<double> standV_mset;
	vector<int> prnnumb; double dxdydz = 0.0;
	int countset = 0;

	ofxslip << flag_epoch_count - 1 << '\t' << constnu << '\t';

	clock_t start1, end1;


	//std::cout << "prnnumb.size " << prnnumb.size() << std::endl;
	//std::cout << "comm_num " << comm_num << std::endl;
	clock_t start3, end3;
	double time33 = 0.0;
	vector<double> time3;
	//std::cout << "circ.size " << circ.size() << std::endl;

//	cout << "  ITDE_OTDE " << ITDE_OTDE << endl;
	ofx2.open("D:\\data10.0\\outt\\slip.out", ios::app);
	ofx1.open("D:\\data10.0\\outt\\TDCPslip.csv", ios::trunc);
	if (ofx1.is_open())
		ofx1 << "epoch,prn,system,freq_index,signal,"
		"wavelength_m,residual_m,residual_cyc,slip_cycles_est,"
		"sigma0_m,threshold_m\n";
	
	if (0)
	{
		
		count_dowhile = 0;
		controlB = 0;
		storesitu.clear(); e_2.clear(); rs_2.clear();

		B1.Clear(); Nbb1.Clear(); X1.Clear(); L1.Clear(); V1.Clear(); VZ.Clear(); Pp1.Clear(); C1.Clear(); Ncc1.Clear(); v__1.Clear(); p__1.Clear(); Sigema0__1.Clear(); Pxx.Clear(); W1.Clear();
		L1.Resize(constnu, 1); V1.Resize(constnu, 1); Pp1.Resize(constnu, constnu); v__1.Resize(comm_num, 1); p__1.Resize(comm_num, comm_num); Sigema0__1.Resize(1, 1); Ncc1.Resize(all_const, all_const);

		do
		{
			listminvalue = 100.0;
			do
			{
				weizhishu = 4 + 1 + count_dowhile - controlB;//+ all_const;  
				B1.Resize(constnu, weizhishu);
				Nbb1.Resize(weizhishu, weizhishu);
				X1.Resize(weizhishu, 1);
				C1.Resize(weizhishu, weizhishu);
				Pxx.Resize(weizhishu, weizhishu);
				W1.Resize(weizhishu, 1);
				VZ.Resize(weizhishu, 1);

				if (controlB == 0)
				{
					storesitu.push_back(5 + count_dowhile);
				}

				for (int shuchu = 0; shuchu < shuliang; shuchu++)
				{
					for (int i = 0; i < all_const; i++)
					{
						for (int j = 0; j < presatinfo[1].m_Size + 1; j++) {
							if (prn_commsat[i] == presatinfo[1].m_Prn[j])
							{
								satpos2[0] = presatinfo[1].m_Xs[3 * j + 0];
								satpos2[1] = presatinfo[1].m_Xs[3 * j + 1];
								satpos2[2] = presatinfo[1].m_Xs[3 * j + 2];

								for (int k = 0; k < 3; ++k)
								{
									rs2 += (RKpos_cur[k] - satpos2[k]) * (RKpos_cur[k] - satpos2[k]);
									xyz2[k] = RKpos_cur[k] - satpos2[k];																// (Xr-Xs)^2 +(Yr-Ys)^2+(Zr-Zs)^2
								}
								rs2 = sqrt(rs2);
								e2.e_x = xyz2[0] / rs2;
								e2.e_y = xyz2[1] / rs2;
								e2.e_z = xyz2[2] / rs2;
								e_2.push_back(e2);
								rs_2.push_back(rs2);
								rs2 = 0.0;
							}
						}
					}
					for (int m = 0; m < all_const; m++)
					{
						if (m == 0)mk_zhy = 0;
						else { mk_zhy += fre_comm[m - 1].size(); }
						for (int i = 0; i < fre_comm[m].size(); i++)
						{
							L1(mk_zhy + i + 1, 1) = (L_Cur[m][i] - rs_2[m]) + (rs_1[m] - L_Pre[m][i]);// - gf_ion[prn_commsat[m]] * dioncoe[i];//#t+1-t
							double cycle_slip_value = (L_Cur[m][i] - rs_2[m]) + (rs_1[m] - L_Pre[m][i]);
							L1(mk_zhy + i + 1, 1) = cycle_slip_value;
							if (fabs(L1(mk_zhy + i + 1, 1)) > 20)
							{
								double a = L1(mk_zhy + i + 1, 1);
								a = fabs(a);
							}
							B1(mk_zhy + i + 1, 1) = (e_2[m].e_x);
							B1(mk_zhy + i + 1, 2) = (e_2[m].e_y);
							B1(mk_zhy + i + 1, 3) = (e_2[m].e_z);
							B1(mk_zhy + i + 1, 4) = 1;						//接收机钟差


							//B1(mk_zhy + i + 1, m + 5) = dioncoe1[m][i];
							char buff = GetSystem_GREC2C3(prn_commsat[m]);
							int a = _SystemIndex[buff];
							//B1(mk_zhy + i + 1, 3 + _SystemIndex[GetSystem_GREC2C3(prn_commsat[m])]) = 1.0;//接收机钟差 一个系统就是多一个； 
						}
					}
				
					for (int i = 0; i < weizhishu; i++)
					{
						C1(i + 1, i + 1) = 1.0;
					}
					if (controlB == 0)
					{
						for (int i = 0; i < constnu; i++)  //# add canshu dia
						{
							if (i == shuchu)
							{
								B1(i + 1, storesitu[count_dowhile]) = 1;
								break;
							}
						}
					}

					for (int i = 0; i < count_sit_min.size(); i++)
					{
						B1(count_sit_min[i] + 1, storesitu[i]) = 1;
					}

					if (LeastSquare(B1, L1, P1, V1, Nbb1, X1) == -1);    // LSQ  Nbb=(btpb)^-1
					//if (LeastSquareconstraint(B1, L1, P1, C1, W1, V1, VZ, Nbb1, Ncc1, X1, Pxx) == -1);//附加约束条件的最小二乘

					//B1.Output_to_File("E:\\B1.out", 1);
					//L1.Output_to_File("E:\\L1.out", 1);

					double sigema0__1 = 0;

					Sigema0 = (V1.Trans() * P1 * V1) + (VZ.Trans() * Pxx * VZ);
					sigema0 = Sigema0.GetValue(1, 1);
					sigema0 = sqrt(sigema0 / (constnu - 4));

					count_num_a = 0;
					if (count_sit_min.size() == 0)
					{
						insigma.push_back(sigema0);
						insigma_mset.insert(sigema0);
					}
					else
					{
						for (int i = 0; i < count_sit_min.size(); i++)
						{
							if (shuchu == count_sit_min[i])
							{
								count_num_a++;
							}
						}
						if (count_num_a != 0)
						{
							insigma.push_back(0.0);
						}
						else
						{
							insigma.push_back(sigema0);
							insigma_mset.insert(sigema0);
						}
					}

					if (controlB == 0)
					{
						for (int i = 0; i < constnu; i++)  // # add canshu dia
						{
							if (i == shuchu)
							{
								B1(i + 1, storesitu[count_dowhile]) = 0;
								break;
							}

						}
					}
				}//验后单位权中误差

				count_insigma = 0;
				if (controlB == 0)
				{
					for (int i = 0; i < insigma.size(); i++)
					{
						double Hypoth_test = (insigma[i] - popt.Zavemin) / (popt.Zsigma);
						if (Hypoth_test > 2.58 && insigma[i] < (3 * listminvalue))
						{
							count_insigma++;
						}
					}
				}
				if (controlB == 1)
				{
					double sigmafinal = 0;
					double sizein = insigma_mset.size();
					for (multiset<double>::iterator it = insigma_mset.begin(); it != insigma_mset.end(); it++)
					{
						sigmafinal += (*it);
					}
					sigmafinal = sigmafinal / sizein;
					ofx2 << sigmafinal << endl;
				}
				if (count_insigma != 0)
				{
					for (int i = 0; i < insigma.size(); i++)
					{
						if (insigma[i] == *(insigma_mset.begin()))
						{
							count_sit_min.push_back(i);
							listminvalue = *(insigma_mset.begin());
							break;
						}
					}
					count_dowhile++;
				}
				else
				{
					controlB++;
					if (controlB == 2)
					{
						X1.Output_to_File("C:\\Users\\15149\\Desktop\\shiy\\X11.out", 1);
						ofsigmaxpre << X1(1, 1) << '\t' << X1(2, 1) << '\t' << X1(3, 1) << '\t' << X1(4, 1) << endl;
					}
				}
				insigma_mset.clear();
				insigma.clear();

			} while (!count_insigma == 0);
		} while (controlB != 2);
		for (int j = 0; j < count_sit_min.size(); j++)
		{
			double temp = X1(5 + j, 1) / (GetRlenth[(count_sit_min[j] / frenumcontrol1)][(count_sit_min[j] % frenumcontrol1)]);
			if (abs(temp) >= 0.5)
			{
				origisit.push_back(count_sit_min[j]);//原始数据中包含的周跳位置
				countorigin += 1;//原始数据中包含的周跳个数
			}
		}
	}
	int countsize = 0;
	double percentl = 0.2;
	double percentr = 0.3;
	for (int i = 0; i < (comm_num - prnnumb.size()); i++)
	{
		double a = ((i + 1.0) / comm_num);
		if (a >= percentl && a < percentr)
			circ.push_back(0);
	}
	std::cout << "DEBUG comm_num=" << comm_num
		<< " all_const=" << all_const
		<< " frenumcontrol1=" << frenumcontrol1
		<< " prnnumb.size=" << prnnumb.size()
		<< std::endl;

	static int total_TP = 0, total_FN = 0, total_FP = 0, total_TN = 0;
	static int total_detect_success = 0;
	static int total_repair_success = 0;
	static int total_detection_and_repair_success = 0;
	static int total_slips_all = 0;
	static int total_repair_samples = 0;
	static double total_repair_error = 0.0;
	static int total_observations = 0;
	static int flag_epoch_count = 0;
	static double previous_rclk = 0.0;
	static vector<double> previous_isb(4, 0.0);

	std::cout << "circ.size " << circ.size() << std::endl;
	std::cout << "comm_num=" << comm_num << ", frenumcontrol1=" << frenumcontrol1 << std::endl;

	vector<vector<double>> L_Curtemp = L_Curtt;

	// ── 抗差最小二乘迭代收敛 ─────────────────────────────────────────────
	if (1)
	{
		for (int p = 0; p < circ.size(); p++)
			std::cout << "Epoch: " << flag_epoch_count << ", p=" << p
			<< ", circ[p]=" << circ[p] << ", countorigin=" << countorigin << std::endl;

		for (int p = 0; p < circ.size(); p++)
		{
			if (circ[p] > countorigin)
			{
				scclip_sit.clear(); scclip_value.clear();
				for (int i = 0; i < countorigin; i++)
				{
					clip_map = origisit[i];
					crow = clip_map / frenumcontrol1;
					ccol = clip_map % frenumcontrol1;
					if ((rand() % porm) == 0)
					{
						cclip_value = rand() % clip_range + 1;
						L_Curtemp[crow][ccol] -= cclip_value * GetRlenth[crow][ccol];
						scclip_sit.push_back(clip_map);
						scclip_value.push_back(-cclip_value);
						ofxslip << scclip_sit[i] + 1 << '!' << scclip_value[i] << '\t';
					}
					else
					{
						cclip_value = rand() % clip_range + 1;
						L_Curtemp[crow][ccol] += cclip_value * GetRlenth[crow][ccol];
						scclip_sit.push_back(clip_map);
						scclip_value.push_back(cclip_value);
						ofxslip << scclip_sit[i] + 1 << '!' << scclip_value[i] << '\t';
					}
				}
				int i1 = 0, counta = 0;
				do {
					counta = 0; int temp = 0;
					do
					{
						temp = 0;
						clip_map = rand() % comm_num;
						for (int i = 0; i < prnnumb.size(); i++)
							if (clip_map == prnnumb[i]) temp++;
					} while (temp);
					crow = clip_map / frenumcontrol1;
					ccol = clip_map % frenumcontrol1;
					for (int j = 0; j < scclip_sit.size(); j++)
					{
						if (clip_map == scclip_sit[j]) { counta++; break; }
					}
					if (counta == 0)
					{
						if ((rand() % porm) == 0)
						{
							cclip_value = rand() % clip_range + 1;
							L_Curtemp[crow][ccol] -= cclip_value * GetRlenth[crow][ccol];
							scclip_sit.push_back(clip_map);
							scclip_value.push_back(-cclip_value);
							ofxslip << scclip_sit[countorigin + i1] + 1 << '!' << scclip_value[countorigin + i1] << '\t';
						}
						else
						{
							cclip_value = rand() % clip_range + 1;
							L_Curtemp[crow][ccol] += cclip_value * GetRlenth[crow][ccol];
							scclip_sit.push_back(clip_map);
							scclip_value.push_back(cclip_value);
							ofxslip << scclip_sit[countorigin + i1] + 1 << '!' << scclip_value[countorigin + i1] << '\t';
						}
						i1++;
					}
				} while (countorigin + i1 < circ[p]);
			}
			else
			{
				countsize++;
			}

			std::cout << "]=" << scclip_sit.size() << std::endl;

			flag_epoch_count++;
		

			// ── 抗差最小二乘迭代收敛（TDCP + IGGIII）──
			if (1)
			{
				double recorddx0 = 0.0, recorddx1 = 0.0, re = 0.0;
				double dxdydz = 0.0;
				int countset = 0;
				double prev_sigema0 = 0.0;
				start1 = clock();

				

				weizhishu = 4 + ISBco;
				B1.Clear(); Nbb1.Clear(); X1.Clear(); L1.Clear();
				V1.Clear(); VZ.Clear(); C1.Clear(); Pxx.Clear(); W1.Clear(); Qvv.Clear();
				L1.Resize(constnu, 1); V1.Resize(constnu, 1);
				Qvv.Resize(constnu, constnu);
				B1.Resize(constnu, weizhishu);
				Nbb1.Resize(weizhishu, weizhishu);
				X1.Resize(weizhishu, 1);

				do
				{
					countset++;
					recorddx0 = recorddx1;
					standV.clear();
					standV_mset.clear();
					e_2.clear();
					rs_2.clear();

					double rs2 = 0.0, xyz2[3] = { 0.0 };
					for (int i = 0; i < all_const; i++)
					{
						rs2 = 0.0;
						for (int j = 0; j < presatinfo[1].m_Size + 1; j++)
						{
							if (prn_commsat[i] == presatinfo[1].m_Prn[j])
							{
								satpos2[0] = presatinfo[1].m_Xs[3 * j + 0];
								satpos2[1] = presatinfo[1].m_Xs[3 * j + 1];
								satpos2[2] = presatinfo[1].m_Xs[3 * j + 2];
								for (int k = 0; k < 3; ++k)
								{
									rs2 += (RKpos_cur[k] - satpos2[k]) * (RKpos_cur[k] - satpos2[k]);
									xyz2[k] = RKpos_cur[k] - satpos2[k];
								}
								rs2 = sqrt(rs2);
								e2.e_x = xyz2[0] / rs2;
								e2.e_y = xyz2[1] / rs2;
								e2.e_z = xyz2[2] / rs2;
								e_2.push_back(e2);
								rs_2.push_back(rs2);
								break;
							}
						}
					}

					for (int m = 0; m < all_const; m++)
					{
						if (m == 0) mk_zhy = 0;
						else        mk_zhy += (int)fre_comm[m - 1].size();

						char sys = GetSystem(prn_commsat[m]);
						
						if (ISBco > 0)
						{
							int sys_col = -1;
							switch (sys)
							{
							case 'G': sys_col = 4;                 break;
							case 'R': sys_col = 4 + systemset[0]; break;
							case 'C': sys_col = 4 + systemset[1]; break;
							case 'B': sys_col = 4 + systemset[2]; break;
							case 'E': sys_col = 4 + systemset[3]; break;
							default: break;
							}
							
						}
						

						for (int i = 0; i < (int)fre_comm[m].size(); i++)
						{
							L1(mk_zhy + i + 1, 1) =
								(L_Curtemp[m][i] - rs_2[m]) + (rs_1[m] - L_Pre[m][i]);

							B1(mk_zhy + i + 1, 1) = e_2[m].e_x;
							B1(mk_zhy + i + 1, 2) = e_2[m].e_y;
							B1(mk_zhy + i + 1, 3) = e_2[m].e_z;
							B1(mk_zhy + i + 1, 4) = 1;
							if (ISBco > 0)
							{
								switch (sys)
								{
								case 'R': B1(mk_zhy + i + 1, 4 + systemset[0]) = 1.0; break;
								case 'C': B1(mk_zhy + i + 1, 4 + systemset[1]) = 1.0; break;
								case 'B': B1(mk_zhy + i + 1, 4 + systemset[2]) = 1.0; break;
								case 'E': B1(mk_zhy + i + 1, 4 + systemset[3]) = 1.0; break;
								default: break;
								}
							}
						}
					}

					int lsq_result = LeastSquareQcc(B1, L1, P1, V1, Nbb1, Qvv, X1);
					

					Sigema0 = (V1.Trans() * P1 * V1);
					sigema0 = sqrt(Sigema0.GetValue(1, 1) / (constnu - 4));

					for (int i = 0; i < constnu; i++)
					{
						standV.push_back(fabs(V1(i + 1, 1) / (sigema0 * sqrt(Qvv(i + 1, i + 1)))));
						standV_mset.insert(fabs(V1(i + 1, 1) / (sigema0 * sqrt(Qvv(i + 1, i + 1)))));
						sV.push_back(V1(i + 1, 1));
					}

					// IGGIII 降权
					double k0 = 0.8, k1 = 1.5;
					for (int i = 0; i < constnu; i++)
					{
						if (standV[i] > k1)
						{
							P1(i + 1, i + 1) = 1e-5;
							int tempset = 0;
							for (int j = 0; j < (int)prnnumb.size(); j++)
								if (i == prnnumb[j]) tempset++;
							if (tempset == 0) prnnumb.push_back(i);
						}
						else if (standV[i] > k0)
						{
							P1(i + 1, i + 1) *= (k0 / standV[i])
								* ((k1 - standV[i]) / (k1 - k0))
								* ((k1 - standV[i]) / (k1 - k0));
						}
					}

					dxdydz = sqrt(X1(1, 1) * X1(1, 1) + X1(2, 1) * X1(2, 1) + X1(3, 1) * X1(3, 1));
					recorddx1 = dxdydz;
					re = fabs(recorddx0 - recorddx1);
					if (re < 1e-5) ofx3 << sigema0 << endl;
					if (re < 1e-5 && fabs(sigema0 - prev_sigema0) < 0.01)
					{
						
						break;
					}
					recorddx0 = recorddx1;
					prev_sigema0 = sigema0;
				} while (!(re < 1e-4));
			}

			double delta_X = X1(1, 1);
			double delta_Y = X1(2, 1);
			double delta_Z = X1(3, 1);
			double rclk_drift = X1(4, 1);


			// ── SDED（星间差分历元差分）+ TDCP 无星间差分输出 ──
			if (flag_epoch_count >= 2)
			{
				const double sample_interval = 15.0;

				// 选参考星（仰角最高）
				int ref_m = -1; double max_ele = -1.0;
				for (int m = 0; m < all_const; m++)
					if (ele[m] > max_ele) { max_ele = ele[m]; ref_m = m; }

				if (ref_m < 0)
				{
					//cout << "警告: 未找到参考星" << endl;
				}
				else
				{
					
					auto corrected_tdcp = [&](int m, int freq_i, int obs_idx,
						double& out_ion, double& out_trop, double& out_satclk) -> double
						{
							double raw = (L_Curtemp[m][freq_i] - rs_2[m]) + (rs_1[m] - L_Pre[m][freq_i]);
							out_ion = 0.0;
							if (obs_idx < (int)ion.size() && obs_idx < (int)ion_prev.size())
								out_ion = ion[obs_idx] - ion_prev[obs_idx];
							out_trop = 0.0;
							if (m < (int)trop.size() && m < (int)trop_prev.size())
								out_trop = trop[m] - trop_prev[m];
							unsigned int prn = prn_commsat[m];
							double clk_cur = 0.0, clk_prev = 0.0;
							for (size_t idx = 0; idx < satclks_prn.size(); ++idx)
								if (satclks_prn[idx] == prn) { clk_cur = satclks[idx]; break; }
							for (size_t idx = 0; idx < satclks_prn_prev.size(); ++idx)
								if (satclks_prn_prev[idx] == prn) { clk_prev = satclks_prev[idx]; break; }
							out_satclk = (clk_cur - clk_prev) * LIGHTSPEED;
							return raw - out_ion - out_trop - out_satclk;
						};

				
					int ref_mk = 0;
					for (int mm = 0; mm < ref_m; mm++) ref_mk += (int)fre_comm[mm].size();

					if (rs_1[ref_m] < 1e6 || rs_2[ref_m] < 1e6)
					{
						//cout << "警告: 参考星几何距离无效，跳过本历元SDED计算" << endl;
					}
					else
					{
						double ref_ion = 0.0, ref_trop = 0.0, ref_satclk = 0.0;
						double ref_tdcp_corr = corrected_tdcp(ref_m, 0, ref_mk, ref_ion, ref_trop, ref_satclk);
						double ref_proj = e_2[ref_m].e_x * delta_X
							+ e_2[ref_m].e_y * delta_Y
							+ e_2[ref_m].e_z * delta_Z;
						

					
						int mk_zhy_loop = 0;
						for (int m = 0; m < all_const; m++)
						{
							if (m == 0) mk_zhy_loop = 0;
							else        mk_zhy_loop += (int)fre_comm[m - 1].size();
							if (m == ref_m) continue;
							if (rs_1[m] < 1e6 || rs_2[m] < 1e6) continue;

							double cur_proj = e_2[m].e_x * delta_X + e_2[m].e_y * delta_Y + e_2[m].e_z * delta_Z;
							double proj_diff = cur_proj - ref_proj;

							for (int i = 0; i < (int)fre_comm[m].size(); i++)
							{
								int obs_idx = mk_zhy_loop + i;
								double cur_ion = 0.0, cur_trop = 0.0, cur_satclk = 0.0;
								corrected_tdcp(m, i, obs_idx, cur_ion, cur_trop, cur_satclk);
							}
						}
					}
				}
				ofx2.flush();

				{
					int mk_t = 0;
					for (int m = 0; m < all_const; m++)
					{
						if (m == 0) mk_t = 0;
						else        mk_t += (int)fre_comm[m - 1].size();
						if (rs_1[m] < 1e6 || rs_2[m] < 1e6) continue;

						double proj_t = e_2[m].e_x * delta_X + e_2[m].e_y * delta_Y + e_2[m].e_z * delta_Z;

						for (int i = 0; i < (int)fre_comm[m].size(); i++)
						{
							int obs_t = mk_t + i;
							double raw_t = (L_Curtemp[m][i] - rs_2[m]) + (rs_1[m] - L_Pre[m][i]);
							double trop_t = 0.0;
							if (m < (int)trop.size() && m < (int)trop_prev.size())
								trop_t = trop[m] - trop_prev[m];
							unsigned int prn_t = prn_commsat[m];
							double clk_c = 0.0, clk_p = 0.0;
							for (size_t idx = 0; idx < satclks_prn.size(); ++idx)
								if (satclks_prn[idx] == prn_t) { clk_c = satclks[idx]; break; }
							for (size_t idx = 0; idx < satclks_prn_prev.size(); ++idx)
								if (satclks_prn_prev[idx] == prn_t) { clk_p = satclks_prev[idx]; break; }
							double satclk_t = (clk_c - clk_p) * LIGHTSPEED;
							double corr_t = raw_t - trop_t - satclk_t;
							double wl_t = GetRlenth[m][i];
							double feat_t = (wl_t > 1e-9 && sample_interval > 1e-9)
								? (corr_t / wl_t / sample_interval) : 0.0;
							double slip_t = 0.0;
							for (int k = 0; k < (int)scclip_sit.size(); k++)
								if (scclip_sit[k] / frenumcontrol1 == m && scclip_sit[k] % frenumcontrol1 == i)
								{
									slip_t = scclip_value[k]; break;
								}

							ofx2 << "Epoch: " << flag_epoch_count
								<< ", PRN=" << prn_commsat[m]
								<< ", freq_index=" << i
								<< ", wavelength=" << wl_t
								<< ", TDCP cycle_slip=" << feat_t
								<< ", Simulated slip=" << slip_t
								<< ", iondelay=" << 0.0
								<< ", tropdelay=" << trop_t
								<< ", satclk_delay=" << satclk_t
								<< ", clock_drift=" << rclk_drift << "\n";
						}
					}
					ofx2.flush();
				}
			} // end if (flag_epoch_count >= 2)

			ion_prev = ion;
			satclks_prev = satclks;
			trop_prev = trop;
			previous_rclk = popt.Rclk;
			satclks_prn_prev = satclks_prn;
			ofxorigincount << prnnumb.size() << '\t';
			vector<double> repairvalue;
			if (1)
			{
				storesitu.clear();
				count_sit_min.clear();
				std::cout << " prnnumb.size() " << prnnumb.size() << std::endl;

				for (int t = 0; t < (int)prnnumb.size(); t++)
				{
					weizhishu = 4 + ISBco + 1;
					B1.Clear(); Nbb1.Clear(); X1.Clear(); L1.Clear(); V1.Clear();
					L1.Resize(constnu, 1); V1.Resize(constnu, 1);
					B1.Resize(constnu, weizhishu);
					Nbb1.Resize(weizhishu, weizhishu);
					X1.Resize(weizhishu, 1);

					for (int m = 0; m < all_const; m++)
					{
						if (m == 0) mk_zhy = 0;
						else        mk_zhy += (int)fre_comm[m - 1].size();
						char sys = GetSystem(prn_commsat[m]);

						for (int i = 0; i < (int)fre_comm[m].size(); i++)
						{
							L1(mk_zhy + i + 1, 1) =
								(L_Curtemp[m][i] - rs_2[m]) + (rs_1[m] - L_Pre[m][i]);
							B1(mk_zhy + i + 1, 1) = e_2[m].e_x;
							B1(mk_zhy + i + 1, 2) = e_2[m].e_y;
							B1(mk_zhy + i + 1, 3) = e_2[m].e_z;
							B1(mk_zhy + i + 1, 4) = 1.0;
							if (ISBco > 0)
							{
								switch (sys)
								{
								case 'R': B1(mk_zhy + i + 1, 4 + systemset[0]) = 1.0; break;
								case 'C': B1(mk_zhy + i + 1, 4 + systemset[1]) = 1.0; break;
								case 'B': B1(mk_zhy + i + 1, 4 + systemset[2]) = 1.0; break;
								case 'E': B1(mk_zhy + i + 1, 4 + systemset[3]) = 1.0; break;
								default: break;
								}
							}
						}
					}

					P1(prnnumb[t] + 1, prnnumb[t] + 1) =
						SQR(0.1 / sin(ele[prnnumb[t] / frenumcontrol1]));
					B1(prnnumb[t] + 1, weizhishu) = 1.0;

					int lsq_result = LeastSquareQcc(B1, L1, P1, V1, Nbb1, Qvv, X1);
					

					double tempx = X1(weizhishu, 1);
					L_Curtemp[prnnumb[t] / frenumcontrol1][prnnumb[t] % frenumcontrol1] -= tempx;
					std::cout << "Debug: repairvalue[" << t << "] = " << tempx << std::endl;
					repairvalue.push_back(tempx);
					ofXX << tempx << '\t';
					ofsitu << tempx << '\t';
				}

				end1 = clock();
				time1.push_back(double(end1 - start1) / CLOCKS_PER_SEC);
				ofx1.open("D:\\data10.0\\outt\\1.out", ios::app);
				ofrecord.open("D:\\data10.0\\outt\\record.out", ios::app);

				// ── 性能统计 ──
				int TP = 0, FN = 0, FP = 0, TN = 0;
				int detect_success_count = 0, repair_success_count = 0;
				int detection_and_repair_success_count = 0;
				double epoch_repair_error = 0.0;
				int repair_count = 0;

				// 清理非法 prnnumb
				for (int t = 0; t < (int)prnnumb.size(); )
				{
					int idx = prnnumb[t];
					if (idx < 0 || idx >= comm_num)
					{
						prnnumb.erase(prnnumb.begin() + t);
						if ((int)repairvalue.size() > t) repairvalue.erase(repairvalue.begin() + t);
					}
					else t++;
				}
				if (repairvalue.size() != prnnumb.size())
					repairvalue.resize(prnnumb.size(), 0.0);

				std::set<int> detected_set(prnnumb.begin(), prnnumb.end());
				for (int i = 0; i < (int)scclip_sit.size(); i++)
				{
					int true_idx = scclip_sit[i];
					if (true_idx < 0 || true_idx >= comm_num) continue;
					double true_slip = scclip_value[i];
					if (detected_set.count(true_idx))
					{
						TP++; detect_success_count++;
						for (int j = 0; j < (int)prnnumb.size(); j++)
						{
							if (prnnumb[j] == true_idx)
							{
								int prn = true_idx / frenumcontrol1;
								int freq = true_idx % frenumcontrol1;
								if (prn >= 0 && prn < (int)GetRlenth.size() &&
									freq >= 0 && freq < (int)GetRlenth[prn].size())
								{
									double wavelength = GetRlenth[prn][freq];
									double repaired_cyc = repairvalue[j] / wavelength;
									double error = fabs(repaired_cyc - true_slip);
									epoch_repair_error += error;
									repair_count++;
									if (error < 0.5)
									{
										repair_success_count++;
										detection_and_repair_success_count++;
									}
								}
								break;
							}
						}
					}
					else FN++;
				}

				FP = 0;
				for (int j = 0; j < (int)prnnumb.size(); j++)
				{
					bool is_real = false;
					for (int i = 0; i < (int)scclip_sit.size(); i++)
						if (scclip_sit[i] == prnnumb[j]) { is_real = true; break; }
					if (!is_real) FP++;
				}
				TN = comm_num - TP - FP - FN;

				total_TP += TP; total_FN += FN; total_FP += FP; total_TN += TN;
				total_detect_success += detect_success_count;
				total_repair_success += repair_success_count;
				total_detection_and_repair_success += detection_and_repair_success_count;
				total_slips_all += (int)scclip_sit.size();
				total_repair_error += epoch_repair_error;
				total_repair_samples += repair_count;
				total_observations += comm_num;

				double recall = (TP + FN) > 0 ? (double)TP / (TP + FN) : 0.0;
				double precision = (TP + FP) > 0 ? (double)TP / (TP + FP) : 0.0;
				double f1 = (precision + recall) > 0
					? 2 * precision * recall / (precision + recall) : 0.0;
				double mae_cycles = repair_count > 0 ? epoch_repair_error / repair_count : 0.0;
				double detection_and_repair_rate = scclip_sit.size() > 0
					? (double)detection_and_repair_success_count / scclip_sit.size() * 100 : 0.0;
				double accuracy = comm_num > 0 ? (double)(TP + TN) / comm_num * 100 : 0.0;

				std::ofstream ofsuccess("D:\\data10.0\\outt\\success_rates.csv", std::ios::app);
				ofsuccess << "历元: " << flag_epoch_count << ", 总模拟周跳数: " << scclip_sit.size() << "\n";
				ofsuccess << "总探测成功率: "
					<< (scclip_sit.size() > 0 ? detect_success_count * 100.0 / scclip_sit.size() : 0)
					<< "% (" << detect_success_count << "/" << scclip_sit.size() << ")\n";
				ofsuccess << "总修复成功率: "
					<< (repair_count > 0 ? repair_success_count * 100.0 / repair_count : 0)
					<< "% (" << repair_success_count << "/" << repair_count << ")\n";
				ofsuccess << "探测与修复总成功率: " << detection_and_repair_rate
					<< "% (" << detection_and_repair_success_count << "/" << scclip_sit.size() << ")\n";
				ofsuccess << "周跳探测准确率: " << accuracy
					<< "% (" << (TP + TN) << "/" << comm_num << ")\n";
				ofsuccess << "召回率: " << recall * 100 << "% | 精确率: " << precision * 100
					<< "% | F1: " << f1 * 100 << "%\n";
				ofsuccess << "修复误差MAE: " << mae_cycles << " 周\n";
				ofsuccess << "--------------------------------\n";
				ofsuccess.close();
			}

			// ── 所有历元最终统计 ──
			std::ofstream oftotal("D:\\data10.0\\outt\\TOTAL_SUCCESS_RATE.csv", std::ios::trunc);
			double final_detect_rate = total_slips_all > 0
				? (double)total_detect_success / total_slips_all * 100 : 0.0;
			double final_repair_rate = total_repair_samples > 0
				? (double)total_repair_success / total_repair_samples * 100 : 0.0;
			double final_dr_rate = total_slips_all > 0
				? (double)total_detection_and_repair_success / total_slips_all * 100 : 0.0;
			double final_accuracy = total_observations > 0
				? (double)(total_TP + total_TN) / total_observations * 100 : 0.0;
			double final_recall = (total_TP + total_FN) > 0
				? (double)total_TP / (total_TP + total_FN) : 0.0;
			double final_precision = (total_TP + total_FP) > 0
				? (double)total_TP / (total_TP + total_FP) : 0.0;
			double final_f1 = (final_precision + final_recall) > 0
				? 2 * final_precision * final_recall / (final_precision + final_recall) : 0.0;
			double final_mae = total_repair_samples > 0
				? total_repair_error / total_repair_samples : 0.0;

			oftotal << "============== 所有历元最终统计 ==============\n";
			oftotal << "总模拟周跳数: " << total_slips_all << "\n";
			oftotal << "总观测数: " << total_observations << "\n";
			oftotal << "总探测成功率: " << final_detect_rate << "% (" << total_detect_success << "/" << total_slips_all << ")\n";
			oftotal << "总修复成功率: " << final_repair_rate << "% (" << total_repair_success << "/" << total_repair_samples << ")\n";
			oftotal << "探测与修复总成功率: " << final_dr_rate << "% (" << total_detection_and_repair_success << "/" << total_slips_all << ")\n";
			oftotal << "周跳探测准确率: " << final_accuracy << "% (" << (total_TP + total_TN) << "/" << total_observations << ")\n";
			oftotal << "召回率(Recall): " << final_recall * 100 << "%\n";
			oftotal << "精确率(Precision): " << final_precision * 100 << "%\n";
			oftotal << "F1-Score: " << final_f1 * 100 << "%\n";
			oftotal << "修复误差MAE: " << final_mae << " 周\n";
			oftotal << "============================================\n";
			oftotal.close();
			

			prnnumb.clear();
			repairvalue.clear();
			scclip_sit.clear();
			scclip_value.clear();
			
			return 0;
		}
	}
}
	


int LeastSquarePPP::Cycle_slipTDCP(ppp_option_t& popt, vector <zhydata>& prepppdata, vector <zhyinfo>& presatinfo, PPPObsData& pppdata, PreciseData& predata, vector<double>& ion_gimsub, vector<int>& commprn)
{
	vector<double> L_P, L_C; vector<vector<double>>L_Pre, L_Cur;
	vector<double> P_P, P_C; vector<vector<double>>P_Pre, P_Cur;
	vector<vector<string>>fre_comm; vector<string> fre;
	vector<int>obstype; vector< vector<int>>Obstype;
	vector<unsigned int> prn_commsat;
	vector<double> rs_1;
	vector<double> rs_2;
	vector<unsigned int>prn_all; vector<int>typeint_all; vector<string>fre_all;
	int YQ = 0; int comm_num = 0;
	typedef struct zhy
	{
		double e_x;
		double e_y;
		double e_z;
	};
	zhy e1, e2;
	vector<zhy>e_1;
	vector<zhy>e_2;
	vector<double>ele;
	double L_pre = 0, L_cur = 0; double P_pre = 0, P_cur = 0; int m = 0;
	double rs1 = 0.0, rs2 = 0.0;
	double xyz1[3] = { 0 }, xyz2[3] = { 0 };
	double satpos1[3] = { 0 }, satpos2[3] = { 0 };
	double RKpos_cur[3] = { 0 };
	double RKpos_pre[3] = { 0 };
	int num = 0;
	num = 11;

	int epoch01 = 7200;//周跳历元
	int starsitu = 5;

	for (int i = 0; i < 3; i++)
	{
		RKpos_cur[i] = CoarseXYZ[i];
		RKpos_pre[i] = CoarseXYZ_pre[i];
	}
	vector<double>rand_num;
	vector<double> getRlenth1, getRfrequ1;
	vector<vector <double>> GetRlenth, GetRfrequ;
	double getRlenth = 0.0, getRfrequ = 0.0;
	int frenumcontrol[5] = { 1,1,1,1,1 };	
	int frenumcontrol1 = popt.frennn;	
	if (frenumcontrol1 == 2) { frenumcontrol[1] = 2; }
	if (frenumcontrol1 == 3) { frenumcontrol[1] = 2; frenumcontrol[2] = 3; }
	int OTDorITD = (popt.schemeset - 3);//0:OTD;1:ITD

	double percentl = 0.00;		
	double percentr = 0.00;		
	double dtt = 0.0;
	
	double dionsigma = (1.0 / 0.09);//0.01; // (1.0 / 0.09);
	double kkkcal = 3.0;


	double Zxave = 0.0, Zavemin = 0.0;
	popt.Zavemin = 0.0;
	for (int i = 0; i < ion_gimsub.size(); i++)
	{
		Zxave += (ion_gimsub[i] / ion_gimsub.size());
	}
	for (int i = 0; i < ion_gimsub.size(); i++)
	{
		Zavemin += ((ion_gimsub[i] - Zxave) * (ion_gimsub[i] - Zxave) / ion_gimsub.size());
	}
	if (abs(Zavemin) <= 0.01) { Zavemin = kkkcal * sqrt(abs(Zavemin)); }
	else if ((abs(Zavemin) > 0.01) && (abs(Zavemin) <= 0.029)) { Zavemin = (kkkcal + 0.2) * sqrt(abs(Zavemin)); }
	else { Zavemin = (kkkcal + 0.5) * sqrt(abs(Zavemin)); }

	popt.Zavemin = abs(Zavemin - (popt.Zsigma * 2.58));

	for (int i = 0; i < prepppdata[0].m_SatCount; i++)
	{
		for (int j = 0; j < prepppdata[1].m_SatCount; j++)
		{
			if (prepppdata[0].m_Prn[i] == prepppdata[1].m_Prn[j] && prepppdata[1].m_Prn[j] < 400)
			{
				int k = 0;
				for (k = 0; k < _Cycle_Slip.size(); ++k)
				{
					if (prepppdata[1].m_Prn[j] == _Cycle_Slip[k])//i=2 k=1
						break;
				}
				if (k == _Cycle_Slip.size())
				{
					unsigned int numb = prepppdata[1].m_Prn[j];
					prn_commsat.push_back(numb);
					for (int k = 0; k < type_int[j].size(); k++)//当前历元卫星频率数量
					{
						for (int l = 0; l < type_int_pre[i].size(); l++)
						{
							if (type_int[j][k] == frenumcontrol[0] && type_int_pre[i][l] == frenumcontrol[0])
							{
								L_cur = prepppdata[1].m_L1[j] * get_wth(numb, 1, predata.m_NavData);
								L_pre = prepppdata[0].m_L1[i] * get_wth(numb, 1, predata.m_NavData);
								P_cur = prepppdata[1].m_P1[j];
								P_pre = prepppdata[0].m_P1[i];
								getRlenth = get_wth(numb, 1, predata.m_NavData);
								getRfrequ = 299792458.0 / getRlenth;
							}
							else if (type_int[j][k] == frenumcontrol[1] && type_int_pre[i][l] == frenumcontrol[1])//第二个频率 
							{
								L_cur = prepppdata[1].m_L2[j] * get_wth(numb, 2, predata.m_NavData);
								L_pre = prepppdata[0].m_L2[i] * get_wth(numb, 2, predata.m_NavData);
								P_cur = prepppdata[1].m_P2[j];
								P_pre = prepppdata[0].m_P2[i];
								getRlenth = get_wth(numb, 2, predata.m_NavData);
								getRfrequ = 299792458.0 / getRlenth;
							}
							else if (type_int[j][k] == frenumcontrol[2] && type_int_pre[i][l] == frenumcontrol[2])
							{
								L_cur = prepppdata[1].m_L3[j] * get_wth(numb, 3, predata.m_NavData);
								L_pre = prepppdata[0].m_L3[i] * get_wth(numb, 3, predata.m_NavData);
								P_cur = prepppdata[1].m_P3[j];
								P_pre = prepppdata[0].m_P3[i];
								getRlenth = get_wth(numb, 3, predata.m_NavData);
								getRfrequ = 299792458.0 / getRlenth;
							}
							else if (type_int[j][k] == frenumcontrol[3] && type_int_pre[i][l] == frenumcontrol[3])
							{
								L_cur = prepppdata[1].m_L4[j] * get_wth(numb, 4, predata.m_NavData);
								L_pre = prepppdata[0].m_L4[i] * get_wth(numb, 4, predata.m_NavData);
								P_cur = prepppdata[1].m_P4[j];
								P_pre = prepppdata[0].m_P4[i];
							}
							else if (type_int[j][k] == frenumcontrol[4] && type_int_pre[i][l] == frenumcontrol[4])
							{
								//if (flag_epoch_count == epoch01 && prn_commsat.size() > starsitu)//flag_epoch_count  代表历元数 
								L_cur = prepppdata[1].m_L5[j] * get_wth(numb, 5, predata.m_NavData);
								L_pre = prepppdata[0].m_L5[i] * get_wth(numb, 5, predata.m_NavData);
								P_cur = prepppdata[1].m_P5[j];
								P_pre = prepppdata[0].m_P5[i];
							}
							if (L_pre != 0) {
								prn_all.push_back(numb); typeint_all.push_back(type_int[j][k]); fre_all.push_back(type_[j][k + 1]);
								L_P.push_back(L_pre); L_C.push_back(L_cur); L_pre = 0;
								P_P.push_back(P_pre); P_C.push_back(P_cur);
								fre.push_back(type_[j][k + 1]); obstype.push_back(type_int[j][k]);
								getRlenth1.push_back(getRlenth);
								getRfrequ1.push_back(getRfrequ);
								break;
							}
						}
					}
					fre_comm.push_back(fre);
					L_Pre.push_back(L_P); L_Cur.push_back(L_C);
					P_Pre.push_back(P_P); P_Cur.push_back(P_C);
					Obstype.push_back(obstype);
					comm_num += fre.size();
					GetRlenth.push_back(getRlenth1);
					GetRfrequ.push_back(getRfrequ1);
					fre.clear(); L_P.clear(); L_C.clear(); obstype.clear(); P_P.clear(); P_C.clear(); getRlenth1.clear(); getRfrequ1.clear();
					break;
				}
			}
		}
	}//确定前后历元相同卫星

	for (int i = 0; i < Obstype.size(); i++)
	{
		if (Obstype[i].size() != frenumcontrol1)
		{
			return 0;
		}
	}

#if (0)
	for (int i = 0; i < prn_commsat.size(); i++)
	{
		QC_Scdia(L_Pre[i], L_Cur[i], P_Pre[i], P_Cur[i], Obstype[i], prn_commsat[i], predata, prepppdata[0]);
	}
#endif 


	//if (flag_epoch_count == epoch01)
	//{
	//	for (int j = 0; j < prn_commsat.size(); j++)
	//	{
	//		cout << prn_commsat[j] << endl;
	//	}
	//}

	for (int i = 0; i < prn_commsat.size(); i++)
	{
		for (int j = 0; j < presatinfo[0].m_Size + 1; j++) {
			if (prn_commsat[i] == presatinfo[0].m_Prn[j] && presatinfo[0].m_Prn[j] < 400) {
				satpos1[0] = presatinfo[0].m_Xs[3 * j + 0];
				satpos1[1] = presatinfo[0].m_Xs[3 * j + 1];
				satpos1[2] = presatinfo[0].m_Xs[3 * j + 2];
				for (int k = 0; k < 3; ++k)
				{
					rs1 += SQR(RKpos_pre[k] - satpos1[k]);// (Xr-Xs)^2 +(Yr-Ys)^2+(Zr-Zs)^2
				}
				rs1 = sqrt(rs1);

				rs_1.push_back(rs1);
				rs1 = 0.0;
			}
		}
	}

	double sigema0 = 0.0;
	MatrixT B1, P1, L1, V1, VZ, Nbb1, X1, Pp1, Sigema0, Sigema1, C1, W1, Ncc1, Qvv, Pxx;
	MatrixT v__1(comm_num, 1);
	MatrixT p__1(comm_num, comm_num); MatrixT Sigema0__1(1, 1);

	int constnu, all_const, weizhishu;
	int mk_zhy = 0;

	all_const = prn_commsat.size();
	constnu = comm_num;

	weizhishu = 4 + all_const;   
	B1.Resize(constnu, weizhishu);
	Nbb1.Resize(weizhishu, weizhishu);
	X1.Resize(weizhishu, 1);

	P1.Resize(constnu, constnu);
	Qvv.Resize(constnu, constnu);
	L1.Resize(constnu, 1);
	V1.Resize(constnu, 1);
	Pp1.Resize(constnu, constnu);

	Pxx.Resize(weizhishu, weizhishu);
	C1.Resize(all_const, weizhishu);
	W1.Resize(weizhishu, 1);
	Ncc1.Resize(all_const, all_const);

	for (int i = 0; i < all_const; i++)
	{
		for (int j = 0; j < presatinfo[1].m_Size + 1; j++)
		{
			if (prn_commsat[i] == presatinfo[1].m_Prn[j])
			{
				ele.push_back(presatinfo[1].Ele[j]);
				break;
			}
		}
	}
	for (int i = 0; i < all_const; i++)
	{
		if (i == 0)mk_zhy = 0;
		else mk_zhy += fre_comm[i - 1].size();
		for (int k = 0; k < fre_comm[i].size(); k++)
		{
			P1(mk_zhy + 1 + k, mk_zhy + 1 + k) = SQR((0.1) / (sin(ele[i])));// 2 * (SQR(0.003) / SQR(sin(ele[i])));SQR(gf_sigemaIon[prn_commsat[k]])
		}
	}
	int shuliang = comm_num;

	typedef struct xy
	{
		int x;
		int y;
	};

	FILE* fp = NULL;
	ofstream ofsitu, ofx, ofx1, ofx2, ofx3, ofxion, ofxslip, ofxorigincount, ofsigmaxpre, ofrecord, ofXX, oftime;
	string base, Zcountbase, ofsitub, ofxb, ofx1b, ofx2b, ofx3b, ofxionb, ofxslipb, ofxorigincountb, ofsigmaxpreb, ofrecordb, ofXXb, oftimeb, X11b;
	base = "C:\\Users\\15149\\Desktop\\shiy\\"; Zcountbase = to_string(popt.cccount);
	base += Zcountbase;
	ofsitub = "\\ofminsitu.out";
	ofxb = "\\sigmafinal.out";
	ofx2b = "\\sigmafinal1.out";
	ofx3b = "\\sigmafinalK.out";
	ofxionb = "\\ion.out";
	ofxslipb = "\\cslip.out";
	ofxorigincountb = "\\ofxorigincount.out";
	ofsigmaxpreb = "\\ofsigmapre.out";
	ofrecordb = "\\ofrecord.out";
	ofXXb = "\\ofXX.out";
	oftimeb = "\\oftime.out";
	X11b = "\\X11.out";

	if (0)//细节输出
	{
		ofxslip.open(base + ofxslipb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\cslip.out"
		ofxion.open(base + ofxionb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ion.out"
		ofsitu.open(base + ofsitub, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ofminsitu.out"
		ofrecord.open(base + ofrecordb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ofrecord.out"
		ofsigmaxpre.open(base + ofsigmaxpreb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ofsigmapre.out"
		oftime.open(base + oftimeb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\oftime.out"
		ofxorigincount.open(base + ofxorigincountb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ofxorigincount.out"
		ofXX.open(base + ofXXb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\ofXX.out"
		ofx.open(base + ofxb, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\sigmafinal.out"
		ofx3.open(base + ofx3b, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\sigmafinalK.out"
		ofx2.open(base + ofx2b, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\sigmafinal1.out"
		ofx1b = "\\repair.txt";
		ofx1.open(base + ofx1b, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\repair.txt"
	}
	else
	{
		string baserepir1, baserepir2, baserepir3, baserepir4, baserepir5, baserepir6, base1;
		base1 = "C:\\Users\\15149\\Desktop\\shiy\\";
		baserepir1 = "\\"; baserepir2 = popt.StaName; baserepir3 = "_repair"; baserepir4 = to_string(popt.schemeset); baserepir5 = ".txt";
		baserepir6 = to_string(popt.cccount);
		ofx1b += baserepir1; ofx1b += baserepir2; ofx1b += baserepir3; ofx1b += baserepir6; ofx1b += baserepir5;
		base1 += baserepir4;
		ofx1.open(base1 + ofx1b, ios::app);//"C:\\Users\\15149\\Desktop\\shiy\\repair.txt"
	}

	vector<double> insigma;
	multiset<double> insigma_mset;
	vector<int> count_sit_min;
	vector<int> storesitu;

	int count_num_a = 0;
	int count_insigma = 0;

	int count_dowhile = 0;
	int controlB = 0;
	double listminvalue = 0.0;
	vector<double> dioncoe2;
	vector<vector<double>> dioncoe1;

	for (int i = 0; i < GetRlenth.size(); i++)
	{
		for (int j = 0; j < GetRlenth[i].size(); j++)
		{
			if (j > 0)
			{
				dioncoe2.push_back(((GetRfrequ[i][0] * GetRfrequ[i][0]) / (GetRfrequ[i][j] * GetRfrequ[i][j])));
			}
			else
			{
				dioncoe2.push_back(1.0);
			}
		}
		if (dioncoe2.size() > 0)
		{
			dioncoe1.push_back(dioncoe2);
			dioncoe2.clear();
		}
	}


	//添加周跳
	int clip_map = 0;			//周跳位置范围
	int clip_range = 100;		//周跳范围
	vector<int> circ;			//周跳个数
	int porm = 2;				//正负
	int cclip_value = 0;		//周跳数值
	int crow = 0, ccol = 0;

	double sigmaofd[4] = { 14365.4679, 43143.07029, 68248.05289, 24485.29627 }; //{1/7200,1/7200,1/7200,1/7200};{ 14365.4679, 43143.07029, 68248.05289, 24485.29627 }; {45.80855642, 61.19817971, 2.872832186, 6.608876568};
	double aveofd[4] = { 0, 0, 0, 0 };

	vector<int> scclip_sit;//添加周跳的序列位置
	vector<int> scclip_value;//添加周跳的值
	vector<vector<double>> recordsucc;
	vector<int> origisit; int countorigin = 0.0;
	circ.clear();
	vector<double> standV; multiset<double> standV_mset;
	vector<int> prnnumb; double dxdydz = 0.0;

	ofxslip << flag_epoch_count - 1 << '\t' << constnu << '\t';

	if (0)//先做一遍求得原始的周跳情况if (flag_epoch_count == epoch01) if (1)
	{
		ofx2.open("C:\\Users\\15149\\Desktop\\shiy\\sigmafinal1.out", ios::app);
		count_dowhile = 0;
		controlB = 0;
		storesitu.clear(); e_2.clear(); rs_2.clear();

		B1.Clear(); Nbb1.Clear(); X1.Clear(); L1.Clear(); V1.Clear(); VZ.Clear(); Pp1.Clear(); C1.Clear(); Ncc1.Clear(); v__1.Clear(); p__1.Clear(); Sigema0__1.Clear(); Pxx.Clear(); W1.Clear();
		L1.Resize(constnu, 1); V1.Resize(constnu, 1); Pp1.Resize(constnu, constnu); v__1.Resize(comm_num, 1); p__1.Resize(comm_num, comm_num); Sigema0__1.Resize(1, 1); Ncc1.Resize(all_const, all_const);

		do
		{
			listminvalue = 100.0;
			do
			{
				weizhishu = 4 + 1 + all_const + count_dowhile - controlB;//  
				B1.Resize(constnu, weizhishu);
				Nbb1.Resize(weizhishu, weizhishu);
				X1.Resize(weizhishu, 1);
				C1.Resize(weizhishu, weizhishu);
				Pxx.Resize(weizhishu, weizhishu);
				W1.Resize(weizhishu, 1);
				VZ.Resize(weizhishu, 1);

				if (controlB == 0)
				{
					storesitu.push_back(5 + all_const + count_dowhile);
				}

				for (int shuchu = 0; shuchu < shuliang; shuchu++)
				{
					for (int i = 0; i < all_const; i++)
					{
						for (int j = 0; j < presatinfo[1].m_Size + 1; j++) {
							if (prn_commsat[i] == presatinfo[1].m_Prn[j])
							{
								satpos2[0] = presatinfo[1].m_Xs[3 * j + 0];
								satpos2[1] = presatinfo[1].m_Xs[3 * j + 1];
								satpos2[2] = presatinfo[1].m_Xs[3 * j + 2];

								for (int k = 0; k < 3; ++k)
								{
									rs2 += (RKpos_cur[k] - satpos2[k]) * (RKpos_cur[k] - satpos2[k]);
									xyz2[k] = RKpos_cur[k] - satpos2[k];																// (Xr-Xs)^2 +(Yr-Ys)^2+(Zr-Zs)^2
								}
								rs2 = sqrt(rs2);
								e2.e_x = xyz2[0] / rs2;
								e2.e_y = xyz2[1] / rs2;
								e2.e_z = xyz2[2] / rs2;
								e_2.push_back(e2);
								rs_2.push_back(rs2);
								rs2 = 0.0;
							}
						}
					}
					for (int m = 0; m < all_const; m++)
					{
						if (m == 0)mk_zhy = 0;
						else { mk_zhy += fre_comm[m - 1].size(); }
						for (int i = 0; i < fre_comm[m].size(); i++)
						{
							L1(mk_zhy + i + 1, 1) = (L_Cur[m][i] - rs_2[m]) + (rs_1[m] - L_Pre[m][i]);// - gf_ion[prn_commsat[m]] * dioncoe[i];//#t+1-t
							if (fabs(L1(mk_zhy + i + 1, 1)) > 20)
							{
								double a = L1(mk_zhy + i + 1, 1);
								a = fabs(a);
							}
							B1(mk_zhy + i + 1, 1) = (e_2[m].e_x);
							B1(mk_zhy + i + 1, 2) = (e_2[m].e_y);
							B1(mk_zhy + i + 1, 3) = (e_2[m].e_z);
							B1(mk_zhy + i + 1, 4) = 1;						//接收机钟差
							B1(mk_zhy + i + 1, m + 5) = dioncoe1[m][i];
							char buff = GetSystem_GREC2C3(prn_commsat[m]);
							int a = _SystemIndex[buff];
							//B1(mk_zhy + i + 1, 3 + _SystemIndex[GetSystem_GREC2C3(prn_commsat[m])]) = 1.0;//接收机钟差 一个系统就是多一个； 
						}
					}
					for (int i = 0; i < (4 + all_const); i++)
					{
						if (i <= 3)
						{
							if (i == 0)
							{
								W1(i + 1, 1) = aveofd[i];
							}
							if (i == 1)
							{
								W1(i + 1, 1) = aveofd[i];
							}
							if (i == 2)
							{
								W1(i + 1, 1) = aveofd[i];
							}
							if (i == 3)
							{
								W1(i + 1, 1) = aveofd[i];
							}
							Pxx(i + 1, i + 1) = sigmaofd[i];
						}
						if (i > 3)
						{
							W1(i + 1, 1) = ion_gimsub[i - 4];
							Pxx(i + 1, i + 1) = dionsigma;
						}
					}
					for (int i = 0; i < weizhishu; i++)
					{
						C1(i + 1, i + 1) = 1.0;
					}
					if (controlB == 0)
					{
						for (int i = 0; i < constnu; i++)  //# add canshu dia
						{
							if (i == shuchu)
							{
								B1(i + 1, storesitu[count_dowhile]) = 1;
								break;
							}
						}
					}

					for (int i = 0; i < count_sit_min.size(); i++)
					{
						B1(count_sit_min[i] + 1, storesitu[i]) = 1;
					}

					//if (LeastSquare(B1, L1, P1, V1, Nbb1, X1) == -1);    // LSQ  Nbb=(btpb)^-1
					if (LeastSquareconstraint(B1, L1, P1, C1, W1, V1, VZ, Nbb1, Ncc1, X1, Pxx) == -1);

					//B1.Output_to_File("E:\\B1.out", 1);
					//L1.Output_to_File("E:\\L1.out", 1);

					double sigema0__1 = 0;

					Sigema0 = (V1.Trans() * P1 * V1) + (VZ.Trans() * Pxx * VZ);
					sigema0 = Sigema0.GetValue(1, 1);
					sigema0 = sqrt(sigema0 / (constnu - 4));

					count_num_a = 0;
					if (count_sit_min.size() == 0)
					{
						insigma.push_back(sigema0);
						insigma_mset.insert(sigema0);
					}
					else
					{
						for (int i = 0; i < count_sit_min.size(); i++)
						{
							if (shuchu == count_sit_min[i])
							{
								count_num_a++;
							}
						}
						if (count_num_a != 0)
						{
							insigma.push_back(0.0);
						}
						else
						{
							insigma.push_back(sigema0);
							insigma_mset.insert(sigema0);
						}
					}

					if (controlB == 0)
					{
						for (int i = 0; i < constnu; i++)  // # add canshu dia
						{
							if (i == shuchu)
							{
								B1(i + 1, storesitu[count_dowhile]) = 0;
								break;
							}

						}
					}
				}

				count_insigma = 0;
				if (controlB == 0)
				{
					for (int i = 0; i < insigma.size(); i++)
					{
						double Hypoth_test = (insigma[i] - popt.Zavemin) / (popt.Zsigma);
						if (Hypoth_test > 2.58 && insigma[i] < (3 * listminvalue))
						{
							count_insigma++;
						}
					}
				}
				if (controlB == 1)
				{
					double sigmafinal = 0;
					double sizein = insigma_mset.size();
					for (multiset<double>::iterator it = insigma_mset.begin(); it != insigma_mset.end(); it++)
					{
						sigmafinal += (*it);
					}
					sigmafinal = sigmafinal / sizein;
					ofx2 << sigmafinal << endl;
				}
				if (count_insigma != 0)
				{
					for (int i = 0; i < insigma.size(); i++)
					{
						if (insigma[i] == *(insigma_mset.begin()))
						{
							count_sit_min.push_back(i);
							listminvalue = *(insigma_mset.begin());
							break;
						}
					}
					count_dowhile++;
				}
				else
				{
					controlB++;
					if (controlB == 2)
					{
						X1.Output_to_File("D:\\data10.0\\out\\X11.out", 1);
						ofsigmaxpre << X1(1, 1) << '\t' << X1(2, 1) << '\t' << X1(3, 1) << '\t' << X1(4, 1) << endl;
					}
				}
				insigma_mset.clear();
				insigma.clear();

			} while (!count_insigma == 0);
		} while (controlB != 2);
		for (int j = 0; j < count_sit_min.size(); j++)
		{
			double temp = X1(5 + all_const + j, 1) / (GetRlenth[(count_sit_min[j] / frenumcontrol1)][(count_sit_min[j] % frenumcontrol1)]);
			if (abs(temp) >= 0.5)
			{
				origisit.push_back(count_sit_min[j]);
				countorigin += 1;
			}
		}
	}
	ofxorigincount << countorigin << endl;
	int countsize = 0;

	for (int i = 0; i < (comm_num - prnnumb.size()); i++)
	{
		double a = ((i + 1.0) / comm_num);
		if (a >= percentl && a < percentr)
		{
			circ.push_back(0);
		}
	}

	clock_t start1, end1;
	vector<double> time1;
	double time11 = 0.0;
	if (1)
	{
		for (int p = 0; p < circ.size(); p++)
		{
			vector<vector<double>> L_Curtemp;
			L_Curtemp = L_Cur;
			if (circ[p] > countorigin)//if (flag_epoch_count == epoch01)
			{
				scclip_sit.clear(); scclip_value.clear();
				for (int i = 0; i < countorigin; i++)
				{
					clip_map = origisit[i];
					crow = clip_map / frenumcontrol1;
					ccol = clip_map % frenumcontrol1;
					if ((rand() % porm) == 0)
					{
						cclip_value = rand() % clip_range + 1;
						L_Curtemp[crow][ccol] -= cclip_value * GetRlenth[crow][ccol];
						scclip_sit.push_back(clip_map);
						scclip_value.push_back(-cclip_value);
						ofxslip << scclip_sit[i] + 1 << '!' << scclip_value[i] << '\t';
					}
					else
					{
						cclip_value = rand() % clip_range + 1;
						L_Curtemp[crow][ccol] += cclip_value * GetRlenth[crow][ccol];
						scclip_sit.push_back(clip_map);
						scclip_value.push_back(cclip_value);
						ofxslip << scclip_sit[i] + 1 << '!' << scclip_value[i] << '\t';
					}
				}
				int i1 = 0, counta = 0;
				do {
					counta = 0; int temp = 0.0;
					do
					{
						temp = 0.0;
						clip_map = rand() % comm_num;
						for (int i = 0; i < prnnumb.size(); i++)
						{
							if (clip_map == prnnumb[i])
							{
								temp++;
							}
						}
					} while (temp);
					crow = clip_map / frenumcontrol1;
					ccol = clip_map % frenumcontrol1;
					for (int j = 0; j < scclip_sit.size(); j++)
					{
						if (clip_map == scclip_sit[j])
						{
							counta++;
							break;
						}
					}
					if (counta == 0)
					{
						if ((rand() % porm) == 0)
						{
							cclip_value = rand() % clip_range + 1;
							L_Curtemp[crow][ccol] -= cclip_value * GetRlenth[crow][ccol];
							scclip_sit.push_back(clip_map);
							scclip_value.push_back(-cclip_value);
							ofxslip << scclip_sit[countorigin + i1] + 1 << '!' << scclip_value[countorigin + i1] << '\t';
						}
						else
						{
							cclip_value = rand() % clip_range + 1;
							L_Curtemp[crow][ccol] += cclip_value * GetRlenth[crow][ccol];
							scclip_sit.push_back(clip_map);
							scclip_value.push_back(cclip_value);
							ofxslip << scclip_sit[countorigin + i1] + 1 << '!' << scclip_value[countorigin + i1] << '\t';
						}
						i1++;
					}
				} while (countorigin + i1 < circ[p]);
			}
			else if (circ[p] <= countorigin)
			{
				countsize++;
				if (countsize == circ.size())
				{
					break;
				}
			}
			if (1)//抗差最小二乘优化权阵if (flag_epoch_count == epoch01) if (1)
			{
				double recorddx0 = 0.0, recorddx1 = 0.0, re = 0.0; dxdydz = 0.0; int countset = 0;
				start1 = clock();
				do
				{
					countset++;
					//	if (countset == 3)
					//	{
						//	break;
					//	}
					recorddx0 = recorddx1;
					standV.clear();
					standV_mset.clear();
					e_2.clear();
					rs_2.clear();
					B1.Clear(); Nbb1.Clear(); X1.Clear(); L1.Clear(); V1.Clear(); VZ.Clear(); C1.Clear(); Pxx.Clear(); W1.Clear(); Qvv.Clear();
					L1.Resize(constnu, 1); V1.Resize(constnu, 1); Qvv.Resize(constnu, constnu);

					if (OTDorITD == 0)
					{
						weizhishu = 4;
					}
					else
					{
						weizhishu = 4 + all_const;
					}
					B1.Resize(constnu, weizhishu);
					Nbb1.Resize(weizhishu, weizhishu);
					X1.Resize(weizhishu, 1);

					for (int i = 0; i < all_const; i++)
					{
						for (int j = 0; j < presatinfo[1].m_Size + 1; j++) {
							if (prn_commsat[i] == presatinfo[1].m_Prn[j])
							{
								satpos2[0] = presatinfo[1].m_Xs[3 * j + 0];
								satpos2[1] = presatinfo[1].m_Xs[3 * j + 1];
								satpos2[2] = presatinfo[1].m_Xs[3 * j + 2];

								for (int k = 0; k < 3; ++k)
								{
									rs2 += (RKpos_cur[k] - satpos2[k]) * (RKpos_cur[k] - satpos2[k]);
									xyz2[k] = RKpos_cur[k] - satpos2[k];																// (Xr-Xs)^2 +(Yr-Ys)^2+(Zr-Zs)^2
								}
								rs2 = sqrt(rs2);
								e2.e_x = xyz2[0] / rs2;
								e2.e_y = xyz2[1] / rs2;
								e2.e_z = xyz2[2] / rs2;
								e_2.push_back(e2);
								rs_2.push_back(rs2);
								rs2 = 0.0;
							}
						}
					}
					for (int m = 0; m < all_const; m++)
					{
						if (m == 0)mk_zhy = 0;
						else { mk_zhy += fre_comm[m - 1].size(); }
						for (int i = 0; i < fre_comm[m].size(); i++)
						{
							L1(mk_zhy + i + 1, 1) = (L_Curtemp[m][i] - rs_2[m]) + (rs_1[m] - L_Pre[m][i]) + dxdydz;// - gf_ion[prn_commsat[m]] * dioncoe[i];//#t+1-t
							if (fabs(L1(mk_zhy + i + 1, 1)) > 20)
							{
								double a = L1(mk_zhy + i + 1, 1);
								a = fabs(a);
							}
							B1(mk_zhy + i + 1, 1) = (e_2[m].e_x);
							B1(mk_zhy + i + 1, 2) = (e_2[m].e_y);
							B1(mk_zhy + i + 1, 3) = (e_2[m].e_z);
							B1(mk_zhy + i + 1, 4) = 1;						//接收机钟差
							if (OTDorITD != 0)
							{
								B1(mk_zhy + i + 1, m + 5) = dioncoe1[m][i];
							}
						}
					}

					if (LeastSquareQcc(B1, L1, P1, V1, Nbb1, Qvv, X1) == -1);

					Sigema0 = (V1.Trans() * P1 * V1);
					sigema0 = Sigema0.GetValue(1, 1);
					sigema0 = sqrt(sigema0 / (constnu - 4));

					//标准化残差并等价权比较
					for (int i = 0; i < constnu; i++)
					{
						standV.push_back(abs(V1(i + 1, 1) / (sigema0 * sqrt(Qvv(i + 1, i + 1)))));
						standV_mset.insert(abs(V1(i + 1, 1) / (sigema0 * sqrt(Qvv(i + 1, i + 1)))));
					}
					double k0 = 0.8;
					double k1 = 1.3;
					for (int i = 0; i < constnu; i++)
					{
						if (standV[i] > k1)
						{
							int tempset = 0;
							P1(i + 1, i + 1) = 1e-6;
							for (int j = 0; j < prnnumb.size(); j++)
							{
								if (i == prnnumb[j])
								{
									tempset++;
								}
							}
							if (tempset == 0)
							{
								prnnumb.push_back(i);
							}
						}
						else if (standV[i] > k0 && standV[i] <= k1)
						{
							P1(i + 1, i + 1) *= (k0 / standV[i]) * (((k1 - standV[i]) / (k1 - k0)) * (((k1 - standV[i]) / (k1 - k0))));
						}
					}
					dxdydz = sqrt((X1(1, 1) * X1(1, 1)) + (X1(2, 1) * X1(2, 1)) + (X1(3, 1) * X1(3, 1)));
					recorddx1 = X1(1, 1);
					re = abs((recorddx0 - recorddx1));
					if ((re < 1e-5))
					{
						ofx3 << sigema0 << endl;
						P1.Output_to_File("C:\\Users\\15149\\Desktop\\shiy\\P1.out", 1);
					}
				} while (!(re < 1e-5));
			}
			ofxorigincount << prnnumb.size() << '\t';
			if (1)//if (flag_epoch_count == epoch01) if (1) countsize == 0
			{
				//ofsitu.open("C:\\Users\\15149\\Desktop\\shiy\\minsitu.out", ios::out);
				storesitu.clear();
				count_sit_min.clear(); e_2.clear(); rs_2.clear();
				B1.Clear(); Nbb1.Clear(); X1.Clear(); L1.Clear(); V1.Clear(); VZ.Clear(); Pp1.Clear(); C1.Clear(); Ncc1.Clear(); v__1.Clear(); p__1.Clear(); Sigema0__1.Clear(); Pxx.Clear(); W1.Clear();
				L1.Resize(constnu, 1); V1.Resize(constnu, 1); Pp1.Resize(constnu, constnu); v__1.Resize(comm_num, 1); p__1.Resize(comm_num, comm_num); Sigema0__1.Resize(1, 1); Ncc1.Resize(all_const, all_const);
				vector<double> repairvalue;

				for (int t = 0; t < prnnumb.size(); t++)
				{
					if (OTDorITD == 0)
					{
						weizhishu = 4 + 1;
					}
					else
					{
						weizhishu = 4 + 1 + all_const;
					}
					B1.Resize(constnu, weizhishu);
					Nbb1.Resize(weizhishu, weizhishu);
					X1.Resize(weizhishu, 1);

					//fp = fopen("C:\\Users\\15149\\Desktop\\shiy\\xigema.out", "a+");

					for (int i = 0; i < all_const; i++)
					{
						for (int j = 0; j < presatinfo[1].m_Size + 1; j++) {

							if (prn_commsat[i] == presatinfo[1].m_Prn[j])
							{
								satpos2[0] = presatinfo[1].m_Xs[3 * j + 0];
								satpos2[1] = presatinfo[1].m_Xs[3 * j + 1];
								satpos2[2] = presatinfo[1].m_Xs[3 * j + 2];

								for (int k = 0; k < 3; ++k)
								{
									rs2 += (RKpos_cur[k] - satpos2[k]) * (RKpos_cur[k] - satpos2[k]);
									xyz2[k] = RKpos_cur[k] - satpos2[k];																// (Xr-Xs)^2 +(Yr-Ys)^2+(Zr-Zs)^2
								}
								rs2 = sqrt(rs2);
								e2.e_x = xyz2[0] / rs2;
								e2.e_y = xyz2[1] / rs2;
								e2.e_z = xyz2[2] / rs2;
								e_2.push_back(e2);
								rs_2.push_back(rs2);
								rs2 = 0.0;
							}
						}
					}
					for (int m = 0; m < all_const; m++)
					{
						if (m == 0)mk_zhy = 0;
						else { mk_zhy += fre_comm[m - 1].size(); }
						for (int i = 0; i < fre_comm[m].size(); i++)
						{
							L1(mk_zhy + i + 1, 1) = (L_Curtemp[m][i] - rs_2[m]) + (rs_1[m] - L_Pre[m][i]);
							if (fabs(L1(mk_zhy + i + 1, 1)) > 20)
							{
								double a = L1(mk_zhy + i + 1, 1);
								a = fabs(a);
							}
							B1(mk_zhy + i + 1, 1) = (e_2[m].e_x);
							B1(mk_zhy + i + 1, 2) = (e_2[m].e_y);
							B1(mk_zhy + i + 1, 3) = (e_2[m].e_z);
							B1(mk_zhy + i + 1, 4) = 1.0;						
							if (OTDorITD != 0)
							{
								B1(mk_zhy + i + 1, m + 5) = dioncoe1[m][i];
							}
						}
					}

					P1(prnnumb[t] + 1, prnnumb[t] + 1) = SQR((0.1) / (sin(ele[(prnnumb[t] / frenumcontrol1)])));
					B1(prnnumb[t] + 1, weizhishu) = 1.0;

					if (LeastSquare(B1, L1, P1, V1, Nbb1, X1) == -1);    // LSQ  Nbb=(btpb)^-1

					Sigema0 = (V1.Trans() * P1 * V1);
					sigema0 = Sigema0.GetValue(1, 1);
					sigema0 = sqrt(sigema0 / (constnu - 4));


					double tempx = 0.0;
					tempx = X1(weizhishu, 1);

					L_Curtemp[(prnnumb[t] / frenumcontrol1)][(prnnumb[t] % frenumcontrol1)] -= tempx;
					repairvalue.push_back(tempx);
					ofXX << X1(weizhishu, 1) << '\t';
					ofsitu << prnnumb[t] << '\t';

					X1.Output_to_File("D:\\data10.0\\out\\X1.out", 1);
				}
				end1 = clock();
				time1.push_back(((double(end1 - start1)) / (CLOCKS_PER_SEC)));

				//ofx1 << flag_epoch_count - 1<< '\t';
				vector<double> recordsucc1;
				double count_compare = 0;
				for (int i = 0; i < scclip_sit.size(); i++)
				{
					for (int j = 0; j < prnnumb.size(); j++)
					{
						if (scclip_sit[i] == prnnumb[j])
						{
							double temp = (scclip_value[i] - (repairvalue[j] / GetRlenth[(scclip_sit[i] / frenumcontrol1)][(scclip_sit[i] % frenumcontrol1)]));//(scclip_value[i] - (X1(5 + all_const + j, 1) / slipwave[(scclip_sit[i] % 3)]));(scclip_value[i] * slipwave[(scclip_sit[i] % 3)]) - (X1(5 + all_const + j, 1))
							//double temp = (scclip_value[i] - (X1(5 + all_const + j, 1) / slipwave[(scclip_sit[i] % frenumcontrol1)]));
							if (abs(temp) < 0.6)
							{
								//ofx1 << 1 << '\t';
								recordsucc1.push_back(1);
							}
							else
							{
								//ofx1 << 0 << '\t';
								recordsucc1.push_back(0);
							}
						}
						else
						{
							count_compare++;
						}
					}
					if (count_compare == prnnumb.size())
					{
						//ofx1 << 0 << '\t';
						recordsucc1.push_back(0);
						count_compare = 0;
					}
					else
					{
						count_compare = 0;
					}
				}
				recordsucc.push_back(recordsucc1);
				for (int i = 0; i < recordsucc1.size(); i++)
				{
					ofrecord << recordsucc1[i] << ' ';
				}
				ofrecord << '\t';

				if (p == circ.size() - 1)
				{
					vector<double> crs;
					double crs1 = 0.0;
					double crs2 = 0.0;
					for (int j = 0; j < recordsucc.size(); j++)
					{
						if (recordsucc[j].size() != 0)
						{
							for (int i = 0; i < recordsucc[j].size(); i++)
							{
								crs2 += recordsucc[j][i];

								if (i == recordsucc[j].size() - 1)
								{
									crs.push_back(crs2);
									crs[j] = crs[j] / (recordsucc[j].size());
									crs2 = 0.0;
								}
							}
							crs1 += crs[j];
						}
					}
					double ktek = crs1 / recordsucc.size();
					ofx1 << ktek << '\t';
					ofx1 << endl;
				}
				prnnumb.clear(); e_2.clear(); rs_2.clear();
			}
			countsize = 0;
		}
		for (int i = 0; i < time1.size(); i++)
		{
			time11 += time1[i];
		}
		oftime << time11 << endl;
		ofrecord << endl;
		ofsitu << endl;
		ofXX << endl;
		ofsitu.close();

		return 0;
	}
}

/// Least square filting option
int LeastSquarePPP::Lsq_Option_Init(ppp_option_t& popt)
{
	m_dt = popt.dt;     // Observation interval

	//% 1.  Option Setting from config
	// Parameter Estimation 
	// 测试
	est_xyz = popt.estpara.xyz;
	//est_xyz = 1;
	est_rclk = popt.estpara.rclk;
	est_trop = popt.estpara.trop;
	est_ion = popt.estpara.ion;

	// Process Option



	_Mode = popt.mode_ppp;     // Static/Kinematic

	//测试
	//_Mode = 0;
	_MathModel = popt.math_model;                   // Math function model
	//_MathModel = MODE_UNCOMBINE;   
	//_MathModel = 2;

	_GIM_P1P2 = popt.GIM_P1P2;
	_IS_DCB = popt.IS_DCB;
	_IS_ISB = popt.IS_ISB;
	_ISBMODE = popt.ISB_MODE;

	_IonModel = ION_MODEL_EST_WHITE;               // Ionosphere model 1: White nois
	_TropModel = 0;                                // Troposphere model 0: Est zwd
	_RclkModel = RCLK_MODEL_WHITE_NOISE;           // Receiver clock model: 0: white noise
	_Pmode = popt.process_mode;                    // Process mode: 0: spp 1: ppp
	_FreqNum = popt.freqn;                         // Frequency number 1
	_System = popt.System;                         // System need to process G,R,E,C
	//popt.mode_filter=2;                              //滤波模式###################
	_FilterMode = popt.mode_filter;                // Filter mode 0: Forward; 1: Backward
	_TimeMode = popt.process_time;                 // Post process; Realtime process

	_traceMode = popt.trace;                       // Log in file
	//% Option Setting from config end

	//% 2. Parameter number setting
	if (_MathModel == MODE_UNCOMBINE)       // Uncombined PPP
	{
		_AmbNum = 1;                       // ambiguity number: m
		_IonNum = 1;                       // ionospere number: m
		_ObsEqNum = 2;                     // observation number: 2m (P3, L3)

		_ParaN_sat = _AmbNum + _IonNum;    // parameter number associate with satellite: 3m or 2m
	}
	if (_MathModel == MODE_UNCOMBINE_ION_CONT)       // 电离层约束
	{
		_AmbNum = 1;                       // ambiguity number: m
		_IonNum = 1;                       // ionospere number: m
		_ObsEqNum = 3;                     // observation number: 2m (P3, L3)

		_ParaN_sat = _AmbNum + _IonNum;    // parameter number associate with satellite: 3m or 2m
	}
	if (_MathModel == MODE_UOFC_COMBINE)       // Uncombined PPP
	{
		_AmbNum = 1;                       // ambiguity number: m
		_IonNum = 0;                       // ionospere number: m
		_ObsEqNum = 2;                     // observation number: 2m (P3, L3)

		_ParaN_sat = _AmbNum + _IonNum;    // parameter number associate with satellite: 3m or 2m
	}
	if (_MathModel == MODE_IF_COMBINE)       // Uncombined PPP
	{
		_AmbNum = 1;                       // ambiguity number: m
		_IonNum = 0;                       // ionospere number: m
		_ObsEqNum = 2;                     // observation number: 2m (P3, L3)

		_ParaN_sat = _AmbNum + _IonNum;    // parameter number associate with satellite: 3m or 2m
	}
	if (_MathModel == 4)       // Uncombined PPP
	{
		_AmbNum = 1;                       // ambiguity number: m
		_IonNum = 1;                       // ionospere number: m
		_ObsEqNum = 6;                     // observation number: 2m (P3, L3)

		_ParaN_sat = _AmbNum * 3 + _IonNum;    // parameter number associate with satellite: 3m or 2m
	}
	_PsuedoRclkNum = _RclkModel;           // RCLK_MODEL_WHITE_NOISE: 0 

	_TropNum = est_trop;                   // wet zenith delay

	// Parameter number in filter
	if (_Mode == MODE_STATIC)	XYZ_Num = 3 * est_xyz;
	if (_Mode == MODE_KINEMATIC)	XYZ_Num = 0;
	XYZ_Rclk_Num = XYZ_Num + _PsuedoRclkNum;
	XYZ_Rclk_Trop_Num = XYZ_Num + _PsuedoRclkNum + _TropNum;
	//% Parameter number setting end

	//% 3. Stochastic model (Priori variance)
	if (popt.sigmaP > 0)
		m_SigemaP = popt.sigmaP;
	else
		m_SigemaP = SIGEMA_P;

	if (popt.sigmaL > 0)
		m_SigemaL = popt.sigmaL;
	else
		m_SigemaL = SIGEMA_L;

	_TropProcessNoise = popt.trop_noise;         // trop delay random-walk noise
	_TropGNnoise = popt.trop_gradientNnoise;     // trop gradient N random-walk noise
	_TropGEnoise = popt.trop_gradientEnoise;     // trop gradient E random-walk noise

	_IonProcessNoise = popt.ion_noise;           // ion delay random-walk noise


	if (_RclkModel == RCLK_MODEL_RANDWALK)
		_RclkNoise = popt.rclk_noise;            // receiver clock random-walk noise
	//% Stochastic model (Priori variance) end

	// Receiver antenna offset
	Antdiff[0] = popt.antenna_offset[0];
	Antdiff[1] = popt.antenna_offset[1];
	Antdiff[2] = popt.antenna_offset[2];

	if (_Mode == MODE_KINEMATIC)
	{
		if (Norm(popt.StaPos_SPP, 3)) { // Obtain from SPP
			CoarseXYZ[0] = popt.StaPos_SPP[0];
			CoarseXYZ[1] = popt.StaPos_SPP[1];
			CoarseXYZ[2] = popt.StaPos_SPP[2];
		}
		else if (Norm(popt.ProXYZ, 3)) { // Obtain from head file
			CoarseXYZ[0] = popt.ProXYZ[0];
			CoarseXYZ[1] = popt.ProXYZ[1];
			CoarseXYZ[2] = popt.ProXYZ[2];
		}
	}
	else if (_Mode == MODE_STATIC)  // Attention: Get approximate position only once for static mode
	{
		if (flag_epoch_count == 0)
		{
			if (!Norm(popt.ProXYZ, 3)) { // Obtain from head file
				CoarseXYZ[0] = popt.ProXYZ[0];
				CoarseXYZ[1] = popt.ProXYZ[1];
				CoarseXYZ[2] = popt.ProXYZ[2];
			}
			else if (!Norm(popt.StaPos_SPP, 3)) { // Obtain from SPP
				CoarseXYZ[0] = popt.StaPos_SPP[0];
				CoarseXYZ[1] = popt.StaPos_SPP[1];
				CoarseXYZ[2] = popt.StaPos_SPP[2];
			}
			else if (Norm(popt.StaPos_CRD, 3)) { // Obtain from SNX/CRD file
				CoarseXYZ[0] = popt.StaPos_CRD[0];
				CoarseXYZ[1] = popt.StaPos_CRD[1];
				CoarseXYZ[2] = popt.StaPos_CRD[2];
			}
		}
	}
	m_Sta[0] = CoarseXYZ[0];
	m_Sta[1] = CoarseXYZ[1];
	m_Sta[2] = CoarseXYZ[2];

	return 1;
}
void LeastSquarePPP::cal_IFB_num(ppp_option_t& popt, PPPObsData& pppdata)
{
	int b2 = 0, b1c = 0, b2a = 0, e5b = 0, e5 = 0, e6 = 0, l5 = 0; int index = 0;
	for (int i = 0; i < pppdata.m_SatCount; ++i)
	{
		for (int j = 0; j < pppdata._num[i]; j++)
		{
			if (pppdata.Prn_Otype[i][j + 1] == "B2")	 b2 = 1;
			if (pppdata.Prn_Otype[i][j + 1] == "B1C")	 b1c = 1;
			if (pppdata.Prn_Otype[i][j + 1] == "B2a")	 b2a = 1;
			if (pppdata.Prn_Otype[i][j + 1] == "E5b")	 e5b = 1;
			if (pppdata.Prn_Otype[i][j + 1] == "E5(ab)") e5 = 1;
			if (pppdata.Prn_Otype[i][j + 1] == "E6")	 e6 = 1;
			if (pppdata.Prn_Otype[i][j + 1] == "L5")	 l5 = 1;
		}
	}
	popt.IFB_num = b2 + b1c + b2a + e5b + e5 + e6 + l5;
	for (int i = 0; i < pppdata.m_SatCount; ++i)
	{
		for (int j = 0; j < pppdata._num[i]; j++)
		{
			index = get_fcbindex(pppdata.Prn_Otype[i][j + 1]);
			switch (index)
			{
			case 1: { _IFB_Index[index] = 1; break; }
			case 2: { _IFB_Index[index] = b2 + 1; break; }
			case 3: { _IFB_Index[index] = b1c + b2 + 1; break; }
			case 4: { _IFB_Index[index] = b2a + b1c + b2 + 1; break; }
			case 5: { _IFB_Index[index] = e5b + b2a + b1c + b2 + 1; break; }
			case 6: { _IFB_Index[index] = e5 + e5b + b2a + b1c + b2 + 1; break; }
			case 7: { _IFB_Index[index] = e6 + e5 + e5b + b2a + b1c + b2 + 1; break; }

			default:break;
			}
		}
	}
}
/// Set ISB index and get system number
int LeastSquarePPP::SetISBindex(PPPObsData& pppdata)
{
	char Sys_bds = 0;
	char sys = 0;
	int sysnum = 1;

	int sys_gps = 0, sys_glo = 0, sys_gal = 0, sys_bds = 0;
	int sys_bds2 = 0, sys_bds3 = 0;
	_GPS_ISB = 0; _GLO_ISB = 0; _GAL_ISB = 0; _BDS_ISB = 0, _ISBNum = 0;
	_BDS2_ISB = 0; _BDS3_ISB = 0;
	for (int i = 0; i < pppdata.m_SatCount; ++i)
	{
		sys = GetSystem_GREC2C3(pppdata.m_Prn[i]);

		switch (sys)
		{
		case 'G': { sys_gps = 1; break; }
		case 'R': { sys_glo = 1; break; }
		case 'E': { sys_gal = 1; break; }
		case 'C': { sys_bds3 = 1; break; }
		case 'B': { sys_bds2 = 1; break; }
		default:break;
		}
	}

	// Only for Single-frequency PPP ------------------------------------------
	/*if ((_FreqNum == 1) && (flag_epoch_count > 1))
	{
		sys_gps = 0, sys_glo = 0, sys_gal = 0, sys_bds = 0;
		for (int i = 0; i < pppdata.m_SatCount; ++i)
		{
			sys = GetSystem(pppdata.m_Prn[i]);

			for (int j = 0; j < _Prn_previous.size(); ++j)
			{
				if (pppdata.m_Prn[i] == _Prn_previous[j])
				{
					switch (sys)
					{
					case 'G': {_GPS_ISB = 1; sys_gps = 1; break; }
					case 'R': {_GLO_ISB = 1; sys_glo = 1; break; }
					case 'E': {_GAL_ISB = 1; sys_gal = 1; break; }
					case 'C': {_BDS_ISB = 1; sys_bds = 1; break; }
					default:break;
					}
					break;
				}
			}
		}
	}*/
	// -------------------------------------------------------------------------

	// Set ISB index ---------------------------------------------------------------
	memset(_SystemIndex, 0, sizeof(char) * 255);
	for (int i = 0; i < pppdata.m_SatCount; ++i)
	{
		sys = GetSystem_GREC2C3(pppdata.m_Prn[i]);

		switch (sys)
		{
		case 'G': { _SystemIndex[sys] = sys_gps; break; }
		case 'R': { _SystemIndex[sys] = sys_gps + sys_glo; break; }
		case 'E': { _SystemIndex[sys] = sys_gps + sys_glo + sys_gal; break; }
		case 'C': { _SystemIndex[sys] = sys_gps + sys_glo + sys_gal + sys_bds3; break; }
		case 'B': { _SystemIndex[sys] = sys_gps + sys_glo + sys_gal + sys_bds3 + sys_bds2; break; }
		default:break;
		}
	}
	// Set ISB index end -----------------------------------------------------------

	_SystemN = sys_gps + sys_glo + sys_gal + sys_bds2 + sys_bds3;  // System number
	_ISBNum = _SystemN - 1;                            // ISB number


	// FOR  BDS_2_3  测试  实验  ISB  系统间偏差
	if (_IS_ISB) {
		//FOR  BDS_2_3   实验  ISB  系统间偏差
		for (int i = 0; i < pppdata.m_SatCount; ++i)
		{
			if (sys = GetSystem(pppdata.m_Prn[i]) == 'C')
			{
				switch (Sys_bds = GetSystem_BDS_2_3(pppdata.m_Prn[i]))
				{
				case 'C': { _SystemIndex[Sys_bds] = 1; break; }
				case 'B': { _SystemIndex[Sys_bds] = 2; break; }
				default:break;
				}
			}
		}
		_SystemN = 2;  // System number
		_ISBNum = _SystemN - 1;
	}
	return _SystemN;
}

/// Initialize least square filter parameter
int LeastSquarePPP::Lsq_ParaInit(ppp_option_t& popt, unsigned int SatN)
{
	unsigned int EqNum = 0, ParaNum = 0;

	// Number of const estimated parameters [X, Y, Z, Rclk, ISB, Trop]
	_ParaN_const = (3 * est_xyz) + est_rclk + _ISBNum + _TropNum;
	if (_MathModel == 4) {
		_ParaN_EstTotal = popt.CUR_mum + SatN + _ParaN_const; // Total estimate parameter number
	}
	else
	{
		_ParaN_EstTotal = _ParaN_sat * SatN + _ParaN_const;
	}
	ParaNum = _ParaN_EstTotal;              // Total estimated parameter number
	if (_MathModel == 4)
	{
		EqNum = popt.CUR_mum * 2 + _PsuedoObsN;
	}
	else
	{
		EqNum = _ObsEqNum * SatN + _PsuedoObsN;
	}
	// Equation number (real-observation + psuedo-observation)

	if (popt.ion_const)
	{
		Q.Resize(EqNum + SatN, EqNum + SatN);                 // Covariance MatrixT of observation
		P.Resize(EqNum + SatN, EqNum + SatN);                 // Weight MatrixT of observation
		B.Resize(EqNum + SatN, ParaNum + popt.IFB_num);               // Design MatrixT
		L.Resize(EqNum + SatN, 1);                     // O-C vector
	}
	else
	{
		Q.Resize(EqNum, EqNum);                 // Covariance MatrixT of observation
		P.Resize(EqNum, EqNum);                 // Weight MatrixT of observation

		B.Resize(EqNum, ParaNum + popt.IFB_num);               // Design MatrixT
		L.Resize(EqNum, 1);
	}
	// O-C vector





	m_SatValidN = SatN;                    // Satellite number
	m_TropZenith = 0.0;                    // Zenith troposphere delay

	return 1;
}

int LeastSquarePPP::KF_ParaInt(ppp_option_t& popt, unsigned int SatN)//kf初始化   zhy_text;
{
	unsigned int EqNum = 0, ParaNum = 0;

	// Number of const estimated parameters [X, Y, Z, Rclk, ISB, Trop]
	_ParaN_const = (3 * est_xyz) + est_rclk + _ISBNum + _TropNum;

	if (_MathModel == 4) {
		_ParaN_EstTotal = popt.CUR_mum + SatN + _ParaN_const; // Total estimate parameter number
	}
	else
	{
		_ParaN_EstTotal = _ParaN_sat * SatN + _ParaN_const;
	} // Total estimate parameter number

	ParaNum = _ParaN_EstTotal;              // Total estimated parameter number 参个数
	if (_MathModel == 4)
	{
		EqNum = popt.CUR_mum * 2 + _PsuedoObsN;
	}
	else
	{
		EqNum = _ObsEqNum * SatN + _PsuedoObsN;
	}
	// Equation number (real-observation + psuedo-observation)
	if (popt.ion_const)
	{
		H_zhy.Resize(EqNum + SatN, ParaNum + popt.IFB_num);               // Design MatrixT
		V_zhy.Resize(EqNum + SatN, 1);
		R_zhy.Resize(EqNum + SatN, EqNum + SatN);
		P_zhy.Resize(ParaNum + popt.IFB_num, ParaNum + popt.IFB_num);
		X_zhy.Resize(ParaNum + popt.IFB_num, 1);
	}
	else
	{
		P_zhy.Resize(ParaNum + popt.IFB_num, ParaNum + popt.IFB_num);                 // Covariance MatrixT of observation
		R_zhy.Resize(EqNum, EqNum);                 // Weight MatrixT of observation
		X_zhy.Resize(ParaNum + popt.IFB_num, 1);
		H_zhy.Resize(EqNum, ParaNum + popt.IFB_num);               // Design MatrixT
		V_zhy.Resize(EqNum, 1);
	}
	// O-C vector

	m_SatValidN = SatN;                    // Satellite number
	m_TropZenith = 0.0;                    // Zenith troposphere delay

	return 1;
}

int LeastSquarePPP::Lsq_CalMatrix(ppp_option_t& popt, GenData& gendata, PPPObsData& pppdata, PreciseData& predata, ClockData& clkdata, int flag_zhy)
{
	//% Log variable information
	// std::cout << "Lsq_CalMatrix: Epoch Time = " << popt.CurrentEpoch.m_Year << "-"

	// 清空当前历元的 ion 向量

	MatrixT XX;
	if (flag_zhy)XX = X_Trans; else XX = X_zhy;
	GIM_t iondata = gendata.GIMdata;
	static gnsstime vblt;
	FILE* vbl_fp = NULL;
	if (_traceMode)
	{
		if (!(popt.CurrentEpoch == vblt))
		{
			vblt = popt.CurrentEpoch;
			vbl_fp = fopen(popt.IOfile.out_f.out_vbl_fpath, "a+");
			if (vbl_fp)
				fprintf(vbl_fp, ">%4d-%02d-%02d %02d:%02d:%02d\n",
					vblt.m_Year, vblt.m_Month, vblt.m_Day, vblt.m_Hour, vblt.m_Min, (int)round(vblt.m_Sec));
		}
	}
	//% Output process variable in file end

	unsigned int SatN = 0;              // Satellite number
	//SatN = pppdata.m_SatCount;
	SatN = m_SatValidN;
	m_Prn.clear();                      // Clear cache in m_Prn
	if (m_Slot.size()) m_Slot.clear();  // m_Slot for GLONASS system
	unsigned int mk_zhy = 0;
	unsigned int vsn = 0;               // Valid satellite count
	std::vector<double> ion_var; // 局部变量，用于存储电离层延迟方差
	for (int i_sat = 0; i_sat < SatN; ++i_sat)
	{
		char sys = 0;    char sys_bds = 0;               // Satellite system flag
		unsigned int prn = 0;
		unsigned int satno = pppdata.m_Prn[i_sat];
		if (satno == 0) continue;

		if (pppdata.m_Valid_flag[i_sat] == 0) continue;

		sys = GetSysPrn(satno, prn);

		int orbn = 0;  // Get orbn for GLONASS satellites
		if (sys == 'R') {
			orbn = Get_Nav_Orbitn(satno, predata.m_NavData);
			if (orbn == -1) continue;
		}
		double freq1, freq2, freq3; double wavelength1, wavelength2, wavelength3, wavelength4, wavelength5;

		//多频

		//unsigned int satno = pppdata.m_Prn[i_sat];
		if (GetSystem(satno) == 'G') {
			orbn = 0;
			wavelength1 = Get_WaveLength(satno, 1, orbn);//L1
			wavelength2 = Get_WaveLength(satno, 2, orbn);//L2
			wavelength3 = Get_WaveLength(satno, 5, orbn);//L5
		}
		if (GetSystem(satno) == 'E')
		{
			orbn = 0;
			wavelength1 = Get_WaveLength(satno, 1, orbn);//E1
			wavelength2 = Get_WaveLength(satno, 2, orbn);//E5A
			wavelength3 = Get_WaveLength(satno, 3, orbn);//E5b
			wavelength4 = Get_WaveLength(satno, 4, orbn);//E5(A+B)
			wavelength5 = Get_WaveLength(satno, 5, orbn);//E6
		}
		if (GetSystem(satno) == 'C')
		{
			orbn = 0;
			wavelength1 = Get_WaveLength(satno, 1, orbn);//B1
			wavelength2 = Get_WaveLength(satno, 3, orbn);//B3
			wavelength3 = Get_WaveLength(satno, 2, orbn);//B2
			wavelength4 = Get_WaveLength(satno, 4, orbn);//B1C
			wavelength5 = Get_WaveLength(satno, 5, orbn);//B2a

		}
		if (GetSystem(satno) == 'R')
		{
			orbn = Get_Nav_Orbitn(satno, predata.m_NavData);
			wavelength1 = Get_WaveLength(satno, 1, orbn);//G1
			wavelength2 = Get_WaveLength(satno, 2, orbn);//G2
			wavelength3 = Get_WaveLength(satno, 3, orbn);//G1a
			wavelength4 = Get_WaveLength(satno, 4, orbn);//G2a
			wavelength5 = Get_WaveLength(satno, 5, orbn);//G3
		}



		freq1 = Get_Frequency(satno, 1, orbn);
		/*wavelength1 = Get_WaveLength(satno, 1, orbn);
		if (sys == 'G' || sys == 'E' || sys == 'R') { freq2 = Get_Frequency(satno, 2, orbn); freq3 = Get_Frequency(satno, 5, orbn); wavelength2 = Get_WaveLength(satno, 2, orbn); wavelength3 = Get_WaveLength(satno, 5, orbn); }
		else if (sys == 'C') { freq2 = Get_Frequency(satno, 3, orbn); wavelength2 = Get_WaveLength(satno, 3, orbn);}*/


		gnsstime obst = popt.CurrentEpoch;                        // Current epoch

		double satpos[3] = { 0.0 }, satvel[3] = { 0.0 };          // Satellite Position, Satellite Velocity
		double satpco[3] = { 0.0 }, satclk = 0.0;                 // Satellite PCO correction, Clock offsets
		double Elev = 0.0, Azimuth = 0.0;                         // Satellite Elevation and Azimuth

		for (int m = 0; m < m_SatInfo.m_Size; ++m)
		{
			if (m_SatInfo.m_Prn[m] != satno) continue;

			satpos[0] = m_SatInfo.m_Xs[3 * m + 0];
			satpos[1] = m_SatInfo.m_Xs[3 * m + 1];
			satpos[2] = m_SatInfo.m_Xs[3 * m + 2];
			satvel[0] = m_SatInfo.m_Vs[3 * m + 0];
			satvel[1] = m_SatInfo.m_Vs[3 * m + 1];
			satvel[2] = m_SatInfo.m_Vs[3 * m + 2];
			satclk = m_SatInfo.m_Clk[m];

			Elev = m_SatInfo.m_Ele[m];
			Azimuth = m_SatInfo.m_Azi[m];
		}

		// Get Sat and Recv PCV
		double SatPCV[NFREQ] = { 0.0 }, RecPCV[NFREQ] = { 0.0 };
		double rec_ant_delta[3] = { Antdiff[1], Antdiff[0], Antdiff[2] };
		double azel[2] = { Azimuth, Elev };

		// Satellite PCO and PCV
		if (popt.ssr_com || popt.process_time == MODE_POST)
		{
			SatAntPCO(gendata.m_pcvs, popt.SunPos, satpos, sys, satno, orbn, satpco);//PCO改正

			satpos[0] += satpco[0];
			satpos[1] += satpco[1];
			satpos[2] += satpco[2];

			SatAntPCV(sys, satno, satpos, m_Sta_tmp, gendata.m_pcvs, SatPCV);
		}
		RecAntModel(sys, &gendata.m_pcvr, rec_ant_delta, azel, 1, RecPCV);
		// Get Sat and Recv PCV end

		// Phase Windup correction
		double PhaseWindup1 = 0.0;
		WindupCorr(popt.SunPos, satpos, m_Sta_tmp, popt.phw1[pppdata.m_Prn[i_sat]]);   // L1 Phase Windup
		PhaseWindup1 = popt.phw1[pppdata.m_Prn[i_sat]];

		// Observations: Code and Phase
		double P1 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr1[i_sat]);
		double P2 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr2[i_sat]);
		double P3 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr3[i_sat]);
		double P4 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr4[i_sat]);
		double P5 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr5[i_sat]);



		double L1 = CarrierPhase_PCV_Correct(wavelength1, RecPCV[0], SatPCV[0], pppdata.m_L1[i_sat]);
		double L2 = CarrierPhase_PCV_Correct(wavelength2, RecPCV[0], SatPCV[0], pppdata.m_L2[i_sat]);
		double L3 = CarrierPhase_PCV_Correct(wavelength3, RecPCV[0], SatPCV[0], pppdata.m_L3[i_sat]);
		double L4 = CarrierPhase_PCV_Correct(wavelength4, RecPCV[0], SatPCV[0], pppdata.m_L4[i_sat]);
		double L5 = CarrierPhase_PCV_Correct(wavelength5, RecPCV[0], SatPCV[0], pppdata.m_L5[i_sat]);
		//if (fabs(P1) < 100 || fabs(P2) < 100 || fabs(L1) < 100 || fabs(L2) < 100)continue;

		//double satdcb = 0.0;
		//if (_IS_DCB)  //DCB改正
		//{
		//	if (GetSystem(satno) == 'C')
		//	{
		//		double ptr_dcb = gendata.m_Dcb.Get_DCB_B1B3(pppdata.m_Prn[i_sat]);
		//		satdcb = ptr_dcb * 1E-9* LIGHTSPEED;
		//		
		//		P1 = P1 + satdcb * BDS_B1_B3_gamma;
		//		P2 = P2 + satdcb * BDS_B1_B3_gamma;
		//	}
		//	if (GetSystem(satno) == 'G')
		//	{
		//		double ptr_dcb = gendata.m_Dcb.Get_DCB_P1P2(pppdata.m_Prn[i_sat]);
		//		satdcb = ptr_dcb * 1E-9 * LIGHTSPEED;
		//		
		//		P1 = P1 + satdcb * BDS_L1_L2_gamma;
		//		P2 = P2 - satdcb * BDS_L1_L2_gamma;
		//	}
		//}
		//double L2 = CarrierPhase_PCV_Correct(wavelength2, RecPCV[0], SatPCV[0], pppdata.m_L2[i_sat]);
		L1 *= wavelength1;
		L2 *= wavelength2;						// Cycle to meter (m)
		L3 *= wavelength3;
		L4 *= wavelength4;
		L5 *= wavelength5;
		double rs = 0.0, rr = 0.0, ss = 0.0;
		for (int i = 0; i < 3; ++i)
		{
			rs += (m_Sta_tmp[i] - satpos[i]) * (m_Sta_tmp[i] - satpos[i]);    // (Xr-Xs)^2 +(Yr-Ys)^2+(Zr-Zs)^2
			rr += m_Sta_tmp[i] * m_Sta_tmp[i];                              // xr^2 + yr^2 + zr^2
			ss += satpos[i] * satpos[i];                                    // xs^2 + ys^2 + zs^2
		}
		rr = sqrt(rr); ss = sqrt(ss);
		rs = sqrt(rs);                                                      // sqrt[ (Xr-Xs)^2 +(Yr-Ys)^2+(Zr-Zs)^2 ]

		// Gravitation delay
		double GravitationDelay = 0.0;
		GravitationDelay = GravitationEffect(ss, rr, rs);

		//% Get Troposphere delay
		double TropDelay = 0.0;                                             // Trop slant total delay
		double Trop_dry_z = 0.0, Trop_wet_z = 0.0;                          // Trop zenith dry and wet delay
		double Trop_dry_s = 0, Trop_wet_s = 0;                              // Trop slant dry and wet delay
		double Trop_dry_func = 0.0, Trop_wet_func = 0.0;                    // Trop dry and wet mapping function 

		Coordinate sta_blh;
		sta_blh.XYZ.setXYZ(m_Sta_tmp[0], m_Sta_tmp[1], m_Sta_tmp[2]);
		sta_blh._xyz2blh();
		double blh[3] = { 0.0 };
		blh[0] = sta_blh.BLH._B; blh[1] = sta_blh.BLH._L; blh[2] = sta_blh.BLH._H;
		int doy = obst._GetDoy();

		//Trop_Saastamoinen(sta_blh.BLH._H, Elev, Trop_dry_z, Trop_wet_z, 0);
		//Trop_Hopfield(sta_blh.BLH._H, Elev, Trop_dry_z, Trop_wet_z, 0);
		Trop_UNB3(blh, doy, Elev, Trop_dry_z, Trop_wet_z, 0);
		GMF(obst.m_Mjd, sta_blh.BLH._B, sta_blh.BLH._L, sta_blh.BLH._H, Elev, Trop_dry_func, Trop_wet_func);

		m_TropZenith = Trop_dry_z + Trop_wet_z;
		Trop_dry_s = Trop_dry_z * Trop_dry_func;

		Trop_wet_z = XX(3 * (est_xyz)+est_rclk + _SystemN, 1) + Trop_wet_z;//   对流层传递 for   zhy_text 2021.10.12

		Trop_wet_s = Trop_wet_z * Trop_wet_func;

		TropDelay = Trop_dry_s + Trop_wet_s;
		//% Get Trop delay end
		double IonDelay = 0.0, var_ion = 0.0, S_TEC = 0;
		//dianlicengyanchi 
//	double IonDelay = 0.0, var_ion = 0.0, S_TEC = 0;

			//IonModel_Grid(obst.m_gtime, freq1, &iondata, blh, azel, 1, &IonDelay, &var_ion, &S_TEC);

			//ion.push_back(IonDelay);
			//ion_var.push_back(var_ion);
			// 调试输出
		//	std::cout << "Sat " << satno << ", Freq[" << freq1 << "] = " << freq1
			//	<< ", IonDelay = " << IonDelay << std::endl;
		//试试用p1p2

		//IonDelay = ION_[i_sat];

		// Sagnac effect (if no correct in satellite position)
		double SagnacDelay = 0.0;
		SagnacDelay = SagnacEffect(satpos, m_Sta_tmp);

		// Relative effect
		double RelativeDelay = 0.0;
		RelativeDelay = RelativeEffect(satpos, satvel);

		// Phase windup effect
		double Phawnd = 0.0, Phawnd1 = 0.0, Phawnd2 = 0.0, Phawnd3 = 0.0, Phawnd4 = 0.0, Phawnd5 = 0.0;
		double PHawnd[5] = { 0 }, Wj[5] = { 0 }, Zb[5] = { 0 };
		Phawnd1 = wavelength1 * PhaseWindup1;
		Phawnd2 = wavelength2 * PhaseWindup1;
		Phawnd3 = wavelength3 * PhaseWindup1;
		Phawnd4 = wavelength4 * PhaseWindup1;
		Phawnd5 = wavelength5 * PhaseWindup1;

		// Stochastic model parameter (for Weight Matrix)
		double fact = 1.0;
		double SinEle = sin(Elev);
		double sigemaP = 0.0, sigemaL = 0.0;
		sigemaP = SQR(m_SigemaP); sigemaL = SQR(m_SigemaL);
		if ((sys == 'R') && (_SystemN > 1))
		{
			fact *= 2;
			sigemaP *= 100;
		}
		if ((sys == 'E') && (_SystemN > 1))
			fact *= 1;
		if ((sys == 'C') && (_SystemN > 1))
		{
			if (prn <= 5)
				fact *= 1000;
			else
				fact *= 2;
		}
		if (_IS_ISB)
		{
			sys_bds = GetSystem_BDS_2_3(pppdata.m_Prn[i_sat]);
			sys = sys_bds;
		}
		//  Stochastic model parameter (for Weight Matrix) end
		double BiaS[5] = { 0 }, IonS[5] = { 0 };
		double isb = 0.0; double ion = 0.0, ion1 = 0.0, ion2 = 0.0, ion3 = 0.0, ion4 = 0.0, bias = 0.0, bias1 = 0.0, bias2 = 0.0, bias3 = 0.0, bias4 = 0.0;
		if (_MathModel == MODE_UNCOMBINE || _MathModel == MODE_UNCOMBINE_ION_CONT)
		{
			ion = XX(_ParaN_const + i_sat + 1, 1);
			bias = XX(_ParaN_const + SatN + i_sat + 1, 1);
		}
		else if (_MathModel == MODE_UOFC_COMBINE)
		{
			ion = IonDelay;
			bias = XX(_ParaN_const + i_sat + 1, 1);
		}
		else if (_MathModel == MODE_IF_COMBINE)
		{
			ion = 0.0;
			bias = XX(_ParaN_const + i_sat + 1, 1);
		}
		else if (_MathModel == 4)
		{
			for (int t4 = 0; t4 < num_zhy[i_sat]; t4++)
			{
				if (type_int[i_sat][t4] == 1)
				{
					IonS[t4] = XX(_ParaN_const + i_sat + 1, 1);
				}
				else if (type_int[i_sat][t4] == 2)
				{
					IonS[t4] = XX(_ParaN_const + i_sat + 1, 1) * SQR(wavelength2 / wavelength1);
				}
				else if (type_int[i_sat][t4] == 3)
				{
					IonS[t4] = XX(_ParaN_const + i_sat + 1, 1) * SQR(wavelength3 / wavelength1);
				}
				else if (type_int[i_sat][t4] == 4)
				{
					IonS[t4] = XX(_ParaN_const + i_sat + 1, 1) * SQR(wavelength4 / wavelength1);
				}
				else if (type_int[i_sat][t4] == 5)
				{
					IonS[t4] = XX(_ParaN_const + i_sat + 1, 1) * SQR(wavelength5 / wavelength1);
				}
			}
			for (int t4 = 0; t4 < num_zhy[i_sat]; t4++)
			{
				BiaS[t4] = XX(_ParaN_const + num_zhy_accu[i_sat] + t4 + 1 + SatN, 1);
			}
		}
		if (_SystemN > 1)
		{
			sys = GetSystem_GREC2C3(pppdata.m_Prn[i_sat]);
			if (_SystemIndex[sys] == 1)isb == 0.0;
			else isb = XX(3 * (est_xyz)+_SystemIndex[sys], 1);
		}
		double rclk = XX(3 * (est_xyz)+est_rclk, 1);
		if (_MathModel == 4)
		{
			for (int t4 = 0; t4 < num_zhy[i_sat]; t4++)
			{
				if (type_int[i_sat][t4] == 1)
				{
					Wj[t4] = P1; Zb[t4] = L1; PHawnd[t4] = Phawnd1;
				}
				else if (type_int[i_sat][t4] == 2)
				{
					Wj[t4] = P2; Zb[t4] = L2; PHawnd[t4] = Phawnd2;
				}
				else if (type_int[i_sat][t4] == 3)
				{
					Wj[t4] = P3; Zb[t4] = L3; PHawnd[t4] = Phawnd3;
				}
				else if (type_int[i_sat][t4] == 4)
				{
					Wj[t4] = P4; Zb[t4] = L4; PHawnd[t4] = Phawnd4;
				}
				else if (type_int[i_sat][t4] == 5)
				{
					Wj[t4] = P5; Zb[t4] = L5; PHawnd[t4] = Phawnd5;
				}
			}
		}
		// Output vriables to files
		if (_traceMode)
		{
			if (vbl_fp)
			{
				fprintf(vbl_fp, " %c%02d %15.3f%15.3f%15.3f %7.3f %8.3f %15.3f%15.3f%15.3f %12.4f %12.4f %7.4f %7.4f %7.4f  %7.4f %8.4f %7.4f, %12.4f\n",
					sys, prn, satpos[0], satpos[1], satpos[2], Elev * R2D, Azimuth * R2D,
					P1, L1, rs, satclk * LIGHTSPEED, GravitationDelay, Trop_dry_z, Trop_wet_z, RelativeDelay, Phawnd1, bias);
			}
		}
		// Output variables to files end ----------------
		unsigned int mk = 0;
		if (_MathModel == 4)
		{
			if (i_sat == 0)
			{
				mk_zhy = vsn * num_zhy[i_sat] * 2 + _PsuedoObsN;
				//mk_zhy = vsn * num_zhy[i_sat]  + _PsuedoObsN;
			}
			else
			{
				mk_zhy = num_zhy[i_sat - 1] * 2 + mk_zhy + _PsuedoObsN;
				//mk_zhy = num_zhy[i_sat - 1] * 2 +1+ mk_zhy + _PsuedoObsN;
			}
		}
		else
		{
			mk = vsn * _ObsEqNum + _PsuedoObsN;
		}
		/// Undifferenced uncombined model
		if (_IS_ISB)
		{
			sys_bds = GetSystem_BDS_2_3(pppdata.m_Prn[i_sat]);
			sys = sys_bds;
		}// FOR  测试  北斗
		if (_MathModel == MODE_IF_COMBINE)//***************************************************************************  无电离层组合
		{
			double C1 = SQR(wavelength2) / (SQR(wavelength2) - SQR(wavelength1));    //这里的波长已经判断完事了
			double C2 = -SQR(wavelength1) / (SQR(wavelength2) - SQR(wavelength1));
			double P_ = C1 * P1 + C2 * P2;
			double L_ = C1 * L1 - Phawnd1 * wavelength1 + C2 * L2 - Phawnd2 * wavelength2;
			// OMC
			L(mk + 1, 1) = P_ - rs + satclk * LIGHTSPEED + RelativeDelay + GravitationDelay - TropDelay - SagnacDelay - ion - rclk - isb;
			L(mk + 2, 1) = L_ - rs + satclk * LIGHTSPEED + RelativeDelay + GravitationDelay - TropDelay - SagnacDelay - bias - rclk - isb;    //ph相位缠绕；
			double a = L(mk + 1, 1);
			double b = L(mk + 2, 1);
			if (flag_epoch_count >= 1313)
			{
				a = b;
			}
			// Coefficient of X,Y,Z
			if (est_xyz)
			{
				B(mk + 1, 1) = (m_Sta_tmp[0] - satpos[0]) / rs;
				B(mk + 2, 1) = (m_Sta_tmp[0] - satpos[0]) / rs;
				B(mk + 1, 2) = (m_Sta_tmp[1] - satpos[1]) / rs;
				B(mk + 2, 2) = (m_Sta_tmp[1] - satpos[1]) / rs;
				B(mk + 1, 3) = (m_Sta_tmp[2] - satpos[2]) / rs;
				B(mk + 2, 3) = (m_Sta_tmp[2] - satpos[2]) / rs;
			}
			//接收机钟差
			B(mk + 1, 3 * est_xyz + est_rclk) = 1.0;
			B(mk + 2, 3 * est_xyz + est_rclk) = 1.0;
			if (_SystemN > 1)//ISB
			{
				B(mk + 1, 3 * est_xyz + _SystemIndex[sys]) = 1.0;
				B(mk + 2, 3 * est_xyz + _SystemIndex[sys]) = 1.0;
			}
			// Coefficient of troposphere
			B(mk + 1, 3 * est_xyz + est_rclk * _SystemN + est_trop) = Trop_wet_func;
			B(mk + 2, 3 * est_xyz + est_rclk * _SystemN + est_trop) = Trop_wet_func;
			// Coefficient of ambiguity
			B(mk + 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + SatN * _IonNum + _PsuedoObsN + vsn + 1) = 1;
			// Stochastic model (Note: Users can modified their own models)
			//sigemaL= SQR(m_SigemaL);
			Q(mk + 1, mk + 1) = fact * (sigemaP * 3) / (SinEle * SinEle);
			Q(mk + 2, mk + 2) = fact * (sigemaL * 3) / (SinEle * SinEle);
		}
		if (_MathModel == MODE_UOFC_COMBINE)//***************************************************************************   半和改正模型
		{
			double L_P = 0.0;
			L_P = P1 + L1; L_P /= 2;
			// OMC
			L(mk + 1, 1) = P1 - rs + satclk * LIGHTSPEED + RelativeDelay + GravitationDelay - TropDelay - SagnacDelay - ion - rclk - isb;
			L(mk + 2, 1) = L_P - rs + satclk * LIGHTSPEED + RelativeDelay + GravitationDelay - TropDelay - SagnacDelay - Phawnd1 / 2 - bias - rclk - isb;    //ph相位缠绕；
			double a = L(mk + 1, 1);
			double b = L(mk + 2, 1);
			// Coefficient of X,Y,Z
			if (est_xyz)
			{
				B(mk + 1, 1) = (m_Sta_tmp[0] - satpos[0]) / rs;
				B(mk + 2, 1) = (m_Sta_tmp[0] - satpos[0]) / rs;
				B(mk + 1, 2) = (m_Sta_tmp[1] - satpos[1]) / rs;
				B(mk + 2, 2) = (m_Sta_tmp[1] - satpos[1]) / rs;
				B(mk + 1, 3) = (m_Sta_tmp[2] - satpos[2]) / rs;
				B(mk + 2, 3) = (m_Sta_tmp[2] - satpos[2]) / rs;
			}
			//接收机钟差
			B(mk + 1, 3 * est_xyz + est_rclk) = 1.0;
			B(mk + 2, 3 * est_xyz + est_rclk) = 1.0;
			if (_SystemN > 1)//ISB
			{
				B(mk + 1, 3 * est_xyz + _SystemIndex[sys]) = 1.0;
				B(mk + 2, 3 * est_xyz + _SystemIndex[sys]) = 1.0;
			}
			// Coefficient of troposphere
			B(mk + 1, 3 * est_xyz + est_rclk * _SystemN + est_trop) = Trop_wet_func;
			B(mk + 2, 3 * est_xyz + est_rclk * _SystemN + est_trop) = Trop_wet_func;
			// Coefficient of ambiguity
			B(mk + 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + SatN * _IonNum + _PsuedoObsN + vsn + 1) = 0.5;
			// Stochastic model (Note: Users can modified their own models)
			sigemaP = SQR(m_SigemaP); sigemaL = SQR((m_SigemaL + m_SigemaP) / 2);
			//sigemaL= SQR(m_SigemaL);
			double  aa = fact * sigemaP / (SinEle * SinEle);
			double  bb = fact * (sigemaL) / (SinEle * SinEle);
			Q(mk + 1, mk + 1) = fact * sigemaP / (SinEle * SinEle);
			Q(mk + 2, mk + 2) = fact * (sigemaL) / (SinEle * SinEle);
		}
		if (_MathModel == MODE_UNCOMBINE)//***********************************************************************************   非组合
		{
			isb = 0; rclk = 0;
			if (_IS_ISB) {
				if (sys == 'B')isb = XX(_SystemIndex[sys] + 3 * (est_xyz), 1);//B2
				else rclk = XX(_SystemIndex[sys] + 3 * (est_xyz), 1);
			}
			else
			{
				if (_SystemN > 1)
				{
					if (_SystemIndex[sys] == 1)isb == 0.0;
					else isb = XX(3 * (est_xyz)+_SystemIndex[sys], 1);
				}
			}
			double rclk = XX(3 * (est_xyz)+est_rclk, 1);
			// OMC
			L(mk + 1, 1) = P1 - rs + satclk * LIGHTSPEED + RelativeDelay + GravitationDelay - TropDelay - SagnacDelay - ion - rclk - isb;
			L(mk + 2, 1) = L1 - rs + satclk * LIGHTSPEED + RelativeDelay + GravitationDelay - TropDelay - SagnacDelay - Phawnd1 - bias + ion - rclk - isb;//ph相位缠绕；
			double a = L(mk + 1, 1);
			double b = L(mk + 2, 1);
			// Coefficient of X,Y,Z
			if (est_xyz)
			{
				B(mk + 1, 1) = (m_Sta_tmp[0] - satpos[0]) / rs;
				B(mk + 2, 1) = (m_Sta_tmp[0] - satpos[0]) / rs;

				B(mk + 1, 2) = (m_Sta_tmp[1] - satpos[1]) / rs;
				B(mk + 2, 2) = (m_Sta_tmp[1] - satpos[1]) / rs;

				B(mk + 1, 3) = (m_Sta_tmp[2] - satpos[2]) / rs;
				B(mk + 2, 3) = (m_Sta_tmp[2] - satpos[2]) / rs;
			}
			if (!(_IS_ISB))
			{
				B(mk + 1, 3 * est_xyz + est_rclk) = 1.0;
				B(mk + 2, 3 * est_xyz + est_rclk) = 1.0;
			}
			if (_SystemN > 1)//ISB
			{
				B(mk + 1, 3 * est_xyz + _SystemIndex[sys]) = 1.0;
				B(mk + 2, 3 * est_xyz + _SystemIndex[sys]) = 1.0;
			}
			// Coefficient of troposphere
			B(mk + 1, 3 * est_xyz + est_rclk * _SystemN + est_trop) = Trop_wet_func;
			B(mk + 2, 3 * est_xyz + est_rclk * _SystemN + est_trop) = Trop_wet_func;
			// Coefficient of Ionosphere
			if (_IonModel == ION_MODEL_EST_WHITE || _IonModel == ION_MODEL_EST_RANDWALK)
			{
				B(mk + 1, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = 1.0;
				B(mk + 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = -1.0;
			}
			// Coefficient of ambiguity
			B(mk + 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + SatN + _PsuedoObsN + vsn + 1) = 1.0;
			// Stochastic model (Note: Users can modified their own models)
			Q(mk + 1, mk + 1) = fact * sigemaP / (SinEle * SinEle);
			Q(mk + 2, mk + 2) = fact * sigemaL / (SinEle * SinEle);
		}
		if (_MathModel == MODE_UNCOMBINE_ION_CONT)//**************************************************************************   电离层约束
		{
			double isb = 0;
			double ion = XX(_ParaN_const + i_sat + 1, 1);
			double bias = XX(_ParaN_const + SatN + i_sat + 1, 1);
			double rclk = XX(3 * (est_xyz)+est_rclk, 1);
			// O-C
			if (_SystemIndex[sys] == 1)
			{
				isb = 0;
			}
			else
			{
				isb = XX(_SystemIndex[sys] + 3, 1);
			}

			L(mk + 1, 1) = P1 - rs + satclk * LIGHTSPEED + RelativeDelay + GravitationDelay - TropDelay - SagnacDelay - ion - rclk + isb;
			L(mk + 2, 1) = L1 - rs + satclk * LIGHTSPEED + RelativeDelay + GravitationDelay - TropDelay - SagnacDelay - Phawnd1 - bias + ion - rclk + isb;	    //ph相位缠绕；
			L(mk + 3, 1) = ION_[vsn] - ion;																														//电离层约束；

			double a = L(mk + 1, 1);
			double b = L(mk + 2, 1);
			double c = L(mk + 3, 1);
			// Coefficient of X,Y,Z
			if (est_xyz)
			{
				B(mk + 1, 1) = (m_Sta_tmp[0] - satpos[0]) / rs;
				B(mk + 2, 1) = (m_Sta_tmp[0] - satpos[0]) / rs;

				B(mk + 1, 2) = (m_Sta_tmp[1] - satpos[1]) / rs;
				B(mk + 2, 2) = (m_Sta_tmp[1] - satpos[1]) / rs;

				B(mk + 1, 3) = (m_Sta_tmp[2] - satpos[2]) / rs;
				B(mk + 2, 3) = (m_Sta_tmp[2] - satpos[2]) / rs;
			}

			// Coefficient of receiver clock and ISB
			/*if (flag_epoch_count == 0)
			{
				B(mk + 1, 3 * est_xyz + _SystemIndex[sys]) = 0.0;
				B(mk + 2, 3 * est_xyz + _SystemIndex[sys]) = 0.0;
			}
			else
			{*/
			B(mk + 1, 3 * est_xyz + est_rclk) = 1.0;
			B(mk + 2, 3 * est_xyz + est_rclk) = 1.0;

			if (_SystemN > 1)
			{
				B(mk + 1, 3 * est_xyz + _SystemIndex[sys]) = 1.0;
				B(mk + 2, 3 * est_xyz + _SystemIndex[sys]) = 1.0;
			}

			// Coefficient of troposphere
			B(mk + 1, 3 * est_xyz + est_rclk * _SystemN + est_trop) = Trop_wet_func;
			B(mk + 2, 3 * est_xyz + est_rclk * _SystemN + est_trop) = Trop_wet_func;

			// Coefficient of Ionosphere
			if (_IonModel == ION_MODEL_EST_WHITE)
			{
				B(mk + 1, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = 1.0;
				B(mk + 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = -1.0;
				B(mk + 3, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = 1;
			}

			// Coefficient of ambiguity
			B(mk + 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + SatN + _PsuedoObsN + vsn + 1) = 1.0;

			// Stochastic model (Note: Users can modified their own models)
			Q(mk + 1, mk + 1) = fact * sigemaP / (SinEle * SinEle);
			Q(mk + 2, mk + 2) = fact * sigemaL / (SinEle * SinEle);
			int zhy_Q = 2;
			if (flag_epoch_count > 10)zhy_Q = 200;
			Q(mk + 3, mk + 3) = 1 + zhy_Q * (flag_epoch_count - 1) * 30;//0.09
		}
		if (_MathModel == 4)//***********************************************************************************  500块
		{
			isb = 0; rclk = 0; double ifb = 0;
			if (_IS_ISB) {
				if (sys == 'B')isb = XX(_SystemIndex[sys] + 3 * (est_xyz), 1);//B2
				else rclk = XX(_SystemIndex[sys] + 3 * (est_xyz), 1);
			}
			else
			{
				if (_SystemN > 1)
				{
					if (_SystemIndex[sys] == 1)isb == 0.0;
					else isb = XX(3 * (est_xyz)+_SystemIndex[sys], 1);
				}
			}

			double rclk = XX(3 * (est_xyz)+est_rclk, 1);
			// OMC
			int ifb_i1 = 0;
			for (int PLS = 0; PLS < num_zhy[i_sat]; PLS++)
			{
				if (get_fcbindex(type_[i_sat][PLS + 1]))
				{
					ifb_i1 = get_fcbindex(type_[i_sat][PLS + 1]); double Cifb = XX(3 * est_xyz + est_rclk * _SystemN + _TropNum + SatN + _PsuedoObsN + popt.CUR_mum + _IFB_Index[ifb_i1], 1);
					ifb = Cifb;
				}

				L(mk_zhy + 2 * PLS + 1, 1) = Wj[PLS] - rs + satclk * LIGHTSPEED + RelativeDelay + GravitationDelay - TropDelay - SagnacDelay - IonS[PLS] - rclk - isb - ifb;
				L(mk_zhy + 2 * PLS + 2, 1) = Zb[PLS] - rs + satclk * LIGHTSPEED + RelativeDelay + GravitationDelay - TropDelay - SagnacDelay - PHawnd[PLS] - BiaS[PLS] + IonS[PLS] - rclk - isb;//ph相位缠绕；
				double a1 = L(mk_zhy + 2 * PLS + 1, 1);
				double b1 = L(mk_zhy + 2 * PLS + 2, 1);

			}
			//L(mk_zhy + 2 * num_zhy[i_sat]+1, 1) = IonDelay;
			//L.Output_to_File("E:\\L.txt", 0);
			// Coefficient of X,Y,Z
			if (est_xyz)
			{
				for (int t1 = 0; t1 < num_zhy[i_sat] * 2; t1++)	// t1 代表的是 行数
				{
					for (int t2 = 0; t2 < 3; t2++)	// t2 代表的是 x y z
					{
						B(mk_zhy + 1 + t1, 1 + t2) = (m_Sta_tmp[t2] - satpos[t2]) / rs;
					}
				}
			}
			if (!(_IS_ISB))
			{
				for (int t2 = 0; t2 < num_zhy[i_sat] * 2; t2++)
				{
					B(mk_zhy + 1 + t2, 3 * est_xyz + est_rclk) = 1.0;
				}
			}
			if (_SystemN > 1)//ISB
			{
				for (int t2 = 0; t2 < num_zhy[i_sat] * 2; t2++)
				{
					B(mk_zhy + 1 + t2, 3 * est_xyz + _SystemIndex[sys]) = 1.0;
				}
			}
			if (popt.IFB_num >= 1)//IFB
			{
				int t1 = 0; int ifb_i = 0;
				for (int t2 = 1; t2 < num_zhy[i_sat] * 2 + 1; t2++)
				{
					if (get_fcbindex(type_[i_sat][t1 + 1])) {
						ifb_i = get_fcbindex(type_[i_sat][t1 + 1]); double Cifb = t2 % 2 == 0 ? 0 : get_fcbcofe(type_[i_sat][t1 + 1]);
						B(mk_zhy + t2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + SatN + _PsuedoObsN + popt.CUR_mum + _IFB_Index[ifb_i]) = Cifb;
					}
					t1 = t2 % 2 == 0 ? (t1 + 1) : t1;//判断t2是不是偶数
				}
			}
			// Coefficient of troposphere
			for (int t2 = 0; t2 < num_zhy[i_sat] * 2; t2++)
			{
				B(mk_zhy + 1 + t2, 3 * est_xyz + est_rclk * _SystemN + est_trop) = Trop_wet_func;
			}

			// Coefficient of Ionosphere
			if (_IonModel == ION_MODEL_EST_WHITE || _IonModel == ION_MODEL_EST_RANDWALK)
			{
				double xishu[10] = { 1, -1,SQR(wavelength2 / wavelength1), -SQR(wavelength2 / wavelength1), SQR(wavelength3 / wavelength1), -SQR(wavelength3 / wavelength1), SQR(wavelength4 / wavelength1), -SQR(wavelength4 / wavelength1), SQR(wavelength5 / wavelength1), -SQR(wavelength5 / wavelength1) };
				for (int t2 = 0; t2 < num_zhy[i_sat]; t2++)
				{
					if (type_int[i_sat][t2] == 1)
					{
						B(mk_zhy + 1 + t2 * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = xishu[0];
						B(mk_zhy + 2 + t2 * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = xishu[1];
					}
					else if (type_int[i_sat][t2] == 2)
					{
						B(mk_zhy + 1 + t2 * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = xishu[2];
						B(mk_zhy + 2 + t2 * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = xishu[3];
					}
					else if (type_int[i_sat][t2] == 3)
					{
						B(mk_zhy + 1 + t2 * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = xishu[4];
						B(mk_zhy + 2 + t2 * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = xishu[5];
					}
					else if (type_int[i_sat][t2] == 4)
					{
						B(mk_zhy + 1 + t2 * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = xishu[6];
						B(mk_zhy + 2 + t2 * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = xishu[7];
					}
					else if (type_int[i_sat][t2] == 5)
					{
						B(mk_zhy + 1 + t2 * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = xishu[8];
						B(mk_zhy + 2 + t2 * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = xishu[9];
					}
				}
				//B(mk_zhy + num_zhy[i_sat] * 2+1, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + vsn + 1) = 1;
				//B.Output_to_File("E:\\B.txt", 0);
			}

			// Coefficient of ambiguity
			for (int t2 = 0; t2 < num_zhy[i_sat]; t2++)
			{

				B(mk_zhy + (1 + t2) * 2, 3 * est_xyz + est_rclk * _SystemN + _TropNum + SatN + _PsuedoObsN + num_zhy_accu[i_sat] + t2 + 1) = 1.0;

			}
			//B.Output_to_File("E:\\B.txt", 0);
			// Stochastic model (Note: Users can modified their own models)
			int add_1 = 0;
			for (int var_dif = 0; var_dif < num_zhy[i_sat]; var_dif++)
			{

				double type_var_l = var_dif > 1 ? Get_dif_fre_var(satno, type_int[i_sat][var_dif], orbn) : 1;
				double type_var_p = var_dif > 1 ? 10 : 1;
				Q(mk_zhy + 1 + 2 * add_1, mk_zhy + 1 + 2 * add_1) = (fact * sigemaP / (SinEle * SinEle)) * SQR(type_var_p);
				Q(mk_zhy + 2 + 2 * add_1, mk_zhy + 2 + 2 * add_1) = (fact * sigemaL / (SinEle * SinEle)) / SQR(type_var_l);
				add_1++;
			}
			//Q.Output_to_File("E:\\Q.txt", 0);
			add_1 = 0;
#if 0


			if (num_zhy[i_sat] == 2)
			{
				Q(mk_zhy + 1, mk_zhy + 1) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 2, mk_zhy + 2) = fact * sigemaL / (SinEle * SinEle);
				Q(mk_zhy + 3, mk_zhy + 3) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 4, mk_zhy + 4) = fact * sigemaL / (SinEle * SinEle);
				//Q(mk_zhy + 5, mk_zhy + 5) = 1.0/var_ion;
			//	Q.Output_to_File("E:\\Q.txt",0);
			}
			else if (num_zhy[i_sat] == 3)
			{
				Q(mk_zhy + 1, mk_zhy + 1) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 2, mk_zhy + 2) = fact * sigemaL / (SinEle * SinEle);
				Q(mk_zhy + 3, mk_zhy + 3) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 4, mk_zhy + 4) = fact * sigemaL / (SinEle * SinEle);
				Q(mk_zhy + 5, mk_zhy + 5) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 6, mk_zhy + 6) = fact * sigemaL / (SinEle * SinEle);
				//Q(mk_zhy + 7, mk_zhy + 7) = 1.0 / var_ion;
			//	Q.Output_to_File("E:\\Q.txt", 0);
			}
			else if (num_zhy[i_sat] == 4)
			{
				Q(mk_zhy + 1, mk_zhy + 1) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 2, mk_zhy + 2) = fact * sigemaL / (SinEle * SinEle);
				Q(mk_zhy + 3, mk_zhy + 3) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 4, mk_zhy + 4) = fact * sigemaL / (SinEle * SinEle);
				Q(mk_zhy + 5, mk_zhy + 5) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 6, mk_zhy + 6) = fact * sigemaL / (SinEle * SinEle);
				Q(mk_zhy + 7, mk_zhy + 7) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 8, mk_zhy + 8) = fact * sigemaL / (SinEle * SinEle);
				//Q(mk_zhy + 10, mk_zhy + 10) = 1.0 / var_ion;
				//	Q.Output_to_File("E:\\Q.txt", 0);
			}
			else if (num_zhy[i_sat] == 5)
			{
				Q(mk_zhy + 1, mk_zhy + 1) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 2, mk_zhy + 2) = fact * sigemaL / (SinEle * SinEle);
				Q(mk_zhy + 3, mk_zhy + 3) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 4, mk_zhy + 4) = fact * sigemaL / (SinEle * SinEle);
				Q(mk_zhy + 5, mk_zhy + 5) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 6, mk_zhy + 6) = fact * sigemaL / (SinEle * SinEle);
				Q(mk_zhy + 7, mk_zhy + 7) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 8, mk_zhy + 8) = fact * sigemaL / (SinEle * SinEle);
				Q(mk_zhy + 9, mk_zhy + 9) = fact * sigemaP / (SinEle * SinEle);
				Q(mk_zhy + 10, mk_zhy + 10) = fact * sigemaL / (SinEle * SinEle);
				//Q(mk_zhy + 11, mk_zhy + 11) = 1.0 / var_ion;
				//	Q.Output_to_File("E:\\Q.txt", 0);
			}

#endif 
		}
		//Q.Output_to_File("E:\\Q.txt", 0);
		m_Prn.push_back(satno);
		m_Slot.push_back(orbn);

		++vsn;
	}

	if (popt.ion_const)
	{
		for (int i = 0; i < SatN; i++)
		{
			B(popt.CUR_mum * 2 + i + 1, 3 * est_xyz + est_rclk * _SystemN + _TropNum + _PsuedoObsN + i + 1) = 1;
			L(popt.CUR_mum * 2 + i + 1, 1) = ion[i];
			Q(popt.CUR_mum * 2 + i + 1, popt.CUR_mum * 2 + i + 1) = 1 / ion_var[i] * 2;
		}
	}

	//L.Output_to_File("E:\\L.txt", 0);
	//B.Output_to_File("E:\\B.txt", 0);
	//Q.Output_to_File("E:\\Q.txt", 0);
	m_SatValidN = vsn;


	P = Q.InvCholesky();//高度角定权；
	// 保存当前历元的 ion 到 ion_prev


	// 调试输出（可选）

	if (vbl_fp) fclose(vbl_fp);
	return 1;
}

/// Generate Least Square MatrixT
int LeastSquarePPP::Lsq_LSMatrix()
{
	if (flag_epoch_count == 1)  // for the 1st epoch   第一个历元
	{
		if (_FreqNum == 1)
		{
			L_x = V_zhy;
			P_x = R_zhy;
			B_x = H_zhy;
			_Prn_previous.clear();
			for (int i = 0; i < m_SatValidN; ++i)   // Save prn in previous epoch
				_Prn_previous.push_back(m_Prn[i]);
		}
	}
	else if ((flag_epoch_count == 2) && (_FreqNum == 1))//  第二个历元
	{
		// Find common satellite between first 2 epoch for sf-uncombined ppp
		vector<unsigned int> Prn_Common, Prn_index_previous, Prn_index_current;
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
						Prn_Common.push_back(m_Prn[i]);
						Prn_index_previous.push_back(j);
						Prn_index_current.push_back(i);
					}
				}
			}
		}

		_CommSatN = Prn_Common.size();                  // Common satellte number with previous epoch
		if (_CommSatN == 0) return 0;

		unsigned int obs_equationN = 0;                 // Number of equations
		obs_equationN = _CommSatN * _ObsEqNum * 2;
		L_LS.Initialize(obs_equationN, 1);               // O-C 
		P_LS.Initialize(obs_equationN, obs_equationN);   // Observation variance


		unsigned int est_paraN_total = 0;
		if (_IonModel == ION_MODEL_EST_WHITE)
			est_paraN_total = _ParaN_const + _CommSatN * _IonNum * 2 + _CommSatN * _AmbNum - _RefSatValidFlag;
		B_LS.Initialize(obs_equationN, est_paraN_total);
		unsigned int IonPsuedoN_Previous = 0;
		double a = 0; double b = 0;
		// Construct LS Calculate MatrixT: L_LS, B_LS, P_LS
		for (int i = 0; i < _CommSatN; ++i) // Previous epoch (1st epoch)
		{
			for (int j = 0; j < _Prn_previous.size(); ++j)
			{
				if (Prn_Common[i] == _Prn_previous[j])
				{
					for (int k = 0; k < _ObsEqNum; ++k)
					{
						L_LS(_PsuedoObsN + i * _ObsEqNum + k + 1, 1) = L_x(IonPsuedoN_Previous + j * _ObsEqNum + k + 1, 1);


						for (int m = 0; m < _ParaN_const; ++m)
						{
							B_LS(_PsuedoObsN + i * _ObsEqNum + k + 1, m + 1) = B_x(IonPsuedoN_Previous + j * _ObsEqNum + k + 1, m + 1);
						}

						// Ion
						if (_IonModel == ION_MODEL_EST_WHITE)
							B_LS(_PsuedoObsN + i * _ObsEqNum + k + 1, _ParaN_const + _PsuedoObsN + i + 1) =
							B_x(IonPsuedoN_Previous + j * _ObsEqNum + k + 1, _ParaN_const + IonPsuedoN_Previous + j + 1);
						//@@@@@@@@@@@

						//@@@@@@@@@@
						P_LS(_PsuedoObsN + i * _ObsEqNum + k + 1, _PsuedoObsN + i * _ObsEqNum + k + 1) =
							P_x(IonPsuedoN_Previous + j * _ObsEqNum + k + 1, IonPsuedoN_Previous + j * _ObsEqNum + k + 1);
						if (Prn_Common[i] == _RefSatPrn && _RefSatValidFlag)
							continue;

						// Amb
						if (Prn_Common[i] == _RefSatPrn && _RefSatValidFlag)
							continue;

						if (_IonModel == ION_MODEL_EST_WHITE)
							B_LS(i * _ObsEqNum + k + 1, _ParaN_const + 2 * _CommSatN * _IonNum + i + 1) =
							B_x(j * _ObsEqNum + k + 1, _ParaN_const + _IonNum * _Prn_previous.size() + j + 1);

						//@@@@@@@@@@@

					 //@@@@@@@@@@@@@@

					}

					break;
				}
			}

		}
		for (int i = 0; i < _CommSatN; ++i) // Current epoch (2nd epoch)
		{
			for (int j = 0; j < m_SatValidN; ++j)
			{
				if (Prn_Common[i] == m_Prn[j])
				{
					for (int k = 0; k < _ObsEqNum; ++k)
					{
						L_LS(_CommSatN * _ObsEqNum + 2 * _PsuedoObsN + i * _ObsEqNum + k + 1, 1) = V_zhy(_PsuedoObsN + j * _ObsEqNum + k + 1, 1);

						for (int m = 0; m < _ParaN_const; ++m)
						{
							B_LS(_CommSatN * _ObsEqNum + 2 * _PsuedoObsN + k + i * _ObsEqNum + 1, m + 1) = H_zhy(_PsuedoObsN + j * _ObsEqNum + k + 1, m + 1);
						}

						// Ion
						if (_IonModel == ION_MODEL_EST_WHITE)
							B_LS(_CommSatN * _ObsEqNum + 2 * _PsuedoObsN + i * _ObsEqNum + k + 1, _ParaN_const + _PsuedoObsN + _CommSatN + i + 1)
							= H_zhy(j * _ObsEqNum + k + 1 + _PsuedoObsN, _ParaN_const + _PsuedoObsN + j + 1);

						// Qxx
						P_LS(_CommSatN * _ObsEqNum + 2 * _PsuedoObsN + i * _ObsEqNum + k + 1, _CommSatN * _ObsEqNum + 2 * _PsuedoObsN + i * _ObsEqNum + k + 1)
							= R_zhy(j * _ObsEqNum + k + 1 + _PsuedoObsN, j * _ObsEqNum + k + 1 + _PsuedoObsN);

						if (Prn_Common[i] == _RefSatPrn && _RefSatValidFlag)
							continue;

						// Amb
						if (_IonModel == ION_MODEL_EST_WHITE)
							B_LS(_CommSatN * _ObsEqNum + 2 * _PsuedoObsN + i * _ObsEqNum + k + 1, _ParaN_const + _PsuedoObsN + 2 * _CommSatN * _IonNum + i + 1)
							= H_zhy(j * _ObsEqNum + k + 1 + _PsuedoObsN, _ParaN_const + _PsuedoObsN + m_SatValidN * _IonNum + j + 1);
					}

					break;
				}
			}
		}
		V_zhy.Initialize(obs_equationN, 1);               // O-C 
		R_zhy.Initialize(obs_equationN, obs_equationN);
		H_zhy.Initialize(obs_equationN, est_paraN_total);
		V_zhy = L_LS;
		R_zhy = P_LS;
		H_zhy = B_LS;
		L_LS.Clear();
		P_LS.Clear();
		B_LS.Clear();
		//H_zhy.Output_to_File("E:\\H", 1);

	}


	return 1;
}

/// Save information for next epoch
/// X,Y,Z,Trop,(Dcb'),Amb
int LeastSquarePPP::Lsq_FiltingMatrix()
{
	_ParaN_iondcb = 0;

	/// Psuedo Observation value ------------------------------------------------
	if (_IonModel == ION_MODEL_EST_WHITE)
		L_x.Resize(XYZ_Rclk_Trop_Num + m_SatValidN * _AmbNum + m_SatValidN * _IonNum, 1);  // No GIM and S-dcb corrected(ion,R-dcb,S-dcb,(1st epoch R-clk))  Time variant

	// X, Y, Z (only static mode)
	if (_Mode == MODE_STATIC) {
		for (int i = 0; i < (3 * est_xyz); ++i)
			L_x(i + 1, 1) = X(i + 1, 1);
	}
	//X.display();
	//Nbb.display();
	//L_x.display();
	//P_x.display();
	// Troposphere
	if (est_trop)
	{
		L_x(XYZ_Rclk_Num + 1, 1) = X(3 * est_xyz + est_rclk * _SystemN + est_trop, 1);
	}

	// Ambiguity
	for (int i = 0; i < m_SatValidN * _AmbNum; ++i)
	{
		if ((_MathModel == MODE_UNCOMBINE || _MathModel == MODE_UNCOMBINE_ION_CONT) && (_FreqNum == 1) && (flag_epoch_count == 1))
		{
			if (_IonModel == ION_MODEL_EST_WHITE)
				L_x(XYZ_Rclk_Trop_Num + i + m_SatValidN + 1, 1) =
				X(_ParaN_const + _PsuedoObsN + _IonNum * _CommSatN * 2 + i + 1, 1);
		}
		else
		{
			L_x(XYZ_Rclk_Trop_Num + i + m_SatValidN + 1, 1) =
				X(_ParaN_const + _PsuedoObsN + _IonNum * m_SatValidN + i + 1, 1);
		}
	}
	//ION  zhy_text
	for (int i = 0; i < m_SatValidN * _IonNum; ++i)
	{
		if ((_MathModel == MODE_UNCOMBINE || _MathModel == MODE_UNCOMBINE_ION_CONT) && (_FreqNum == 1) && (flag_epoch_count == 1))
		{
			if (_IonModel == ION_MODEL_EST_WHITE)
				L_x(XYZ_Rclk_Trop_Num + i + 1, 1) =
				X(_ParaN_const + _PsuedoObsN + _IonNum * _CommSatN + i + 1, 1);
		}
		else
		{
			L_x(XYZ_Rclk_Trop_Num + i + 1, 1) =
				X(_ParaN_const + _PsuedoObsN + i + 1, 1);
		}
	}
	/// Psuedo observation end --------------------------------------------------

	// Psuedo observation power start -------------------------------------------
	P_x.Resize(XYZ_Rclk_Trop_Num + _AmbNum * m_SatValidN + m_SatValidN, XYZ_Rclk_Trop_Num + _AmbNum * m_SatValidN + m_SatValidN);

	unsigned int m = 0, fact = 1;
	if (((_MathModel == MODE_UNCOMBINE || _MathModel == MODE_UNCOMBINE_ION_CONT)) && (_FreqNum == 1) && (flag_epoch_count == 1))
		fact = 2;

	// Qxx filter
	if (flag_epoch_count == 1) {
		for (int i = 0; i < Nbb.GetRow(); ++i)
		{
			if (_Mode == MODE_STATIC)
			{
				// skip receiver clk and isb parameter
				if ((i >= (3 * est_xyz + _PsuedoRclkNum)) && (i < 3 * est_xyz + est_rclk * _SystemN))
					continue;
			}
			if (_Mode == MODE_KINEMATIC)
			{
				// skip position, receiver clk and isb parameter
				if (i < 3 * est_xyz + est_rclk * _SystemN)
					continue;
			}
			if ((i <= _ParaN_const - 1 + _IonNum * _CommSatN) && (i > _ParaN_const - 1)) continue;
			//if (_IonModel == ION_MODEL_EST_WHITE)
			//{
				// skip ionosphere parameter
				//if ((i >= (3 * est_xyz + est_rclk * _SystemN + _TropNum)) && (i < (3 * est_xyz + _PsuedoObsN + est_rclk * _SystemN + _TropNum + m_SatValidN * fact)))
				//	continue;
			//}

			unsigned int n = 0;

			for (int j = 0; j < Nbb.GetCol(); ++j)
			{
				if (_Mode == MODE_STATIC)
				{
					// skip receiver clk and isb parameter
					if ((j >= (3 * est_xyz + _PsuedoRclkNum)) && (j < 3 * est_xyz + est_rclk * _SystemN))
						continue;
				}
				if (_Mode == MODE_KINEMATIC)
				{
					// skip position, receiver clk and isb parameter
					if (j < 3 * est_xyz + est_rclk * _SystemN)
						continue;
				}
				if ((j <= _ParaN_const - 1 + _IonNum * _CommSatN) && (j > _ParaN_const - 1)) continue;
				//if (_IonModel == ION_MODEL_EST_WHITE)
				//{
					// skip ionosphere parameter
				//	if ((j >= (3 * est_xyz + est_rclk * _SystemN + _TropNum)) && (j < (3 * est_xyz + _PsuedoObsN + est_rclk * _SystemN + _TropNum + m_SatValidN * fact)))
				//		continue;
				//}

				P_x(m + 1, n + 1) = Nbb(i + 1, j + 1);

				n++;
			}
			m++;
		}
	}
	else {

		for (int i = 0; i < Nbb.GetRow(); ++i)
		{
			if (_Mode == MODE_STATIC)
			{
				// skip receiver clk and isb parameter
				if ((i >= (3 * est_xyz + _PsuedoRclkNum)) && (i < 3 * est_xyz + est_rclk * _SystemN))
					continue;
			}
			if (_Mode == MODE_KINEMATIC)
			{
				// skip position, receiver clk and isb parameter
				if (i < 3 * est_xyz + est_rclk * _SystemN)
					continue;
			}

			//if (_IonModel == ION_MODEL_EST_WHITE)
			//{
				// skip ionosphere parameter
				//if ((i >= (3 * est_xyz + est_rclk * _SystemN + _TropNum)) && (i < (3 * est_xyz + _PsuedoObsN + est_rclk * _SystemN + _TropNum + m_SatValidN * fact)))
				//	continue;
			//}

			unsigned int n = 0;

			for (int j = 0; j < Nbb.GetCol(); ++j)
			{
				if (_Mode == MODE_STATIC)
				{
					// skip receiver clk and isb parameter
					if ((j >= (3 * est_xyz + _PsuedoRclkNum)) && (j < 3 * est_xyz + est_rclk * _SystemN))
						continue;
				}
				if (_Mode == MODE_KINEMATIC)
				{
					// skip position, receiver clk and isb parameter
					if (j < 3 * est_xyz + est_rclk * _SystemN)
						continue;
				}

				//if (_IonModel == ION_MODEL_EST_WHITE)
				//{
					// skip ionosphere parameter
				//	if ((j >= (3 * est_xyz + est_rclk * _SystemN + _TropNum)) && (j < (3 * est_xyz + _PsuedoObsN + est_rclk * _SystemN + _TropNum + m_SatValidN * fact)))
				//		continue;
				//}

				P_x(m + 1, n + 1) = Nbb(i + 1, j + 1);

				n++;
			}
			m++;
		}
	}
	// Psuedo observation power end ---------------------------------------------

	/// Random Walk Noise -------------------------------------------------------对流层随机游走
	double delta_t = 0, ion_noise = 0.0, trop_noise = 0.0;
	if (m_dt)	delta_t = m_dt;
	else	delta_t = 30.0;
	delta_t = abs(delta_t);


	if (est_trop)
	{
		if (_Mode == MODE_STATIC)
		{
			P_x(XYZ_Rclk_Num + 1, XYZ_Rclk_Num + 1) += _TropProcessNoise * _TropProcessNoise * delta_t;
		}
		else if (_Mode == MODE_KINEMATIC)
		{
			P_x(_PsuedoRclkNum + 1, _PsuedoRclkNum + 1) += _TropProcessNoise * _TropProcessNoise * delta_t;
		}
	}
	//P_ION zhy_text

	if (est_trop)
	{
		if (_Mode == MODE_STATIC)
		{
			for (int i = 0; i < m_SatValidN * _IonNum; ++i)
			{
				P_x(XYZ_Rclk_Num + 1 + i, XYZ_Rclk_Num + 1 + i) +=
					SQR(_IonProcessNoise) * delta_t;
				//0.04* 0.04 * delta_t;//电离层随机游走  不对    完事还有协方差没加上

			}
		}
		else if (_Mode == MODE_KINEMATIC)
		{
			for (int i = 0; i < m_SatValidN * _IonNum; ++i)
			{
				P_x(XYZ_Rclk_Num + 1 + i, XYZ_Rclk_Num + 1 + i) +=
					_TropProcessNoise * _TropProcessNoise * delta_t;

			}
		}

	}
	// --------------------------------------------------------------------------
	//L_x.display();
	//P_x.display();

	//L_x.Output_to_File("E:\\L_x.out", 1);m     
	//P_x.Output_to_File("E:\\P_x.out", 1);
	return 1;
}

int LeastSquarePPP::Lsq_StateUpdate(ppp_option_t& popt)
{

	if (est_xyz)
	{
		m_Sta[0] = X_Trans(1, 1);
		m_Sta[1] = X_Trans(2, 1);
		m_Sta[2] = X_Trans(3, 1);
	}

	if (est_trop)
	{
		// Zenith trop delay
		m_TropZenith += X_Trans(3 * est_xyz + est_rclk * _SystemN + est_trop, 1);
	}
	double var0, var, var1;
	if (_MathModel == MODE_IF_COMBINE || _MathModel == 4);
	{
		if ((_Cycle_Sliptest.size() + _Cycle_Sliptest13.size()) > 0)
		{
			for (int i = 0; i < m_Prn.size(); i++)
			{
				tLast[m_Prn[i]] = popt.CurrentEpoch;
				gf_pre[m_Prn[i]] = gf_cur[m_Prn[i]];
				gf_pre13[m_Prn[i]] = gf_cur13[m_Prn[i]];
				for (int k = 0; k < _Cycle_Sliptest.size(); k++)
				{
					if (_Cycle_Sliptest[k] == m_Prn[i])
					{
						mw_count[m_Prn[i]] = 0;
					}
				}
				for (int k = 0; k < _Cycle_Sliptest13.size(); k++)
				{
					if (_Cycle_Sliptest13[k] == m_Prn[i])
					{
						mw_count[m_Prn[i]] = 0;
					}
				}

				if (mw_count[m_Prn[i]] > 0)
				{
					int j = mw_count[m_Prn[i]];
					mw_pre[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_pre13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					var0 = mw_var[m_Prn[i]];
					var1 = SQR(mw_cur[m_Prn[i]] - mw_avg[m_Prn[i]]) - var0;
					var1 = var0 + var1 / j;
					mw_avg[m_Prn[i]] = (mw_avg[m_Prn[i]] * j + mw_cur[m_Prn[i]]) / (j + 1);
					mw_avg13[m_Prn[i]] = (mw_avg13[m_Prn[i]] * j + mw_cur13[m_Prn[i]]) / (j + 1);
					mw_count[m_Prn[i]]++;
					mw_var[m_Prn[i]] = var1;
				}
				else
				{
					mw_pre[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_avg[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_pre13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					mw_avg13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					mw_count[m_Prn[i]]++;
					mw_var[m_Prn[i]] = 0.25;
				}
			}
		}
		else
		{
			for (int i = 0; i < m_Prn.size(); i++)
			{
				tLast[m_Prn[i]] = popt.CurrentEpoch;
				gf_pre[m_Prn[i]] = gf_cur[m_Prn[i]];
				gf_pre13[m_Prn[i]] = gf_cur13[m_Prn[i]];
				if (mw_count[m_Prn[i]] > 0)
				{
					int j = mw_count[m_Prn[i]];
					mw_pre[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_pre13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					var0 = mw_var[m_Prn[i]];
					var1 = SQR(mw_cur[m_Prn[i]] - mw_avg[m_Prn[i]]) - var0;
					var1 = var0 + var1 / j;
					mw_avg[m_Prn[i]] = (mw_avg[m_Prn[i]] * j + mw_cur[m_Prn[i]]) / (j + 1);
					mw_avg13[m_Prn[i]] = (mw_avg13[m_Prn[i]] * j + mw_cur13[m_Prn[i]]) / (j + 1);
					mw_count[m_Prn[i]]++;
					mw_var[m_Prn[i]] = var1;
				}
				else
				{
					mw_pre[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_avg[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_pre13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					mw_avg13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					mw_count[m_Prn[i]]++;
					mw_var[m_Prn[i]] = 0.25;
				}
			}
		}
	}

	GetSatNum(m_Prn, _GpsSatNum, _GloSatNum, _BdsSatNum, _GalSatNum);

	_Prn_previous.clear();
	for (int i = 0; i < m_SatValidN; ++i)
		_Prn_previous.push_back(m_Prn[i]);

	return 1;
}

int LeastSquarePPP::Lsq_StateUpdateSPP(ppp_option_t& popt)
{
	double var0, var, var1;
	if (_MathModel == MODE_IF_COMBINE || _MathModel == 4);
	{
		if ((_Cycle_Sliptest.size() + _Cycle_Sliptest13.size()) > 0)
		{
			for (int i = 0; i < m_Prn.size(); i++)
			{
				tLast[m_Prn[i]] = popt.CurrentEpoch;
				gf_pre[m_Prn[i]] = gf_cur[m_Prn[i]];
				gf_pre13[m_Prn[i]] = gf_cur13[m_Prn[i]];
				for (int k = 0; k < _Cycle_Sliptest.size(); k++)
				{
					if (_Cycle_Sliptest[k] == m_Prn[i])
					{
						mw_count[m_Prn[i]] = 0;
					}
				}
				for (int k = 0; k < _Cycle_Sliptest13.size(); k++)
				{
					if (_Cycle_Sliptest13[k] == m_Prn[i])
					{
						mw_count[m_Prn[i]] = 0;
					}
				}

				if (mw_count[m_Prn[i]] > 0)
				{
					int j = mw_count[m_Prn[i]];
					mw_pre[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_pre13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					var0 = mw_var[m_Prn[i]];
					var1 = SQR(mw_cur[m_Prn[i]] - mw_avg[m_Prn[i]]) - var0;
					var1 = var0 + var1 / j;
					mw_avg[m_Prn[i]] = (mw_avg[m_Prn[i]] * j + mw_cur[m_Prn[i]]) / (j + 1);
					mw_avg13[m_Prn[i]] = (mw_avg13[m_Prn[i]] * j + mw_cur13[m_Prn[i]]) / (j + 1);
					mw_count[m_Prn[i]]++;
					mw_var[m_Prn[i]] = var1;
				}
				else
				{
					mw_pre[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_avg[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_pre13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					mw_avg13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					mw_count[m_Prn[i]]++;
					mw_var[m_Prn[i]] = 0.25;
				}
			}
		}
		else
		{
			for (int i = 0; i < m_Prn.size(); i++)
			{
				tLast[m_Prn[i]] = popt.CurrentEpoch;
				gf_pre[m_Prn[i]] = gf_cur[m_Prn[i]];
				gf_pre13[m_Prn[i]] = gf_cur13[m_Prn[i]];
				if (mw_count[m_Prn[i]] > 0)
				{
					int j = mw_count[m_Prn[i]];
					mw_pre[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_pre13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					var0 = mw_var[m_Prn[i]];
					var1 = SQR(mw_cur[m_Prn[i]] - mw_avg[m_Prn[i]]) - var0;
					var1 = var0 + var1 / j;
					mw_avg[m_Prn[i]] = (mw_avg[m_Prn[i]] * j + mw_cur[m_Prn[i]]) / (j + 1);
					mw_avg13[m_Prn[i]] = (mw_avg13[m_Prn[i]] * j + mw_cur13[m_Prn[i]]) / (j + 1);
					mw_count[m_Prn[i]]++;
					mw_var[m_Prn[i]] = var1;
				}
				else
				{
					mw_pre[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_avg[m_Prn[i]] = mw_cur[m_Prn[i]];
					mw_pre13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					mw_avg13[m_Prn[i]] = mw_cur13[m_Prn[i]];
					mw_count[m_Prn[i]]++;
					mw_var[m_Prn[i]] = 0.25;
				}
			}
		}
	}

	return 1;
}

int LeastSquarePPP::Lsq_Filter(ppp_option_t& popt, gnsstime& t, GenData& gendata, ObsEpochData& obsdata, PreciseData* predata, ClockData* clkdata, vector<zhydata>& prepppdata, vector<zhyinfo>& presatinfo, int zhy_falg, deque<vector<int>>& compareprn, deque<vector<double>>& iongim123)
{
	//% SPP part =========================================================
	vector<double> ion_gim, ion_gimsub, STEC, isb, trop, isbsub; vector<int> compareprn1, commprn; int tempcount = 0;
	for (int q = 0; q < 4; q++) { popt.ISBsystem_pre[q] = popt.ISBsystem[q]; popt.ISB_pre[q] = popt.ISB[q]; popt.isbsit_pre[q] = popt.isbsit[q]; }
	vector<double> trop_prev;
	vector<double> satclks;
	vector<double> satclks_prev;
	vector<double> trop_spp; 
	static std::deque<std::vector<double>> trop_deque;
	static std::deque<std::vector<double>> satclks_deque; 
	vector<unsigned int> satclks_prn;
	vector<unsigned int> satclks_prn_prev;
	static std::deque<std::vector<unsigned int>> satclks_prn_deque; 
	if (CodeSPP(popt, obsdata, gendata, predata, clkdata,
		m_SatInfo, ion_gim, STEC, trop,
		satclks, satclks_prev,
		satclks_prn, satclks_prn_prev) == 0)
		
	{
		_Prn_previous.clear(); num_zhy_pre.clear();
		return 0;
	}
	trop_spp = trop; 
	satclks_prn_deque.push_back(satclks_prn);
	for (int q = 0; q < m_SatInfo.m_Size; q++)
	{
		compareprn1.push_back(m_SatInfo.m_Prn[q]);
		switch ((m_SatInfo.m_Prn[q] / 100))
		{
		case 0:_GpsSatNum++; break;
		case 1:_GloSatNum++; break;
		case 2:_GalSatNum++; break;
		case 3:_BdsSatNum++; break;
		default:break;
		}
	}
	for (int q = 0; q < 4; q++)
	{
		if ((popt.ISBsystem[q] == popt.ISBsystem_pre[q]) && (popt.ISBsystem_pre[q]) == 'C' || popt.ISBsystem_pre[q] == 'R' || popt.ISBsystem_pre[q] == 'E' || popt.ISBsystem_pre[q] == 'B')
		{
			isb.push_back(popt.ISB[popt.isbsit[q] - 1] - popt.ISB_pre[popt.isbsit_pre[q] - 1]); tempcount++;
		}
	}


	satclks_deque.push_back(satclks); 
	satclks_prn_deque.push_back(satclks_prn);  
	trop_deque.push_back(trop_spp);
	compareprn.push_back(compareprn1);
	iongim123.push_back(ion_gim);

	if (satclks_deque.size() >= 2) {
		satclks_prev = satclks_deque[satclks_deque.size() - 2];
		satclks = satclks_deque.back();
		satclks_prn_prev = satclks_prn_deque[satclks_prn_deque.size() - 2];  
		satclks_prn = satclks_prn_deque.back();                          
	}
	else {
		satclks = satclks_deque.back();
		satclks_prn = satclks_prn_deque.back(); 
		satclks_prev.assign(satclks.size(), 0.0);
		satclks_prn_prev.assign(satclks_prn.size(), 0);  
	}

	
	if (trop_deque.size() >= 2) {
		trop_prev = trop_deque[trop_deque.size() - 2];
		trop = trop_deque.back();
		
	}
	else {
		trop = trop_deque.back();
		trop_prev.assign(trop.size(), 0.0);
	
	}
	while (compareprn.size() > 2) compareprn.pop_front();
	while (iongim123.size() > 2) iongim123.pop_front();
	while (trop_deque.size() > 2) trop_deque.pop_front();
	while (satclks_deque.size() > 2) satclks_deque.pop_front();
	while (satclks_prn_deque.size() > 2) satclks_prn_deque.pop_front();  


	ofstream ofiongimsub;
	if (0)
	{
		ofiongimsub.open("C:\\Users\\15149\\Desktop\\shiy\\ofiongimsub.out", ios::app);
		if (flag_epoch_count >= 2)
		{
			for (int j = 0; j < ion_gimsub.size(); j++)
			{
				if (j == 0)
				{
					ofiongimsub << setw(6) << setfill(' ') << flag_epoch_count;
				}
				ofiongimsub << setw(18) << setfill(' ') << fixed << setprecision(8) << ion_gimsub[j];
				if (j == (ion_gimsub.size() - 1))
				{
					ofiongimsub << endl;
				}
			}
		}
	}

	//% SPP end ===========================================================
	printf("%4d", flag_epoch_count);
	//////=================== PPP Part ====================================

	if (!Lsq_Option_Init(popt))  return 0;    

	// Get Observations for ppp
	int SatCount = obsdata.m_SatCount;
	zhydata prepppdata_ = { 0 }; zhyinfo presatinfo_ = { 0 };

	PPPObsData pppdata(SatCount);//
	pppdata.m_obst = obsdata.m_EpochTime;


	SatCount = pppdata.GetBlockData(popt, &gendata, obsdata, predata, clkdata, -1);

	for (int i = 0; i < pppdata.m_SatCount; ++i)
	{
		double dmp[3] = { 0 };
		BDSMultipathCorr(pppdata.m_Prn[i], *(m_SatInfo.m_Ele + i), dmp);//北斗二星端多路径改正
		*pppdata.m_Pr1 += dmp[0];
		*pppdata.m_Pr2 += dmp[2];
		*pppdata.m_Pr3 += dmp[1];
	}
	gnsstime tutc = obsdata.m_EpochTime;
	double sun[3] = { 0 }, moon[3] = { 0 }, gmst = 0;
	GetSunMoonPos(tutc, sun, moon, &gmst);               // Get sun and moon position
	for (int i = 0; i < 3; ++i) {
		popt.SunPos[i] = sun[i];
		popt.MoonPos[i] = moon[i];
	}
	RC_Repair(popt, pppdata, predata->m_NavData, m_SatInfo.m_Ele, flag_epoch_count);
	
	for (int ll = 0; ll < 3; ll++)
	{
		CoarseXYZ_pre[ll] = CoarseXYZ[ll];
	}

	if (flag_epoch_count >= 1 && est_xyz)
	{
		if (_Mode == 0)//动态
		{
			for (int zhy_ = 0; zhy_ < 3; zhy_++)
			{
				CoarseXYZ[zhy_] = popt.StaPos_SPP[zhy_];
			}
		}
		else if (_Mode == 1)//静态
		{
			for (int zhy_ = 0; zhy_ < 3; zhy_++)
			{
				CoarseXYZ[zhy_] = X_Trans(zhy_ + 1, 1);
			}
		}
		else return -1;
	}


	TideDisp(tutc, CoarseXYZ, 1, &gendata.m_Erp, gendata.OceanDisp, _TideCorrection);
	for (int i = 0; i < 3; ++i)	m_Sta_tmp[i] = CoarseXYZ[i] + _TideCorrection[i];

	// Delete Sat lower elevation threshold

	for (int i = 0; i < pppdata.m_SatCount; i++)
	{
		int NUM = *(pppdata._num + i);
		int ACCU = *(pppdata._num_accu + i);
		num_zhy.push_back(NUM);
		num_zhy_accu.push_back(ACCU);
	}

	for (int i = 0; i < pppdata.m_SatCount; ++i)
	{
		int j = 0;
		for (j = 0; j < m_SatInfo.m_Size; ++j) {
			if (pppdata.m_Prn[i] == m_SatInfo.m_Prn[j])
				break;
		}
		if (j == m_SatInfo.m_Size) {
			if (_MathModel == 4)
			{
				popt.CUR_mum = popt.CUR_mum - num_zhy[i];
				DeleteSatData_num(pppdata, i--);
			}
			else
			{
				DeleteSatData(pppdata, i--);
			}
		}
	}

	DelSat_ElevThresold(popt, m_SatInfo, pppdata);
	for (int j = 0; j < m_SatInfo.m_Size; j++)
	{
		presatinfo_.m_Prn[j] = *(m_SatInfo.m_Prn + j);
		presatinfo_.m_Xs[j * 3 + 0] = *(m_SatInfo.m_Xs + j * 3);
		presatinfo_.m_Xs[j * 3 + 1] = *(m_SatInfo.m_Xs + j * 3 + 1);
		presatinfo_.m_Xs[j * 3 + 2] = *(m_SatInfo.m_Xs + j * 3 + 2);
		presatinfo_.Ele[j] = *(m_SatInfo.m_Ele + j);
	}

	for (int i = 0; i < pppdata.m_SatCount; i++)
	{
		prepppdata_.m_L1[i] = *(pppdata.m_L1 + i);
		prepppdata_.m_L2[i] = *(pppdata.m_L2 + i);
		prepppdata_.m_L3[i] = *(pppdata.m_L3 + i);
		prepppdata_.m_L4[i] = *(pppdata.m_L4 + i);
		prepppdata_.m_L5[i] = *(pppdata.m_L5 + i);
		prepppdata_.m_P1[i] = *(pppdata.m_Pr1 + i);
		prepppdata_.m_P2[i] = *(pppdata.m_Pr2 + i);
		prepppdata_.m_P3[i] = *(pppdata.m_Pr3 + i);
		prepppdata_.m_P4[i] = *(pppdata.m_Pr4 + i);
		prepppdata_.m_P5[i] = *(pppdata.m_Pr5 + i);
		prepppdata_.m_Prn[i] = *(pppdata.m_Prn + i);

	}
	presatinfo_.m_Size = m_SatInfo.m_Size;
	prepppdata_.m_SatCount = pppdata.m_SatCount;
	prepppdata.push_back(prepppdata_);//存入观测值；
	presatinfo.push_back(presatinfo_);//存入卫星坐标；
	SetISBindex(pppdata);

	num_zhy_pre = num_zhy;
	type_pre = type_;
	type_int_pre = type_int;

	type_ = pppdata.Prn_Otype;//传递观测类型
	type_int = pppdata.Prn_Otype_int;
	type_pre = type_;
	type_int_pre = type_int;

#if 0
	if (flag_epoch_count >= 1) QC_zhy_CycleSlip(prepppdata, pppdata);
#endif 

	
#if 1

	if (_MathModel == MODE_IF_COMBINE || _MathModel == 4)                                                    
	{
		QC_detect(popt, prepppdata, pppdata, *predata, gendata); Lsq_StateUpdateSPP(popt);  
	}
#endif
	/*std::cout << "TDCP 前: ion.size=" << ion.size() << ", ion_prev.size=" << ion_prev.size()
		<< ", iongim123.size=" << iongim123.size() << std::endl;*/

#if 1
	if (flag_epoch_count >= 1 && num_zhy_pre.size() != 0) {
		
		if (prepppdata.size() < 2) {
			std::cout << "[Cycle_slip] 跳过：prepppdata.size()="
				<< prepppdata.size() << " < 2，历元="
				<< flag_epoch_count << std::endl;
		}
		else if (popt.schemeset < 3) {
			Cycle_slip(popt, t, prepppdata, presatinfo, pppdata, *predata,
				ion_gimsub, commprn, isb, trop, trop_prev,
				satclks, satclks_prev,
				satclks_prn, satclks_prn_prev); 
		}
		else {
			Cycle_slipTDCP(popt, prepppdata, presatinfo, pppdata, *predata,
				ion_gimsub, commprn);
		}
	}
#endif

	_Cycle_Sliptest.clear();
	_Cycle_Sliptest13.clear();
	_Cycle_Sliptestgf.clear();
	_Cycle_Sliptestgf13.clear();

	if (popt.process_mode == MODE_SPP)
	{
		m_SatValidN = ion_gimsub.size();
		LsqSPP_Output(popt);
		flag_epoch_count++;
		_GpsSatNum = 0; _GloSatNum = 0; _BdsSatNum = 0; _GalSatNum = 0;
		return 1;
	}

	type_ = pppdata.Prn_Otype;
	type_int = pppdata.Prn_Otype_int;

	SatCount = pppdata.m_SatCount;
	if (SatCount == 0) {
		// Log in file
		FILE* logfp = fopen(popt.IOfile.out_f.out_log_fpath, "a+");
		if (logfp)
			fprintf(logfp, "  Warrning: Get no data in GetBlockData(). \n");
		if (logfp) fclose(logfp);
		// log end

		return 0;
	}
	if (RankDetect(popt, SatCount))	return 0;

	
	unsigned int itern = 1, QC_flag = 0; int itern2 = 4;//
	for (int i = 0; (i < itern) && (itern <= ITERNMAX); ++i)
	{
		SatCount = pppdata.m_SatCount;
		if (RankDetect(popt, SatCount)) return 0;//

		cal_IFB_num(popt, pppdata);

		SetISBindex(pppdata);                 // Get system number and set system index

		Lsq_ParaInit(popt, SatCount);     // Parameter of PPP intialize

		KF_ParaInt(popt, SatCount);    // Parameter of PPP intialize  KFint (X,H,V,P,R);

		if ((flag_epoch_count == 0) && (_FreqNum == 1)) // SF PPP Start from 2nd epoch
		{
			flag_epoch_count++;
			//return 2;
		}
		popt.ion_const = 0;
		//IonModel_Grid(obst.m_gtime, freq, &iondata, blh3, azel, 3, &IonDelay, &var_ion);
		UdState2KF(&popt, m_SatValidN, pppdata, ion_gim, *predata);//kf updata zhy_text 2021/10/6; // Generate Matrix for EKF//

		_Cycle_Slip.clear();

		Lsq_CalMatrix(popt, gendata, pppdata, *predata, *clkdata, 0);//当前历元 系数阵 权证 pco pcv

		H_zhy = B; R_zhy = P.InvCholesky(); V_zhy = L;//赋值

		if (RankDetect(popt, m_SatValidN)) return 0;
		//V_zhy.Output_to_File("E:\\V_zhy_0.out", 0);
		//H_zhy.Output_to_File("E:\\H_zhy.out", 0);
#if	0	

		{
			X_zhy.Output_to_File("E:\\X_zhy.out", 0);
			H_zhy.Output_to_File("E:\\H_zhy.out", 0);
			P_zhy.Output_to_File("E:\\P_zhy.out", 0);
			R_zhy.Output_to_File("E:\\R_zhy.out", 0);
			V_zhy.Output_to_File("E:\\V_zhy.out", 0);
		}

#endif //	qw
		bool KEY = false;
		for (int kf = 0; kf < itern2; kf++)
		{

			XIgema0 = Kfilter(X_zhy, H_zhy, P_zhy, R_zhy, V_zhy, X_Trans, P_Trans, flag_epoch_count);
#if 1		


			if (est_xyz)
			{
				for (int i = 0; i < 3; i++)   CoarseXYZ[i] = X_Trans(i + 1, 1);
			}



			TideDisp(tutc, CoarseXYZ, 1, &gendata.m_Erp, gendata.OceanDisp, _TideCorrection);
			for (int i = 0; i < 3; ++i)	m_Sta_tmp[i] = CoarseXYZ[i] + _TideCorrection[i];
			
			std::vector<double> trop_temp = trop; 
			std::vector<double> satclks_temp = satclks; 
			Lsq_CalMatrix(popt, gendata, pppdata, *predata, *clkdata, 1);
			trop = trop_temp;
			satclks = satclks_temp; 
			H_zhy = B; R_zhy = P.InvCholesky();
			V_zhy = L; X_zhy = X_Trans; P_zhy = P_Trans;
		
			if (QC_flag)
			{
				//if (flag_epoch_count < 0)break;
				switch (Residual_editing(V_zhy, R_zhy, pppdata))
				{
				case 0: {
					break;
				}
				case 1: {

					itern2++;
				}
				default:
					break;
				}
				//itern2++;
			}
			//R_zhy.Output_to_File("E:\\R_zhy.out", 1);
			if (itern2 > 1) break;
#endif // 0

		}
		ion_gim.clear();
#if 0 


		{
			X_Trans.Output_to_File("E:\\X_Trans.out", 0);
			P_Trans.Output_to_File("E:\\P_Trans.out", 0);
		}

#endif
	}
	if (itern > ITERNMAX) return 0;

	_previousEpoch = popt.CurrentEpoch;

	//Lsq_FiltingMatrix();

	Lsq_StateUpdate(popt);   //更新MW,平滑GF
	if (popt.mode_filter == 2 && zhy_falg)
	{
		if (_OutputFile) {
			LsqPPP_Output(popt, pppdata);
		}
	}
	else if (popt.mode_filter != 2)
	{
		if (_OutputFile)
			LsqPPP_Output(popt, pppdata);
	}

	_Cycle_Sliptest.clear();
	_Cycle_Sliptest13.clear();
	_Cycle_Sliptestgf.clear();
	_Cycle_Sliptestgf13.clear();

	LsqEpochClear();                // Clear parameter in current epoch
	num_zhy_pre = num_zhy;
	type_pre = type_;
	type_int_pre = type_int;
	num_zhy_accu_pre = num_zhy_accu;
	_ParaN_const_pre = _ParaN_const;
	_SystemN_pre = _SystemN;
	popt.PRE_num = popt.CUR_mum;

	type_int.clear();
	type_.clear();
	num_zhy.clear();
	num_zhy_accu.clear();
	flag_epoch_count++;

	////// PPP end ============================================================

	// ========================================================================
	// Plot results for realtime ppp ------------------------------------------
	//if (_Pmode == MODE_PPP && _TimeMode == MODE_REALTIME && popt.Gnuplot)

	if (popt.mode_filter == 2 && zhy_falg) {
		if (_Pmode == MODE_PPP && popt.Gnuplot)
			if (flag_epoch_count % 20 == 0)
				Plot_SUPREME_Post_Results(&popt);
	}
	else if (popt.mode_filter != 2)
	{
		if (_Pmode == MODE_PPP && popt.Gnuplot)
			if (flag_epoch_count % 20 == 0)
				Plot_SUPREME_Post_Results(&popt);
	}

	return 1;
}

int LeastSquarePPP::LsqEpochClear()
{
	_CommSatN = 0;
	_PsuedoObsNum = 0;

	Q.Clear();
	P.Clear();
	B.Clear();
	L.Clear();


	X_zhy.Clear();
	P_zhy.Clear();
	R_zhy.Clear();
	V_zhy.Clear();


	m_Prn.clear();
	m_SatInfo.Clear();
	//_Cycle_Slip.clear();

	return 1;
}

int LeastSquarePPP::Clear()
{
	m_Prn.clear();
	m_SatInfo.Clear();

	B.Clear(); P.Clear(); Q.Clear(); L.Clear();
	B_LS.Clear(); P_LS.Clear(); L_LS.Clear(); Nbb.Clear();
	V.Clear(); X.Clear();
	B_x.Clear(); P_x.Clear(); L_x.Clear();
	Ion_X->dgf_value.clear();
	Ion_X->gf_value.clear();
	Ion_X->time.clear();
	Ion_X->ele.clear();
	num_zhy.clear();
	num_zhy_accu.clear();
	num_zhy_accu_pre.clear();
	_Prn_previous.clear();
	_Prn_common.clear();
	BIAS_.clear();
	ION_.clear();

	return 1;
}

/// Least Square Solution ---------------
/* Parameter:
*            nn = bt*p*b
*            x = (bt*p*b)_-1 * bt*p*l
*            v = b*x-l
* -------------------------------------- */
int LeastSquare(MatrixT& b, MatrixT& l, MatrixT& p, MatrixT& v, MatrixT& nn, MatrixT& x)
{//间接平差
	MatrixT bt(b.GetCol(), b.GetRow());
	MatrixT btp(b.GetCol(), b.GetRow());
	MatrixT pp(b.GetRow(), b.GetRow());

	pp = p;

	bt = b.Trans();
	btp = bt * pp;
	nn = btp * b;                // Nbb = BTPB

	nn = nn.InvCholesky();     // Nbb_-1 = BTPB_-1

	MatrixT nbt(nn.GetRow(), bt.GetCol());
	MatrixT nbtp(nn.GetRow(), b.GetRow());

	nbt = nn * bt;

	nbtp = nbt * pp;
	x = nbtp * l;        // BTPB_-1 * BTPL

	v = b * x - l;       // V = B * X - L
	return 1;
}

int LeastSquareQcc(MatrixT& b, MatrixT& l, MatrixT& p, MatrixT& v, MatrixT& nn, MatrixT& cc, MatrixT& x)
{
	//间接平差
	MatrixT bt(b.GetCol(), b.GetRow()); MatrixT btp(b.GetCol(), b.GetRow()); MatrixT pp(b.GetRow(), b.GetRow());

	pp = p;
	bt = b.Trans();
	btp = bt * pp;
	nn = btp * b;                // Nbb = BTPB
	nn = nn.InvCholesky();     // Nbb_-1 = BTPB_-1

	MatrixT nbt(nn.GetRow(), bt.GetCol()); MatrixT nbtp(nn.GetRow(), b.GetRow());

	nbt = nn * bt;
	nbtp = nbt * pp;

	x = nbtp * l;        // BTPB_-1 * BTPL
	v = b * x - l;       // V = B * X - L

	pp = p.InvCholesky();		//Q
	cc = pp - (b * nn * bt);		//Qvv

	return 1;
}

int LeastSquareconstraint(MatrixT& b, MatrixT& l, MatrixT& p, MatrixT& c, MatrixT& w, MatrixT& v, MatrixT& vz, MatrixT& nn, MatrixT& cc, MatrixT& x, MatrixT& Pxx)
{
	MatrixT bt(b.GetCol(), b.GetRow()); MatrixT ct(c.GetCol(), c.GetRow());
	MatrixT btp(b.GetCol(), b.GetRow());
	MatrixT pp(b.GetRow(), b.GetRow());

	pp = p;

	bt = b.Trans();
	btp = bt * pp;
	nn = btp * b + Pxx;                // Nbb = BTPB+Pxx

	nn = nn.InvCholesky();		// Nbb_-1

	pp = p.InvCholesky();		//Q
	cc = pp - (b * nn * bt);		//Qvv

	x = nn * (btp * l + Pxx * w);
	//x = (nn - nn * ct * cc * c * nn) * btp * l - nn * ct * cc * w;        // (Nbb_-1-Nbb_-1*CT*Ncc_-1*C*Nbb_-1)*W-Nbb_-1*CT*Ncc_-1*Wx
	v = b * x - l;       // V = B * X - L
	vz = c * x - w;
	return 1;
}

double Kfilter(MatrixT& X, MatrixT& H, MatrixT& Pp, MatrixT& R, MatrixT& V, MatrixT& X_Trans, MatrixT& P_Trans, unsigned int flag_epoch_count)
{
	MatrixT Kk(H.GetCol(), H.GetRow());
	MatrixT I; I = I.Eye(H.GetCol());
	MatrixT F(H.GetCol(), H.GetCol());
	MatrixT Q(R.GetCol(), R.GetCol());
	MatrixT PP, xigema0; double Xigema0 = 0.0; double vk = 0.0, rk = 0.0, ZHY = 0.0;
	int nub = V.GetRow() - X.GetRow(); double c0 = 3, c1 = 8, a = 0.0;
	bool key = true;
	if (flag_epoch_count > 40 && !(key))
	{
		for (int i = 0; i < nub; i++)
		{
			vk += sqrt(V(i + 1, 1));
			rk += sqrt(R(i + 1, i + 1));
		}
		ZHY = vk / rk;
		if (ZHY <= c0) a = 1;
		if (ZHY > c0 && ZHY < c1)a = (c0 / fabs(vk)) * pow(((c1 - fabs(vk)) / (c1 - c0)), 2);
		if (ZHY >= c1) a = 10;
		Pp = Pp * a;
	}

	Q = H * Pp * H.Trans() + R;
	//Q.Output_to_File("E:\\Q_zhy.out", 1);
	Q = Q.Inv();
	//Q.Output_to_File("E:\\Q_zhy_1.out", 1);
	Kk = Pp * H.Trans() * Q;
	//Kk.Output_to_File("E:\\Kk_zhy.out", 1);
	X_Trans = X + Kk * V;
	F = (I - Kk * H);
	P_Trans = F * Pp * F.Trans();
	P_Trans = P_Trans + Kk * R * Kk.Trans();

	MatrixT v(V.GetRow(), 1);
	v = H * Kk * V - V;
	//v.Output_to_File("E:\\v.dat", 0);

	PP = (R + H * P_Trans * H.Trans()).Inv();
	xigema0 = (v.Trans() * R.Inv() * v) / nub;
	Xigema0 = xigema0.GetValue(1, 1);
	//Pk = Pp + Kk * R * Kk.Trans() - Kk * H * Pp * H.Trans() * Kk.Trans();
	Q.Clear();
	F.Clear();
	Kk.Clear();
	I.Clear();
	return sqrt(Xigema0);
}
int LeastSquarePPP::Residual_editing(MatrixT& v, MatrixT& r, PPPObsData& pppdata)//郭斐_改进的抗差卡尔曼滤波方案 zhy_text 2011.11.11
{
	
	int nub = v.GetRow(), porl = 0;
	double k0 = 3, k1 = 5, v_ = 0.0, value = 0.0; vector<double>V_i, V_i1;
	double k0_ = 3, k1_ = 5; double all = 0.0;
	if (nub == 0) return 0;
	double Xigema0 = XIgema0;
	
	MatrixT ppp = H_zhy * P_Trans * H_zhy.Trans();
							   
	MatrixT D_V = (ppp)*SQR(Xigema0);
	
	for (int i = 0; i < nub; i++)
	{
		double a = v(i + 1, 1);
		double b = sqrt(D_V(i + 1, i + 1));
		v_ = fabs(v(i + 1, 1));
		V_i1.push_back(v_);
		v_ /= sqrt(D_V(i + 1, i + 1));
		V_i.push_back(v_);
	}
	double maxValue = *max_element(V_i.begin(), V_i.end());
	auto itMax = max_element(V_i.begin(), V_i.end()) - V_i.begin();
	int imax = (itMax == 0) ? 1 : itMax;
	if (flag_epoch_count < 20) {
		for (int j = 0; j < nub; j++)
		{
			if (fabs(V_i1[j]) > 90)
			{
				double a = 1e8;
				r(j * 2 + 2, j * 2 + 2) *= a;
			}

		}
	}
	else
	{
		if (maxValue <= k1_ && maxValue > k0_)
		{
			double a = 1 / (k0_ / maxValue * SQR((k1_ - maxValue) / (k1_ - k0_)));
			r(imax, imax) *= a;
		}
		else if (maxValue > k1_)
		{
			double a = 1e8;
			r(imax, imax) *= a;
		}

	}
	return 1;
}




int LeastSquarePPP::QC_zhy_CycleSlip(vector<zhydata>& data, PPPObsData& pppdata)
{
	int index = 0; int  KEY = 0;
	double key = 0.0;

	for (int i = 0; i < data[0].m_SatCount; i++)
	{
		for (int j = 0; j < data[1].m_SatCount; j++)
		{
			if (data[0].m_Prn[i] == data[1].m_Prn[j])
			{
				key = ((data[1].m_L1[j] * BDS_WAVELENGTH_B1 - data[1].m_P1[j]) - (data[0].m_L1[i] * BDS_WAVELENGTH_B1 - data[0].m_P1[i]));
				key /= BDS_WAVELENGTH_B1;
				key = fabs(key);
				KEY = int(key);
				if (data[1].m_Prn[j] <= 305)
				{
					if (KEY > 15)  //  
					{
						DeleteSatData(pppdata, j);
						printf("  剔除%d号卫星  ", data[1].m_Prn[j]);
						deletecurrentobsdata(data, j);
					}
				}
				else
				{
					if (KEY > 15)  //  
					{
						DeleteSatData(pppdata, j);

						printf("  剔除%d号卫星  ", data[1].m_Prn[j]);
						deletecurrentobsdata(data, j);
					}
				}
			}
		}
	}

	return 1;
}



int LeastSquarePPP::QC_detect(ppp_option_t& popt, vector<zhydata>& prepppdata,
	PPPObsData& pppdata, PreciseData& predata, GenData& gendata)
{
	int dt = popt.dt;
	ofstream ofclip, ofclipsitandval, mwgfresult;

	
	for (int i = 0; i < pppdata.m_SatCount; i++)
	{
		double DT = fabs(popt.CurrentEpoch - tLast[pppdata.m_Prn[i]]);
		if (DT > dt)
		{
			mw_avg[pppdata.m_Prn[i]] = 0.0;
			mw_avg13[pppdata.m_Prn[i]] = 0.0;
			mw_count[pppdata.m_Prn[i]] = 0;
			mw_count13[pppdata.m_Prn[i]] = 0;
		}
		if (DT >= 2 * dt)
			ION_Initialize(Ion_X[pppdata.m_Prn[i]]);

		ofclip << *(pppdata.m_Prn + i) << setw(4);
		if (i == 0) m_Prn.clear();
		m_Prn.push_back(*(pppdata.m_Prn + i));
	}

	double percentl = 0.00, percentr = 0.10;
	int    frenummwgf = 3;
	int    clip_range = 100;
	int    comm_num = pppdata.m_SatCount * frenummwgf;
	int    pre_num = 0, cur_num = 0;

	vector<unsigned int> prnnametemp, cleartemp;
	vector<int>          recordmiss;
	vector<int>          scclip_sit;
	vector<int>          scclip_value;

	
	if (flag_epoch_count == 0)
	{
		for (int i = 0; i < comm_num; i++)
		{
			double a = (i + 1.0) / comm_num;
			if (a >= percentl && a <= percentr) {
				pre_num = i + 1;
				cur_num = cur_num > pre_num ? cur_num : pre_num;
			}
			if ((i % frenummwgf) == 0) {
				vector<double> temp(frenummwgf, 0.0);
				unsigned int tempprnl = *(pppdata.m_Prn + (i / frenummwgf));
				prnnameL.push_back(tempprnl);
				L_Cur.push_back(temp);
			}
		}
		for (int i = 0; i < (int)L_Cur.size(); i++) {
			ofclip << "|";
			for (int j = 0; j < frenummwgf; j++) ofclip << L_Cur[i][j] << setw(5);
		}
		ofclip << endl;
	}
	else
	{
		for (int i = 0; i < comm_num; i++)
		{
			double a = (i + 1.0) / comm_num;
			if (a >= percentl && a <= percentr) {
				pre_num = i + 1;
				cur_num = cur_num > pre_num ? cur_num : pre_num;
			}
			if ((i % frenummwgf) == 0) {
				unsigned int tempprnl = *(pppdata.m_Prn + (i / frenummwgf));
				prnnametemp.push_back(tempprnl);
				int kmax = (int)prnnameL.size(); int count = 0;
				for (int k = 0; k < kmax; k++) {
					if (tempprnl == prnnameL[k]) break;
					else count++;
					if (count == kmax) {
						vector<double> temp(frenummwgf, 0.0);
						cleartemp.push_back(tempprnl);
						prnnameL.push_back(tempprnl);
						L_Cur.push_back(temp);
					}
				}
			}
		}
		int kmax = (int)prnnameL.size();
		for (int k = 0; k < kmax; k++) {
			int count = 0;
			for (int i = 0; i < (int)prnnametemp.size(); i++) {
				if (prnnameL[k] == prnnametemp[i]) break;
				else count++;
				if (count == (int)prnnametemp.size())
					recordmiss.push_back(prnnameL[k]);
			}
		}
	}

	if (flag_epoch_count >= 2)
	{
		vector<unsigned int>   prnnameLtemp;
		vector<vector<double>> L_Curtemp;
		for (int i = 0; i < (int)prnnametemp.size(); i++) {
			for (int j = 0; j < (int)prnnameL.size(); j++) {
				if (prnnameL[j] == prnnametemp[i]) {
					prnnameLtemp.push_back(prnnameL[j]);
					L_Curtemp.push_back(L_Cur[j]);
					break;
				}
			}
		}
		prnnameL.clear(); L_Cur.clear();
		prnnameL = prnnameLtemp; L_Cur = L_Curtemp;

		ofclipsitandval << pppdata.m_SatCount << '\t' << "|";
		for (int i = 0; i < (int)scclip_sit.size(); i++)
			ofclipsitandval << scclip_sit[i] << '!' << scclip_value[i] << '\t';
		ofclipsitandval << endl;
		for (int i = 0; i < (int)L_Cur.size(); i++) {
			ofclip << "|";
			for (int j = 0; j < frenummwgf; j++) ofclip << L_Cur[i][j] << setw(5);
		}
		ofclip << endl;
	}

	
	vector<vector<double>> mwNvalue, gfNvalue;
	vector<double> test1, test2;
	mwNvalue.push_back(test1); mwNvalue.push_back(test2);
	gfNvalue.push_back(test1); gfNvalue.push_back(test2);

	detslp_MW(popt, pppdata, predata, gendata, L_Cur, mwNvalue, frenummwgf);
	detslp_GF(popt, pppdata, predata, gendata, L_Cur, gfNvalue, frenummwgf);

	
	vector<double> clipsitandvalue(comm_num, 0.0);
	vector<double> clipsitandvalue12(comm_num, 0.0);
	vector<double> clipsitandvalue13(comm_num, 0.0);


	const double slipwavemg[3] = { 0.1920394, 0.2363325, 0.2483494 };

	if (!_Cycle_Sliptest.empty() && !_Cycle_Sliptestgf.empty())
	{
		for (int i = 0; i < (int)_Cycle_Sliptest.size(); i++)
			for (int j = 0; j < (int)_Cycle_Sliptestgf.size(); j++)
				if (_Cycle_Sliptest[i] == _Cycle_Sliptestgf[j])
				{
					MatrixT B, L, V, X, Nbb, P;
					int settemp = 0;
					B.Resize(2, 2); L.Resize(2, 1); X.Resize(2, 1);
					V.Resize(2, 1); Nbb.Resize(2, 2); P.Resize(2, 2);
					B(1, 1) = 1; B(1, 2) = -1;
					B(2, 1) = 1; B(2, 2) = -(slipwavemg[1] / slipwavemg[0]);
					L(1, 1) = mwNvalue[0][i]; L(2, 1) = gfNvalue[0][j];
					P(1, 1) = 1.0; P(2, 2) = 1.0;
					if (LeastSquare(B, L, P, V, Nbb, X) == -1);

					for (int k = 0; k < pppdata.m_SatCount; k++)
						if (_Cycle_Sliptest[i] == pppdata.m_Prn[k]) {
							for (int m = 0; m < (int)_Cycle_Sliptest13.size(); m++)
								if (_Cycle_Sliptest[i] == _Cycle_Sliptest13[m]) { settemp++; break; }
							for (int m = 0; m < (int)_Cycle_Sliptestgf13.size(); m++)
								if (_Cycle_Sliptest[i] == _Cycle_Sliptestgf13[m]) { settemp++; break; }
							if (settemp == 0) {
								clipsitandvalue12[k * frenummwgf + 0] = 0;
								clipsitandvalue12[k * frenummwgf + 1] = -mwNvalue[0][i];
								continue;
							}
							clipsitandvalue12[k * frenummwgf + 0] = X(1, 1);
							clipsitandvalue12[k * frenummwgf + 1] = X(2, 1);
						}
					break;
				}
		for (auto& v : clipsitandvalue12) v = round(v);
	}

	// L1-L3
	if (!_Cycle_Sliptest13.empty() && !_Cycle_Sliptestgf13.empty())
	{
		for (int i = 0; i < (int)_Cycle_Sliptest13.size(); i++)
			for (int j = 0; j < (int)_Cycle_Sliptestgf13.size(); j++)
				if (_Cycle_Sliptest13[i] == _Cycle_Sliptestgf13[j])
				{
					MatrixT B, L, V, X, Nbb, P;
					int settemp = 0;
					B.Resize(2, 2); L.Resize(2, 1); X.Resize(2, 1);
					V.Resize(2, 1); Nbb.Resize(2, 2); P.Resize(2, 2);
					B(1, 1) = 1; B(1, 2) = -1;
					B(2, 1) = 1; B(2, 2) = -(slipwavemg[2] / slipwavemg[0]);
					L(1, 1) = mwNvalue[1][i]; L(2, 1) = gfNvalue[1][j];
					P(1, 1) = 1.0; P(2, 2) = 1.0;
					if (LeastSquare(B, L, P, V, Nbb, X) == -1);

					for (int k = 0; k < pppdata.m_SatCount; k++)
						if (_Cycle_Sliptest13[i] == pppdata.m_Prn[k]) {
							for (int m = 0; m < (int)_Cycle_Sliptest.size(); m++)
								if (_Cycle_Sliptest13[i] == _Cycle_Sliptest[m]) { settemp++; break; }
							for (int m = 0; m < (int)_Cycle_Sliptestgf.size(); m++)
								if (_Cycle_Sliptest13[i] == _Cycle_Sliptestgf[m]) { settemp++; break; }
							if (settemp == 0) {
								clipsitandvalue13[k * frenummwgf + 0] = 0;
								clipsitandvalue13[k * frenummwgf + 2] = -mwNvalue[1][i];
								continue;
							}
							clipsitandvalue13[k * frenummwgf + 0] = X(1, 1);
							clipsitandvalue13[k * frenummwgf + 2] = X(2, 1);
						}
					break;
				}
		for (auto& v : clipsitandvalue13) v = round(v);
	}

	if (flag_epoch_count >= 2)
	{
		double recordrepair = 0.0;
		for (int i = 0; i < (int)scclip_sit.size(); i++) {
			if (fabs(clipsitandvalue12[scclip_sit[i]] - scclip_value[i]) < 0.5) { recordrepair += 1.0; continue; }
			if (fabs(clipsitandvalue13[scclip_sit[i]] - scclip_value[i]) < 0.5)   recordrepair += 1.0;
		}
		if (!scclip_sit.empty()) recordrepair /= scclip_sit.size();
		mwgfresult << recordrepair << endl;
	}

	// ──────────────────────────────────────────────────────────
	// 方法1：TurboEdit（MW+GF）CSV 输出
	// ──────────────────────────────────────────────────────────
	{
		// [FIX-LP-3] 使用 trunc 模式只在第一次运行时清空文件
		static bool te_out_header = false;
		auto te_mode = te_out_header ? std::ios::app : (std::ios::out | std::ios::trunc);
		std::ofstream of_te("D:\\data10.0\\outt\\TurboEdit_detect_repair.csv", te_mode);
		if (!te_out_header) {
			of_te << "历元,PRN,系统,频点组合,"
				<< "ΔMW差值(周),ΔGF差值(m),"
				<< "ΔN1估计(周),ΔN2估计(周),"
				<< "ΔN1整周,ΔN2整周,"
				<< "修复量L1(m),修复量L2/L3(m),"
				<< "是否模拟周跳,模拟真值(周),周跳判定\n";
			te_out_header = true;
		}

		const double slipw[3] = { 0.1920394, 0.2363325, 0.2483494 };

		for (int i = 0; i < (int)_Cycle_Sliptest.size(); i++) {
			unsigned int sat = _Cycle_Sliptest[i];
			char sys = GetSystem(sat);
			int k = -1;
			for (int j = 0; j < pppdata.m_SatCount; j++)
				if (pppdata.m_Prn[j] == sat) { k = j; break; }
			if (k < 0) continue;

			double dN1 = clipsitandvalue12[k * frenummwgf + 0];
			double dN2 = clipsitandvalue12[k * frenummwgf + 1];
			double dmw = (i < (int)mwNvalue[0].size()) ? mwNvalue[0][i] : 0.0;
			double dgf = 0.0;
			for (int jj = 0; jj < (int)_Cycle_Sliptestgf.size(); jj++)
				if (_Cycle_Sliptestgf[jj] == sat && jj < (int)gfNvalue[0].size())
				{
					dgf = gfNvalue[0][jj]; break;
				}

			double rep1_m = dN1 * slipw[0], rep2_m = dN2 * slipw[1];
			std::string verdict;
			if (dN1 != 0 && dN2 != 0) verdict = "slip_L1andL2";
			else if (dN1 != 0)    verdict = "slip_L1only";
			else if (dN2 != 0)    verdict = "slip_L2only";
			else                verdict = "no_slip";
			if (verdict == "no_slip" && (fabs(dmw) > 0.5 || fabs(dgf) > 0.03))
				verdict = "subint_or_iono";

			of_te << flag_epoch_count << "," << sat << "," << sys << "," << "L1-L2,"
				<< std::setprecision(6) << dmw << "," << dgf << ","
				<< dN1 << "," << dN2 << "," << dN1 << "," << dN2 << ","
				<< rep1_m << "," << rep2_m << "," << "0,0," << verdict << "\n";
		}

		for (int i = 0; i < (int)_Cycle_Sliptest13.size(); i++) {
			unsigned int sat = _Cycle_Sliptest13[i];
			char sys = GetSystem(sat);
			int k = -1;
			for (int j = 0; j < pppdata.m_SatCount; j++)
				if (pppdata.m_Prn[j] == sat) { k = j; break; }
			if (k < 0) continue;

			double dN1 = clipsitandvalue13[k * frenummwgf + 0];
			double dN3 = clipsitandvalue13[k * frenummwgf + 2];
			double dmw13 = (i < (int)mwNvalue[1].size()) ? mwNvalue[1][i] : 0.0;
			double dgf13 = 0.0;
			for (int jj = 0; jj < (int)_Cycle_Sliptestgf13.size(); jj++)
				if (_Cycle_Sliptestgf13[jj] == sat && jj < (int)gfNvalue[1].size())
				{
					dgf13 = gfNvalue[1][jj]; break;
				}

			double rep1_m = dN1 * slipw[0], rep3_m = dN3 * slipw[2];
			std::string verdict;
			if (dN1 != 0 && dN3 != 0) verdict = "slip_L1andL3";
			else if (dN1 != 0)    verdict = "slip_L1only";
			else if (dN3 != 0)    verdict = "slip_L3only";
			else                verdict = "no_slip";
			if (verdict == "no_slip" && (fabs(dmw13) > 0.5 || fabs(dgf13) > 0.03))
				verdict = "subint_or_iono";

			of_te << flag_epoch_count << "," << sat << "," << sys << "," << "L1-L3,"
				<< std::setprecision(6) << dmw13 << "," << dgf13 << ","
				<< dN1 << "," << dN3 << "," << dN1 << "," << dN3 << ","
				<< rep1_m << "," << rep3_m << "," << "0,0," << verdict << "\n";
		}

		of_te.flush(); of_te.close();

		{
			static int  te_total_all = 0, te_detect_all = 0;
			static bool sum_hdr = false;
			int detected = (int)_Cycle_Sliptest.size() + (int)_Cycle_Sliptest13.size();
			int tested = pppdata.m_SatCount;
			te_total_all += tested; te_detect_all += detected;
			auto sm = (sum_hdr ? std::ios::app : (std::ios::out | std::ios::trunc));
			std::ofstream of_sum("D:\\data10.0\\outt\\TurboEdit_summary.csv", sm);
			if (!sum_hdr) {
				of_sum << "历元,本历元卫星数,本历元探测周跳数,累计卫星数,累计周跳数,累计探测率(%)\n";
				sum_hdr = true;
			}
			double rate = te_total_all > 0 ? (double)te_detect_all / te_total_all * 100.0 : 0.0;
			of_sum << flag_epoch_count << "," << tested << "," << detected << ","
				<< te_total_all << "," << te_detect_all << ","
				<< std::fixed << std::setprecision(2) << rate << "\n";
			of_sum.close();
		}
	}


	

	flag_epoch_count++;
	return 0;

	


	//detslp_GF13(popt, pppdata, predata, gendata);
#if (1)
	for (int i = 0; i < pppdata.m_SatCount; i++)
	{
		calculate_dion(pppdata.m_Prn[i]);//军哥

		ION_del_one_first(Ion_X[pppdata.m_Prn[i]]);
	}
#endif
	return 0;
}

int LeastSquarePPP::detslp_MW(ppp_option_t& popt, PPPObsData& pppdata, PreciseData& predata, GenData& gendata, vector<vector<double>> L_Cur, vector<vector<double>>& mwNvalue, int frenummwgf)
{
	double L1, L2, P1, P2; int g_or_c = 2, orbn = 0; bool IS = false; bool IS13 = false;
	double SatPCV[NFREQ] = { 0.0 }, RecPCV[NFREQ] = { 0.0 }; double satpos[3] = { 0.0 }; double PhaseWindup = 0.0;
	int satN = pppdata.m_SatCount; int bLowElev = 0; double el, thres; double rec_ant_delta[3] = { Antdiff[1], Antdiff[0], Antdiff[2] };
	double wavelength1 = 0, wavelength2 = 0, wavelength3 = 0, wavelength4 = 0, wavelength5 = 0;
	double fact = 1.0;

	for (int i = 0; i < satN; i++)
	{
		IS = false; bool IS13 = false;
		unsigned int satno = pppdata.m_Prn[i];
		if (GetSystem(satno) == 'G') {
			orbn = 0;
			wavelength1 = Get_WaveLength(satno, 1, orbn);//L1
			wavelength2 = Get_WaveLength(satno, 2, orbn);//L2
			wavelength3 = Get_WaveLength(satno, 5, orbn);//L5
		}
		if (GetSystem(satno) == 'E')
		{
			orbn = 0;
			wavelength1 = Get_WaveLength(satno, 1, orbn);//E1
			wavelength2 = Get_WaveLength(satno, 2, orbn);//E5A
			wavelength3 = Get_WaveLength(satno, 3, orbn);//E5b
			wavelength4 = Get_WaveLength(satno, 4, orbn);//E5(A+B)
			wavelength5 = Get_WaveLength(satno, 5, orbn);//E6
		}
		if (GetSystem(satno) == 'C')
		{
			orbn = 0;
			wavelength1 = Get_WaveLength(satno, 1, orbn);//B1
			wavelength2 = Get_WaveLength(satno, 3, orbn);//B3
			wavelength3 = Get_WaveLength(satno, 2, orbn);//B2
			wavelength4 = Get_WaveLength(satno, 4, orbn);//B1C
			wavelength5 = Get_WaveLength(satno, 5, orbn);//B2a

		}
		if (GetSystem(satno) == 'R')
		{
			orbn = Get_Nav_Orbitn(satno, predata.m_NavData);
			wavelength1 = Get_WaveLength(satno, 1, orbn);//G1
			wavelength2 = Get_WaveLength(satno, 2, orbn);//G2
			wavelength3 = Get_WaveLength(satno, 3, orbn);//G1a
			wavelength4 = Get_WaveLength(satno, 4, orbn);//G2a
			wavelength5 = Get_WaveLength(satno, 5, orbn);//G3
		}

		double azel[2] = { m_SatInfo.m_Azi[i], m_SatInfo.m_Ele[i] };
		satpos[0] = m_SatInfo.m_Xs[3 * i + 0];
		satpos[1] = m_SatInfo.m_Xs[3 * i + 1];
		satpos[2] = m_SatInfo.m_Xs[3 * i + 2];
		SatAntPCV(GetSystem(satno), satno, satpos, m_Sta_tmp, gendata.m_pcvs, SatPCV);
		RecAntModel(GetSystem(satno), &gendata.m_pcvr, rec_ant_delta, azel, 1, RecPCV);
		WindupCorr(popt.SunPos, satpos, m_Sta_tmp, PhaseWindup);

		double P1 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr1[i]);
		double P2 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr2[i]);
		double P3 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr3[i]);
		double P4 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr4[i]);
		double P5 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr5[i]);
		double L1 = CarrierPhase_PCV_Correct(wavelength1, RecPCV[0], SatPCV[0], pppdata.m_L1[i]) + L_Cur[i][0];
		double L2 = CarrierPhase_PCV_Correct(wavelength2, RecPCV[0], SatPCV[0], pppdata.m_L2[i]) + L_Cur[i][1];
		double L3 = 0.0;
		if (frenummwgf == 3) { L3 = CarrierPhase_PCV_Correct(wavelength3, RecPCV[0], SatPCV[0], pppdata.m_L3[i]) + L_Cur[i][2]; }
		else { L3 = CarrierPhase_PCV_Correct(wavelength3, RecPCV[0], SatPCV[0], pppdata.m_L3[i]); }
		double L4 = CarrierPhase_PCV_Correct(wavelength4, RecPCV[0], SatPCV[0], pppdata.m_L4[i]);
		double L5 = CarrierPhase_PCV_Correct(wavelength5, RecPCV[0], SatPCV[0], pppdata.m_L5[i]);
		L1 -= wavelength1 * PhaseWindup;
		L2 -= wavelength2 * PhaseWindup;
		L3 -= wavelength3 * PhaseWindup;

		mw_cur[satno] = (L1 - L2) - (wavelength2 - wavelength1) / (wavelength1 + wavelength2) * (P1 / wavelength1 + P2 / wavelength2);
		if ((L3 != 0 && P3 != 0) && frenummwgf == 3)
		{
			mw_cur13[satno] = (L1 - L3) - (wavelength3 - wavelength1) / (wavelength1 + wavelength3) * (P1 / wavelength1 + P3 / wavelength3);
		}
		if (flag_epoch_count == 0)continue;

		if (mw_cur[satno] == 0 || mw_avg[satno] == 0) continue;
		if ((L3 != 0 && P3 != 0) && frenummwgf == 3)
		{
			if (mw_cur13[satno] == 0 || mw_avg13[satno] == 0) continue;
		}

		double key = fabs(mw_cur[satno] - mw_avg[satno]);
		double key13 = 0.0;
		if ((L3 != 0 && P3 != 0) && frenummwgf == 3) { key13 = fabs(mw_cur13[satno] - mw_avg13[satno]); }

		double elev = m_SatInfo.m_Ele[i] * R2D;
		el = elev;
		if (elev < popt.ele_mask) {
			el = popt.ele_mask;
			bLowElev = 1;
		}
		if (bLowElev) continue;
		double dtmp = el * D2R;
		if (el >= 20.0) thres = ThresMW;
		else thres = -ThresMW * 0.1 * dtmp + 3 * ThresMW;

		if (key > min(thres, 6.0)) IS = true;
		if (!IS) {
			if (key >= 1.0) {
				
				_Cycle_Sliptest.push_back(satno);
				mwNvalue[0].push_back((mw_cur[satno]) - mw_pre[satno]);
			
			}
		}
		else
		{
			
			_Cycle_Sliptest.push_back(satno);
			mwNvalue[0].push_back((mw_cur[satno]) - mw_pre[satno]);
		}
		if (frenummwgf == 3)
		{
			if (key13 > min(thres, 6.0)) IS13 = true;
			if (!IS13) {
				if (key13 >= 1.0) {
					
					_Cycle_Sliptest13.push_back(satno);
					mwNvalue[1].push_back((mw_cur13[satno]) - mw_pre13[satno]);
				}
			}
			else
			{
				_Cycle_Sliptest13.push_back(satno);
				mwNvalue[1].push_back((mw_cur13[satno]) - mw_pre13[satno]);
			}
		}
	}
	return 0;
}

int LeastSquarePPP::detslp_GF(ppp_option_t& popt, PPPObsData& pppdata, PreciseData& predata, GenData& gendata, vector<vector<double>> L_Cur, vector<vector<double>>& gfNvalue, int frenummwgf)
{
	//ofstream gf01,gf02;
//gf01.open("C:\\Users\\15149\\Desktop\\shiy\\gf12.out", ios::app);
//gf02.open("C:\\Users\\15149\\Desktop\\shiy\\gf13.out", ios::app);
	int satN = pppdata.m_SatCount; int g_or_c = 2, orbn = 0;
	double SatPCV[NFREQ] = { 0.0 }, RecPCV[NFREQ] = { 0.0 }; double satpos[3] = { 0.0 }; double PhaseWindup = 0.0;
	double rec_ant_delta[3] = { Antdiff[1], Antdiff[0], Antdiff[2] };
	double thres = 0.0;
	for (int i = 0; i < satN; i++)
	{
		double elev = m_SatInfo.m_Ele[i]; unsigned int satno = pppdata.m_Prn[i];
		if (GetSystem(satno) == 'G' || GetSystem(satno) == 'E') { g_or_c = 2; orbn = 0; }
		if (GetSystem(satno) == 'C') { g_or_c = 3; orbn = 0; }
		if (GetSystem(satno) == 'R') { orbn = Get_Nav_Orbitn(satno, predata.m_NavData); g_or_c = 2; }
		double wavelength1 = Get_WaveLength(satno, 1, orbn);  
		double wavelength2 = Get_WaveLength(satno, g_or_c, orbn);
		double wavelength3 = Get_WaveLength(satno, 5, orbn);
		double azel[2] = { m_SatInfo.m_Azi[i], m_SatInfo.m_Ele[i] };
		satpos[0] = m_SatInfo.m_Xs[3 * i + 0];
		satpos[1] = m_SatInfo.m_Xs[3 * i + 1];
		satpos[2] = m_SatInfo.m_Xs[3 * i + 2];
		SatAntPCV(GetSystem(satno), satno, satpos, m_Sta_tmp, gendata.m_pcvs, SatPCV);
		RecAntModel(GetSystem(satno), &gendata.m_pcvr, rec_ant_delta, azel, 1, RecPCV);
		WindupCorr(popt.SunPos, satpos, m_Sta_tmp, PhaseWindup);

		double P1 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr1[i]);
		double P2 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr2[i]);
		double P3 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr3[i]);
		double L1 = CarrierPhase_PCV_Correct(wavelength1, RecPCV[0], SatPCV[0], pppdata.m_L1[i]) + L_Cur[i][0];
		double L2 = CarrierPhase_PCV_Correct(wavelength2, RecPCV[0], SatPCV[0], pppdata.m_L2[i]) + L_Cur[i][1];
		double L3 = 0.0;
		if (frenummwgf == 3) { L3 = CarrierPhase_PCV_Correct(wavelength3, RecPCV[0], SatPCV[0], pppdata.m_L3[i]) + L_Cur[i][2]; }
		else { L3 = CarrierPhase_PCV_Correct(wavelength3, RecPCV[0], SatPCV[0], pppdata.m_L3[i]); }
		L1 -= wavelength1 * PhaseWindup;
		L2 -= wavelength2 * PhaseWindup;
		L3 -= wavelength3 * PhaseWindup;
		double fre2 = wavelength2 / LIGHTSPEED;
		double fre1 = wavelength1 / LIGHTSPEED;
		double fre3 = wavelength3 / LIGHTSPEED;
		int listhead[8] = { 345,320,336,329,330,332,319,335 };
		int satnum = 0; double gf13set = 0.0;
		for (int j = 0; j < 8; j++)
		{
			if (listhead[j] == satno)
			{
				satnum = j;
				break;
			}
		}
		if ((flag_epoch_count >= 2) && frenummwgf == 3)
		{
			gf13set = setvaluecal(flag_epoch_count, satnum);
		}

		gf_cur[satno] = wavelength1 * L1 - wavelength2 * L2;
		gf_cur[satno] /= wavelength1;
		if (frenummwgf == 3)
		{
			gf_cur13[satno] = wavelength1 * L1 - wavelength3 * L3;
			gf_cur13[satno] /= wavelength1;
		}


		if (flag_epoch_count == 0 || gf_pre[satno] == 0)
		{
			continue;
		}
		if (_Cycle_Sliptest.size() > 1)
		{
			for (int cy = 0; cy < _Cycle_Sliptest.size(); cy++)
			{
				if (_Cycle_Sliptest[cy] == satno)
				{
					gf_cur[satno] == 0;
					continue;
				}
			}
		}
		if (frenummwgf == 3)
		{
			if (_Cycle_Sliptest13.size() > 1)
			{
				for (int cy = 0; cy < _Cycle_Sliptest13.size(); cy++)
				{
					if (_Cycle_Sliptest13[cy] == satno)
					{
						gf_cur13[satno] == 0;
						continue;
					}
				}
			}
		}

		if (elev < popt.ele_mask * R2D) elev = popt.ele_mask * R2D;
		if (elev >= 15.0)thres = THresGF;
		else thres = -THresGF / 15.0 * elev + 2 * THresGF;

		if (0)
		{
			if (Ion_X[satno].a == 0)
			{
				if (fabs(gf_cur[satno] - gf_pre[satno]) > min(thres, 1.5))
				{
					gf_ion[satno] = 0.001; gf_sigemaIon[satno] = 0.0001;
					_Cycle_Sliptestgf.push_back(satno);
					gfNvalue[0].push_back((gf_cur[satno]) - gf_pre[satno]);
				}
				else
				{
					int cy_ = 1;
					if (_Cycle_Sliptest.size() > 1)
					{
						for (int cy = 0; cy < _Cycle_Sliptest.size(); cy++)
						{
							if (_Cycle_Sliptest[cy] == satno) { cy_ = 0; }
						}
					}
					if (cy_ == 1)
					{
						Ion_X[satno].gf_value.push_back(gf_cur[satno] / (SQR(wavelength2) / SQR(wavelength1) - 1));
						Ion_X[satno].dgf_value.push_back((gf_cur[satno] - gf_pre[satno]) / (SQR(wavelength2) / SQR(wavelength1) - 1));
						Ion_X[satno].time.push_back(flag_epoch_count);
						//SatInfo
						for (int elesize = 0; elesize < m_SatInfo.m_Size; elesize++)
						{
							if (m_SatInfo.m_Prn[elesize] == satno)
							{
								Ion_X[satno].ele.push_back(m_SatInfo.m_Ele[elesize]);
							}
						}


						//Ion_X[satno].time1.push_back(flag_epoch_count);
						Ion_X[satno].num++;
						//Ion_X[satno].num1++;
					}
					double a = 1 / (SQR(wavelength2) / SQR(wavelength1) - 1);
					double b = SQR(fre1) / (SQR(fre2) - SQR(fre1));

					gf_ion[satno] = (gf_cur[satno] - gf_pre[satno]) / (SQR(wavelength2) / SQR(wavelength1) - 1);
					gf_sigemaIon[satno] = 0.006 / (SQR(wavelength2) / SQR(wavelength1) - 1);
				}
				if (fabs(gf_cur13[satno] - gf_pre13[satno]) > min(thres, 1.5))
				{
					_Cycle_Sliptestgf13.push_back(satno);
					gfNvalue[1].push_back((gf_cur13[satno]) - gf_pre13[satno]);
				}
			}
			else
			{

				double ion_predict = fabs(Ion_X[satno].a * (flag_epoch_count)+Ion_X[satno].b * (flag_epoch_count * flag_epoch_count) + Ion_X[satno].c * SQR(flag_epoch_count) * (flag_epoch_count));
				double ion_gf = fabs(gf_cur[satno] - gf_pre[satno]);
				/*if (satno==320)
				{
					MatrixT a(1, 3);
					a(1, 1) = flag_epoch_count;
					a(1, 2) = gf_cur[satno];
					a(1, 3) = ion_predict;
					a.Output_to_File_no("E:\\320.txt", 1);
				}*/
				//if ((ion_gf-ion_predict)>min(thres, 1.5)|| fabs(gf_cur[satno] - gf_pre[satno]) > min(thres, 1.5))
				if (fabs(gf_cur[satno] - gf_pre[satno]) > min(thres, 1.5))
				{
					_Cycle_Sliptestgf.push_back(satno);
					gfNvalue[0].push_back((gf_cur[satno]) - gf_pre[satno]);
				}
				else
				{
					int cy_ = 1;
					if (_Cycle_Sliptest.size() > 1)
					{
						for (int cy = 0; cy < _Cycle_Sliptest.size(); cy++)
						{
							if (_Cycle_Sliptest[cy] == satno) { cy_ = 0; }
						}
					}
					if (cy_ == 1)
					{
						Ion_X[satno].gf_value.push_back(gf_cur[satno] / (SQR(wavelength2) / SQR(wavelength1) - 1));
						Ion_X[satno].dgf_value.push_back((gf_cur[satno] - gf_pre[satno]) / (SQR(wavelength2) / SQR(wavelength1) - 1));
						Ion_X[satno].time.push_back(flag_epoch_count);
						//Ion_X[satno].time1.push_back(flag_epoch_count);


						for (int elesize = 0; elesize < m_SatInfo.m_Size; elesize++)
						{
							if (m_SatInfo.m_Prn[elesize] == satno)
							{
								Ion_X[satno].ele.push_back(m_SatInfo.m_Ele[elesize]);
							}
						}

						Ion_X[satno].num++;
						//Ion_X[satno].num1++;
					}
					gf_ion[satno] = (gf_cur[satno] - gf_pre[satno]) / (SQR(wavelength2) / SQR(wavelength1) - 1);
					gf_sigemaIon[satno] = 0.006 / (SQR(wavelength2) / SQR(wavelength1) - 1);
				}
				if (fabs(gf_cur13[satno] - gf_pre13[satno]) > min(thres, 1.5))
				{
					_Cycle_Sliptestgf13.push_back(satno);
					gfNvalue[1].push_back((gf_cur13[satno]) - gf_pre13[satno]);
				}
			}
		}

		if (fabs(gf_cur[satno] - gf_pre[satno]) > 0.03)
		{
			_Cycle_Sliptestgf.push_back(satno);
			gfNvalue[0].push_back((gf_cur[satno]) - gf_pre[satno]);
		}
		if (frenummwgf == 3)
		{
			if (fabs(gf13set - (gf_cur13[satno] - gf_pre13[satno])) > 1.0)
			{
				_Cycle_Sliptestgf13.push_back(satno);
				gfNvalue[1].push_back((gf_cur13[satno] - gf_pre13[satno]) - gf13set);
			}
		}
	}

	return 0;
}

double LeastSquarePPP::setvaluecal(int eopchtime, int satnum)
{
	if (satnum == 0)//345
	{
		eopchtime -= 1;
		return 9e-07 * (eopchtime) * (eopchtime)-0.0138 * (eopchtime)-9.2427;
	}
	else if (satnum == 1)//320
	{
		eopchtime -= 1;
		return -7e-07 * (eopchtime) * (eopchtime)-0.0099 * (eopchtime)+46.881;
	}
	else if (satnum == 2)//336
	{
		eopchtime -= 1;
		return 1e-06 * (eopchtime) * (eopchtime)-0.0046 * (eopchtime)-69.742;
	}
	else if (satnum == 3)//329
	{
		eopchtime -= 1;
		return -4e-11 * (eopchtime) * (eopchtime) * (eopchtime)+9e-08 * (eopchtime) * (eopchtime)-0.0009 * (eopchtime)+9.7808;
	}
	else if (satnum == 4)//330
	{
		eopchtime -= 1;
		return -2e-11 * (eopchtime) * (eopchtime) * (eopchtime)-2e-07 * (eopchtime) * (eopchtime)-0.0003 * (eopchtime)-31.125;
	}
	else if (satnum == 5)//332
	{
		eopchtime -= 1;
		return 3e-07 * (eopchtime) * (eopchtime)-0.0131 * (eopchtime)-21.604;
	}
	else if (satnum == 6)//319
	{
		eopchtime -= 466;
		return -1e-06 * (eopchtime) * (eopchtime)+0.0003 * (eopchtime)+79.278;
	}
	else if (satnum == 7)//335
	{
		eopchtime -= 1394;
		return -3e-11 * (eopchtime) * (eopchtime) * (eopchtime)+2e-07 * (eopchtime) * (eopchtime)+0.0005 * (eopchtime)+40.494;
	}
	else
	{
		return 0.0;
	}
}

void LeastSquarePPP::calculate_dion(unsigned int satno)
{
	if (Ion_X[satno].num < 10)return;
	MatrixT B_i(Ion_X[satno].num, 3);
	MatrixT L_i(Ion_X[satno].num, 1);
	MatrixT P_i; P_i.Eye(Ion_X[satno].num);
	MatrixT V_i(Ion_X[satno].num, 1);
	MatrixT nn(3, 3);
	MatrixT x_i(3, 1);


	MatrixT B_i1(Ion_X[satno].num, 3);
	MatrixT L_i1(Ion_X[satno].num, 1);
	MatrixT P_i1; P_i1.Eye(Ion_X[satno].num);
	MatrixT V_i1(Ion_X[satno].num, 1);
	MatrixT nn1(3, 3);
	MatrixT x_i1(3, 1);
	//////////////////////////////////////////////////////////
	//*******************************************************
	int numb = 0;
	//*******************************************************
	// //////////////////////////////////////////////////////////
	//fprintf(fp, "%d,%15.9f \n", flag_epoch_count, sigema0__1);

	double s_30 = 0, s_30_ = 0;
	double s_60 = 0, s_60_ = 0; int  k = 0, k1 = 0;
	double sum = 0;
	vector <double>edi_5, edi_5_;
	int interval = 2;
	double interval1 = 0;
	double sum_30 = 0;
	//for (int interval = 2; interval <= 90; interval++)
	//{
#if (0)
	for (int i = 0; i < Ion_X[satno].num; i++)
	{
		B_i(i + 1, 1) = 1; B_i(i + 1, 2) = Ion_X[satno].time[i]; B_i(i + 1, 3) = SQR(Ion_X[satno].time[i]);

		L_i(i + 1, 1) = Ion_X[satno].dgf_value[i];
	}

	for (int i = 0; i < Ion_X[satno].num - interval; i++)
	{
		sum_30 += Ion_X[satno].dgf_value[i + interval];
	}
	double mean_30 = sum_30 / (Ion_X[satno].num - interval);
	for (int i = 0; i < Ion_X[satno].num - interval; i++)
	{
		s_30 += SQR(Ion_X[satno].dgf_value[i + interval]);
		s_30_ += SQR(Ion_X[satno].dgf_value[i + interval] - mean_30);
		k1++;
	}
	for (int i = 0; i < Ion_X[satno].num; i++)
	{
		if (i >= interval)
		{
			sum += (Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval]);
		}
	}
	double mean_60 = sum / (Ion_X[satno].num - interval);
	for (int i = 0; i < Ion_X[satno].num; i++)
	{
		if (i == 0) continue;
		double b = (i >= interval) ? SQR(Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval] - mean_60) : 0;
		double a = (i >= interval) ? SQR(Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval]) : 0;
		if (i >= interval)
		{
			//k+= Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval];
			edi_5_.push_back(Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval] - mean_60);
			edi_5.push_back(Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval]);
		}
		s_60 += a;
		s_60_ += b;
	}


#endif

#if(0)
	for (int i = 90; i < 390; i++)
	{
		sum_30 += Ion_X[satno].dgf_value[i];
	}
	double mean_30 = sum_30 / 300;
	for (int i = 90; i < 390; i++)
	{
		s_30 += SQR(Ion_X[satno].dgf_value[i]);
		s_30_ += SQR(Ion_X[satno].dgf_value[i] - mean_30);
		k1++;
	}
	for (int i = 90; i < 390; i++)
	{

		sum += (Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval]);

	}
	double mean_60 = sum / (300);
	for (int i = 90; i < 390; i++)
	{


		if (i == 0) continue;
		double b = SQR(Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval] - mean_60);
		double a = SQR(Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval]);

		//k+= Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval];
		edi_5_.push_back(Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval] - mean_60);
		edi_5.push_back(Ion_X[satno].gf_value[i] - Ion_X[satno].gf_value[i - interval]);

		s_60 += a;
		s_60_ += b;

	}
#endif
#if (0)
	for (int i = 0; i < Ion_X[satno].num; i++)
	{
		B_i(i + 1, 1) = 1; B_i(i + 1, 2) = Ion_X[satno].time[i]; B_i(i + 1, 3) = SQR(Ion_X[satno].time[i]);

		L_i(i + 1, 1) = Ion_X[satno].dgf_value[i];

		s_30 += SQR(Ion_X[satno].dgf_value[i]);
		if (i == 0) continue;
		s_60 += (i % 2 == 0) ? (SQR(Ion_X[satno].dgf_value[i] - Ion_X[satno].dgf_value[i - 2])) : 0;

	}
#endif
	sigema30[satno] = s_30 / (Ion_X[satno].num - interval);
	sigema60[satno] = s_60 / (Ion_X[satno].num - interval);
	//sigema60[satno] = s_60 / (Ion_X[satno].num * 0.5); 
	//double apsd60 = sqrt(s_60 / (Ion_X[satno].num - interval ) / interval / interval)*1000;//看好除那个   1s 跟30s 不一样
	//double apsd30 = sqrt(s_30 / (Ion_X[satno].num - interval)/1) * 1000;
	//正常
	double psd30_ = (s_30_ / (Ion_X[satno].num - interval));
	double psd60_ = (s_60_ / (Ion_X[satno].num - interval));

	double apsd30_ = sqrt(s_30_ / (Ion_X[satno].num - interval) / 30);
	//double apsd60_ = sqrt(s_60_ / 300 / interval) * 1000;

#if(0)
	if (flag_epoch_count < 120 * 2)
	{
		PSD[satno] = sqrt(fabs(sigema60[satno] - sigema30[satno]) / 30);
	}
	else
	{
		PSD[satno] = sqrt(fabs(sigema60[satno] - sigema30[satno])) / 30;
	}
#endif
	PSD[satno] = sqrt(fabs(sigema60[satno] - sigema30[satno]) / (interval - 1));
	//double psd1 = sqrt(fabs(s_60_ / (Ion_X[satno].num - interval) - (s_30_ / (Ion_X[satno].num - interval))) / (interval - 1));

	///double psd1 = sqrt(fabs((s_60_ / 300) - (s_30_ /300)) / (interval - 1));
	double psd1 = sqrt(fabs((psd60_)-(psd30_)) / (30));  //正常数据

	PSD[satno] = apsd30_;

	for (int i = 0; i < Ion_X[satno].num; i++)
	{

		B_i1(i + 1, 1) = Ion_X[satno].time[i]; B_i1(i + 1, 2) = SQR(Ion_X[satno].time[i]); B_i1(i + 1, 3) = SQR(Ion_X[satno].time[i]) * Ion_X[satno].time[i];

		L_i1(i + 1, 1) = Ion_X[satno].gf_value[i];
	}



	LeastSquare(B_i, L_i, P_i, V_i, nn, x_i);
	LeastSquare(B_i1, L_i1, P_i1, V_i1, nn1, x_i1);
	//B_i.Output_to_File("E:\\B_i", 1);
	//L_i.Output_to_File("E:\\L_i", 1);
	//x_i.Output_to_File("E:\\x_i", 1);
	Ion_X[satno].a = x_i1(1, 1);
	Ion_X[satno].b = x_i1(2, 1);
	Ion_X[satno].c = x_i1(3, 1);

	//Ion_X[satno].a1 = x_i1(1, 1);
	//Ion_X[satno].b1 = x_i1(2, 1);
	//Ion_X[satno].c1 = x_i1(3, 1);
	double p_30 = 0;
	double p_60 = 0;
	double psd = 0;
	double p_30_mean = 0;
	double p_60_mean = 0;
	vector<double>dion_pre;
	vector<double>ion11;
	vector<double>ion22;

#if (0)
	if (satno == numb)  //输出PDS参数
	{
		FILE* fp = NULL;
		fp = fopen("E:\\psd.out", "a+");
		fprintf(fp, "%15.9f", psd1);
		fprintf(fp, "%15.9f \n", psd);
		fclose(fp);
	}

	//if (satno == numb)
	//{
	//	FILE* fp = NULL;
	//	fp = fopen("E:\\gfion.out", "a+");

	//	for (int i = 0; i < Ion_X[satno].gf_value.size(); i++)
	//	{
	//		fprintf(fp, "%15.9f", Ion_X[satno].gf_value[i]);

	//		fprintf(fp, "%15.9f\n", dion_pre[i]);
	//	}

	//	//fprintf(fp, "%15.9f\n", Ion_X[satno].gf_value[Ion_X[satno].gf_value.size()-1]);
	//	//fprintf(fp, "%15.9f \n", x_i1(1, 1) + x_i1(2, 1) * flag_epoch_count + x_i1(3, 1) * SQR(flag_epoch_count));
	//	fclose(fp);
	//}
	//if (satno == numb)
	//{
	//	FILE* fp = NULL;
	//	fp = fopen("E:\\apsd.out", "a+");
	///*	fprintf(fp, "%15.9f", apsd30);
	//	fprintf(fp, "%15.9f", apsd60);*/
	//	//fprintf(fp, "%15.9f", apsd30_);
	//	//fprintf(fp, "%15.9f \n", apsd60_);
	//	fclose(fp);
	//}

	//if (satno == numb)
	//{
	//	FILE* fp = NULL;
	//	fp = fopen("E:\\ele.out", "a+");
	//	for (int i = 0; i < Ion_X[satno].ele.size(); i++)
	//	{
	//		fprintf(fp, "%15.9f\n", Ion_X[satno].ele[i]*R2D);
	//	}
	//	fclose(fp);
	//}

	//if (satno == numb)
	//{
	//	FILE* fp = NULL;
	//	FILE* fp1 = NULL;
	//	fp = fopen("E:\\EDI.out", "a+");
	//	//fp1 = fopen("E:\\dgfion1.out", "a+");
	//	double d = Ion_X[satno].dgf_value[Ion_X[satno].dgf_value.size() - 1];
	//	
	//	double c = edi_5_[edi_5_.size()-1];
	//	for (int i = 0; i < Ion_X[satno].dgf_value.size(); i++)
	//	{
	//		/*if (i<interval)
	//		{
	//			continue;
	//		}*/
	//		fprintf(fp, "%15.9f\n", Ion_X[satno].dgf_value[i]);
	//	}
	//	for (int i = 0; i < edi_5.size(); i++)
	//	{
	//		fprintf(fp1, "%15.9f\n", edi_5[i]);
	//	}
	//	//fprintf(fp, "%15.9f", d);
	//	//fprintf(fp, "%15.9f \n", c); 
	//	fclose(fp);
	//	//fclose(fp1);
	//	//d = c;
	//}
#endif
	sum_30 = 0;
	s_30 = 0, s_30_ = 0;
	s_60 = 0, s_60_ = 0;  k = 0, k1 = 0;
	sum = 0;
	edi_5.clear();
	edi_5_.clear();

	//}
}

int LeastSquarePPP::detslp_GF13(ppp_option_t& popt, PPPObsData& pppdata, PreciseData& predata, GenData& gendata)
{
	int satN = pppdata.m_SatCount; int g_or_c = 2, orbn = 0;
	double SatPCV[NFREQ] = { 0.0 }, RecPCV[NFREQ] = { 0.0 }; double satpos[3] = { 0.0 }; double PhaseWindup = 0.0;
	double rec_ant_delta[3] = { Antdiff[1], Antdiff[0], Antdiff[2] };
	double thres = 0.0;
	for (int i = 0; i < satN; i++)
	{
		double elev = m_SatInfo.m_Ele[i]; unsigned int satno = pppdata.m_Prn[i];
		if (GetSystem(satno) == 'G' || GetSystem(satno) == 'E') { g_or_c = 2; orbn = 0; }
		if (GetSystem(satno) == 'C') { g_or_c = 3; orbn = 0; }
		if (GetSystem(satno) == 'R') { orbn = Get_Nav_Orbitn(satno, predata.m_NavData); g_or_c = 2; }
		double wavelength1 = Get_WaveLength(satno, 1, orbn);  //测试 
		double wavelength2 = Get_WaveLength(satno, g_or_c, orbn);
		double wavelength3 = Get_WaveLength(satno, 5, orbn);
		double azel[2] = { m_SatInfo.m_Azi[i], m_SatInfo.m_Ele[i] };
		satpos[0] = m_SatInfo.m_Xs[3 * i + 0];
		satpos[1] = m_SatInfo.m_Xs[3 * i + 1];
		satpos[2] = m_SatInfo.m_Xs[3 * i + 2];
		SatAntPCV(GetSystem(satno), satno, satpos, m_Sta_tmp, gendata.m_pcvs, SatPCV);
		RecAntModel(GetSystem(satno), &gendata.m_pcvr, rec_ant_delta, azel, 1, RecPCV);
		WindupCorr(popt.SunPos, satpos, m_Sta_tmp, PhaseWindup);

		double P1 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr1[i]);
		double P2 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr2[i]);
		double P3 = PseudoRange_PCV_Correct(RecPCV[0], SatPCV[0], pppdata.m_Pr3[i]);
		double L1 = CarrierPhase_PCV_Correct(wavelength1, RecPCV[0], SatPCV[0], pppdata.m_L1[i]);
		double L2 = CarrierPhase_PCV_Correct(wavelength2, RecPCV[0], SatPCV[0], pppdata.m_L2[i]);
		double L3 = CarrierPhase_PCV_Correct(wavelength3, RecPCV[0], SatPCV[0], pppdata.m_L3[i]);

		if (P3 != 0 && L3 != 0) {
			L1 -= wavelength1 * PhaseWindup;
			L2 -= wavelength2 * PhaseWindup;
			L3 -= wavelength3 * PhaseWindup;
			gf_cur13[satno] = wavelength1 * L1 - wavelength3 * L3;
		}
		else
		{
			continue;
		}
		if (flag_epoch_count == 0)continue;

		if (P3 != 0 && L3 != 0)
		{
			if (gf_cur13[satno] == 0 || gf_pre13[satno] == 0)
			{
				continue;
			}

		}
		if (elev < popt.ele_mask * R2D) elev = popt.ele_mask * R2D;
		if (elev >= 15.0)thres = THresGF;
		else thres = -THresGF / 15.0 * elev + 2 * THresGF;

		if (fabs(gf_cur13[satno] - gf_pre13[satno]) > min(thres, 1.5))
		{

			FILE* fp = NULL;
			fp = fopen("E:\\gf13.out", "a+");


			//fprintf(fp, "%d%d \n", flag_epoch_count, satno);
			fprintf(fp, "%d \n", satno);
			fclose(fp);
			
			_Cycle_Slip.push_back(satno);
		}


	}
	if (_Cycle_Slip.size() > 1)//删除相同的元素
	{
		sort(_Cycle_Slip.begin(), _Cycle_Slip.end());
		_Cycle_Slip.erase(unique(_Cycle_Slip.begin(), _Cycle_Slip.end()), _Cycle_Slip.end());
	}

	return 0;
}

void deletecurrentobsdata(vector<zhydata>& data, int index)
{
	if (data[1].m_SatCount == 0)return;

	data[1].m_L1[index] = data[1].m_L1[data[1].m_SatCount - 1];
	data[1].m_L1[data[1].m_SatCount - 1] = 0;
	data[1].m_P1[index] = data[1].m_P1[data[1].m_SatCount - 1];
	data[1].m_P1[data[1].m_SatCount - 1] = 0;
	data[1].m_Prn[index] = data[1].m_Prn[data[1].m_SatCount - 1];
	data[1].m_Prn[data[1].m_SatCount - 1] = 0;
	--data[1].m_SatCount;
}
int LeastSquarePPP::LSQ_qc(MatrixT& b, MatrixT& l, MatrixT& p, MatrixT& v, MatrixT& nn, MatrixT& x, int ind, unsigned int prn)
{
	LeastSquare(b, l, p, v, nn, x);
	int n = b.GetCol(); int m = b.GetRow();
	MatrixT vtpv = v.Trans() * p * v;
	double val = 0; double thres = 0; int info = 0;
	if (ind == 0)
	{
		val = vtpv.GetValue(1, 1) / (m - n);
#if (1)
		thres = chisqr_arr[m - n - 1] / (m - n);
#else
		thres = 3;
#endif
	}
	else
	{
		val = vtpv.GetValue(1, 1) / double(m);
		thres = 35;
	}

	if (val > thres)
	{
		info = 1;
		FILE* fp = NULL;
		fp = fopen("E:\\vtpv.out", "a+");


		fprintf(fp, "%d , %d ,%15.9f \n", flag_epoch_count, prn, val);
		fclose(fp);
		/*b.Output_to_File("E:\\b.dat",1);
		l.Output_to_File("E:\\l.dat", 1);
		p.Output_to_File("E:\\p.dat", 1);
		v.Output_to_File("E:\\v.dat", 1);
		x.Output_to_File("E:\\x.dat", 1);
		nn.Output_to_File("E:\\nn.dat", 1);*/
	}



	return info;
}
void LeastSquarePPP::QC_Scdia(vector<double> L_P, vector<double> L_C, vector<double> P_P, vector<double> P_C, vector<int> Obstype, unsigned int prn, PreciseData& predata, zhydata data)
{
	int nf = L_P.size();
	int m = 2 * nf + 1;
	//int n = 2 + nf; 
	int n = 2;
	MatrixT Z, R, H, vp, nn, xp;
	Z.Resize(m, 1);
	R.Resize(m, m);
	H.Resize(m, n); int i = 0;
	for (int f = 0; f < nf; f++)
	{
		Z(i + 1, 1) = L_P[f] - L_C[f];  R(i + 1, i + 1) = 1 / (2 * SQR(0.003));
		H(i + 1, 1) = 1;
		H(i + 1, 2) = 0 - get_ion_xishu(Obstype[f], prn, predata.m_NavData);
		//H(i + 1, 3 + f) = get_wth(prn,Obstype[f],  predata.m_NavData);
		i++;
		Z(i + 1, 1) = P_P[f] - P_C[f];  R(i + 1, i + 1) = 1 / (2 * SQR(0.3));

		H(i + 1, 1) = 1;
		H(i + 1, 2) = get_ion_xishu(Obstype[f], prn, predata.m_NavData);
		i++;

	}
	{
		Z(i + 1, 1) = gf_ion[prn];	R(i + 1, i + 1) = 1 / SQR(gf_sigemaIon[prn]);	H(i + 1, 2) = 1;
		i++;
	}

	//H.Output_to_File("E:\\H.out", 1);
	//Z.Output_to_File("E:\\Z.out", 1);
	//R.Output_to_File("E:\\R.out", 1);

	int ind = LSQ_qc(H, Z, R, vp, nn, xp, 1, prn);

	//vp.Output_to_File("E:\\vp.out", 1);
	//nn.Output_to_File("E:\\nn.out", 1);
	//xp.Output_to_File("E:\\xp.out", 1);

	if (ind == 0)
	{
		return;
	}
#if (0)
	n = 2 + nf;
	H.Resize(m, n);

	for (int f = 0; f < nf; f++)
	{
		Z(i + 1, 1) = L_P[f] - L_C[f];  R(i + 1, i + 1) = 1 / (2 * SQR(0.003));
		H(i + 1, 1) = 1;
		H(i + 1, 2) = 0 - get_ion_xishu(Obstype[f], prn, predata.m_NavData);
		H(i + 1, 3 + f) = get_wth(prn, Obstype[f], predata.m_NavData);
		i++;
		Z(i + 1, 1) = P_P[f] - P_C[f];  R(i + 1, i + 1) = 1 / (2 * SQR(0.3));

		H(i + 1, 1) = 1;
		H(i + 1, 2) = get_ion_xishu(Obstype[f], prn, predata.m_NavData);
		i++;

	}
	{
		Z(i + 1, 1) = gf_ion[prn];	R(i + 1, i + 1) = 1 / SQR(gf_sigemaIon[prn]);	H(i + 1, 2) = 1;
		i++;
	}
	int ind1 = LSQ_qc(H, Z, R, vp, nn, xp, 0, prn);

	MatrixT Nx, Np, F;
	Nx.Resize(nf, 1);
	Np.Resize(nf, nf);
	F.Resize(nf, 2);
	for (int k = 0; k < nf; k++)
	{
		for (int j = 0; j < nf; j++)
		{
			Nx(k + 1, 1) = xp(k + 3, 1);
			Np(k + 1, j + 1) = nn(k + 3, j + 3);
		}
	}

	double s[2]; bool pass;
	//Nx.Output_to_File("E:\\Nx.out",1);
	//Np.Output_to_File("E:\\Np.out", 1);
	lambda(nf, 2, Nx.m_mat.data(), Np.m_mat.data(), F.m_mat.data(), s, 0.001, pass);

	FILE* fp = NULL;
	fp = fopen("E:\\F.out", "a+");
	fprintf(fp, "%d, \n", flag_epoch_count);
	fclose(fp);
	F.Output_to_File("E:\\F.out", 1);
	int index_ = 0; int kk = 0;
	if (pass) {
		for (int i = 0; i < data.m_SatCount; i++)
		{
			if (data.m_Prn[i] == prn)
			{
				for (int k = 0; k < i; k++)
				{
					index_ += num_zhy_pre[k];
					kk++;
				}
			}
		}
		vector<double>a;
		for (int i = 0; i < nf; i++)
		{
			a.push_back(F.m_mat.data()[i]);
			X_Trans(_ParaN_const_pre + data.m_SatCount + 1 + i, 1) += F.m_mat.data()[i];
		}

	}
	else
#endif
	{
		
		_Cycle_Slip.push_back(prn);
	}



	return;
}