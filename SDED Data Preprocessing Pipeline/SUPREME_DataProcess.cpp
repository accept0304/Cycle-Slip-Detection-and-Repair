#include "SUPREME_DataProcess.h"
#include "SUPREME_Anntenna.h"
#include "SUPREME_CommonFunction.h"
#include "SUPREME_OrbClk.h"
#include "SUPREME_QC.h"
#include "SUPREME_Troposphere.h"
#include "iostream"
#include <string>

static double Carrier1ObsOffset[MAX_GNSS_PRN_NUM] = { 0 };
static double Carrier2ObsOffset[MAX_GNSS_PRN_NUM] = { 0 };

//bool cmp(string a, string b); bool equal(string a, string b);
void SetCodePri(int sys, int freq, const char *pri)
{
	if (freq <= 0 || MAXFREQ<freq) return;
	if (sys&PRO_SYS_GPS) strcpy(CODE_PRIS[0][freq - 1], pri);
	if (sys&PRO_SYS_GLO) strcpy(CODE_PRIS[1][freq - 1], pri);
	if (sys&PRO_SYS_GAL) strcpy(CODE_PRIS[2][freq - 1], pri);
	if (sys&PRO_SYS_QZS) strcpy(CODE_PRIS[3][freq - 1], pri);
	if (sys&PRO_SYS_SBS) strcpy(CODE_PRIS[4][freq - 1], pri);
	if (sys&PRO_SYS_BDS) strcpy(CODE_PRIS[5][freq - 1], pri);
	if (sys&PRO_SYS_IRN) strcpy(CODE_PRIS[6][freq - 1], pri);
}

SppOutputPara::SppOutputPara()
{
	m_Size = 0;
	Spos = NULL; Sclk = NULL;
	Trop = NULL; Elev = NULL; Azim = NULL; Resp = NULL;
}

SppOutputPara::SppOutputPara(unsigned int satn)
{
	m_Size = satn;

	Rclk_spp = 0.0; RclkCov_spp = 0.0;
	StaPos_spp[0] = 0.0; StaPos_spp[1] = 0.0; StaPos_spp[2] = 0.0;
	StaPosCov_spp[0] = 0.0; StaPosCov_spp[1] = 0.0; StaPosCov_spp[2] = 0.0;

	Spos = (double*)(malloc(sizeof(double)*satn * 3));    // satellite position
	Sclk = (double*)(malloc(sizeof(double)*satn));        // satellite clock offset
	Trop = (double*)(malloc(sizeof(double)*satn));        // troposphere delay
	Elev = (double*)(malloc(sizeof(double)*satn));        // satellite elevation
	Azim = (double*)(malloc(sizeof(double)*satn));        // satellite azimuth
	Resp = (double*)(malloc(sizeof(double)*satn));        // residuals

	memset((double*)Spos, 0.0, sizeof(double)*satn * 3);
	memset((double*)Sclk, 0.0, sizeof(double)*satn);
	memset((double*)Trop, 0.0, sizeof(double)*satn);
	memset((double*)Elev, 0.0, sizeof(double)*satn);
	memset((double*)Azim, 0.0, sizeof(double)*satn);
	memset((double*)Resp, 0.0, sizeof(double)*satn);
}

SppOutputPara::~SppOutputPara()
{
	m_Size = 0;
	if (Spos) free(Spos);
	if (Sclk) free(Sclk);
	if (Trop) free(Trop);
	if (Elev) free(Elev);
	if (Azim) free(Azim);
	if (Resp) free(Resp);
}

void SppOutputPara::DeleSat(unsigned int i)
{
	Spos[3 * i + 0] = Spos[3 * (m_Size - 1) + 0];
	Spos[3 * i + 1] = Spos[3 * (m_Size - 1) + 1];
	Spos[3 * i + 2] = Spos[3 * (m_Size - 1) + 2];
	Sclk[i] = Sclk[m_Size - 1];
	Trop[i] = Trop[m_Size - 1];
	Elev[i] = Elev[m_Size - 1];
	Azim[i] = Azim[m_Size - 1];
	Resp[i] = Resp[m_Size - 1];
	--m_Size;
}

PPPObsData::PPPObsData(unsigned int satn)
{
	Initialize(satn);
}

PPPObsData::~PPPObsData()
{
	Clear();
}

int PPPObsData::Initialize(unsigned int satn)
{
	m_Valid_flag = (unsigned int*)(malloc(sizeof(unsigned int*)*satn));
	m_Prn = (unsigned int*)(malloc(sizeof(unsigned int)*satn));
	m_Pr1 = (double*)(malloc(sizeof(double)*satn));
	m_Pr2 = (double*)(malloc(sizeof(double)*satn));
	m_Pr3 = (double*)(malloc(sizeof(double) * satn));
	m_Pr4 = (double*)(malloc(sizeof(double) * satn));
	m_Pr5 = (double*)(malloc(sizeof(double) * satn));
	m_L1 = (double*)(malloc(sizeof(double) * satn));
	m_L2 = (double*)(malloc(sizeof(double) * satn));
	m_L3 = (double*)(malloc(sizeof(double) * satn));
	m_L4 = (double*)(malloc(sizeof(double) * satn));
	m_L5 = (double*)(malloc(sizeof(double) * satn));
	m_D1 = (double*)(malloc(sizeof(double) * satn));
	_num = (int*)(malloc(sizeof(double) * satn));
	_num_accu = (int*)(malloc(sizeof(double) * satn));
	
	memset((unsigned int*)m_Valid_flag, 0, sizeof(unsigned int)*satn);
	memset((unsigned int*)m_Prn, 0, sizeof(unsigned int)*satn);
	memset((double*)m_Pr1, 0.0, sizeof(double)*satn);
	memset((double*)m_Pr2, 0.0, sizeof(double)*satn);
	memset((double*)m_Pr3, 0.0, sizeof(double) * satn);
	memset((double*)m_Pr4, 0.0, sizeof(double) * satn);
	memset((double*)m_Pr5, 0.0, sizeof(double) * satn);
	memset((double*)m_L1, 0.0, sizeof(double)*satn);
	memset((double*)m_L2, 0.0, sizeof(double)*satn);
	memset((double*)m_L3, 0.0, sizeof(double) * satn);
	memset((double*)m_L4, 0.0, sizeof(double) * satn);
	memset((double*)m_L5, 0.0, sizeof(double) * satn);
	memset((double*)m_D1, 0.0, sizeof(double) * satn);
	memset((int*)_num, 0.0, sizeof(int) * satn);
	memset((int*)_num_accu, 0.0, sizeof(int) * satn);
	
	cur_num = 0;
	return 1;
}

int PPPObsData::Clear()
{
	if (m_Valid_flag) free(m_Valid_flag);
	if (m_Prn) free(m_Prn);
	if (m_Pr1) free(m_Pr1);	if (m_Pr2) free(m_Pr2);	if (m_Pr3) free(m_Pr3);
	if (m_Pr4) free(m_Pr4); if (m_Pr5) free(m_Pr5);
	if (m_L1) free(m_L1);	if (m_L2) free(m_L2); 
	if (m_L3) free(m_L3); if (m_L4) free(m_L4); if (m_L5) free(m_L5);
	if (m_D1) free(m_D1); if (_num){ free(_num); }if (_num_accu){ free(_num_accu); }
	
	
	return 1;
}

int PPPObsData::P_detec(GenData* gendata, ppp_option_t& popt, unsigned int satno)
{
	double p1p2, p1c1, p2c2;
	double p1_ = 0, p2_ = 0, p3_ = 0, p4_ = 0, p5_ = 0; string type; string nul = "  ";
	DCB_Item* PtrDcb = gendata->m_Dcb.Find_mgex_DCB_Data(satno); 
	if (PtrDcb == NULL)return  1;
	int num = epoch_otype.size();
	if (GetSystem(satno) == 'G') 
	{
		if (epoch_otype[0] != "C1P") 
		{
			type = "C1C" + nul + epoch_otype[0];
			double value = find_Dcb_value(type, PtrDcb);
			p1_ = p1 + value;
			p1c1 = Get_DCB_P1C1(satno);
			p1_ += p1c1 * DCB_Parameter;  
			p1 = p1_;
		}
		else { p1_ = p1; }
		if (epoch_otype[1] != "C2P")
		{
			double value = 0, value1 = 0, value2 = 0;
			if (epoch_otype[1] != "C2W")
			{
				type =   "C2W"+nul +epoch_otype[1];
				value = find_Dcb_value(type, PtrDcb);
				p2_ = p2 + value;
			}
			type = "C2C" + nul + epoch_otype[1];
			value1 = find_Dcb_value(type, PtrDcb);
			p2_ = p2 + value+ value1;
			p2c2 = Get_DCB_P2C2(satno);
			p2_ += p2c2 * DCB_Parameter;   
			p2 = p2_;
			p1p2 = Get_DCB_P1P2(satno);
			p2_ += p1p2 * DCB_Parameter;   
		}
		else { p2_ = p2; }
		if(num==3)
		{
			type= "C1C" + nul + epoch_otype[2];
			double value = find_Dcb_value(type, PtrDcb);  
			p3_ = p3 + value;   
			p1c1 = Get_DCB_P1C1(satno);
			p3_ += p1c1 * DCB_Parameter;      
			p1p2 = 0 - Get_DCB_P1P2(satno);
			p3 = p3 + (value + (G_B_12 / G_B_13) * p1p2);
#if 1																		
			double ifgc = gendata->IFCL_bias.get_g_ifcb(m_obst, satno);
			
			
			l3 = l3 + ifgc * ((0 - 1) / G_B_13);
#endif
		}
	}
	else if ((GetSystem(satno) == 'C')) 
	{
		if (epoch_otype[0] == "C2I")
		{
			type = epoch_otype[0] + nul + "C6I";
			double value = find_Dcb_value(type, PtrDcb);
			p1_ = p1 + value;  //2->6
		}
		p2_ = p2;
		for (int j = 0; j < num; j++)
		{
			if (strstr(epoch_otype[j].c_str(), "7")) 
			{
				type = "C2I" + nul + epoch_otype[j];
				double value1 = find_Dcb_value(type, PtrDcb);
				p3_ = p3 + value1;	//7->2
				type = "C2I  C6I";
				double value= find_Dcb_value(type, PtrDcb);
				p3_ -= value;		//2->6
				p3 = p3 + (value1 + (C_B_12 / C_B_13) * (0-value));
			}
			if (strstr(epoch_otype[j].c_str(), "1")) 
			{
				type = epoch_otype[j] + nul + "C6I";
				double value = find_Dcb_value(type, PtrDcb);
				type = "C2I  C6I";
				double value1 = find_Dcb_value(type, PtrDcb);
				p4_ = p4 - value;    //1->6
				p4 = p4 + (value1 - value + (C_B_12 / C_B_14) * (0-value1));
			}
			if (strstr(epoch_otype[j].c_str(), "5"))
			{
				double value = 0, value1 = 0, value2 = 0;
				if (strstr(epoch_otype[j].c_str(), "X")) 
				{
					type = "C1X" + nul + epoch_otype[j];
					value=find_Dcb_value(type, PtrDcb);
					type = "C1X  C6I";
					value1 = find_Dcb_value(type, PtrDcb);
					value -= value1;
				}
				else if (strstr(epoch_otype[j].c_str(), "P"))
				{
					type = "C1P" + nul + epoch_otype[j];
					value = find_Dcb_value(type, PtrDcb);
					type = "C1P  C6I";
					value1 = find_Dcb_value(type, PtrDcb);
					value -= value1;
				}
				else if (strstr(epoch_otype[j].c_str(), "D"))
				{
					type = "C1D" + nul + epoch_otype[j];
					value = find_Dcb_value(type, PtrDcb);
					type = "C1D  C6I";
					value1 = find_Dcb_value(type, PtrDcb);
					value -= value1;
				}
				p5_ = p5 + value;    //1->6
				type= "C2I  C6I";
				value2 = find_Dcb_value(type, PtrDcb);
				value += value2;
				p5 = p5 + (value + (C_B_12 / C_B_15) * (0 - value2));
			}
		}
	}
#if (1)
	else if ((GetSystem(satno) == 'E')) 
	{
		p1_ = p1;
		for (int j = 0; j < num; j++)
		{
			if (strstr(epoch_otype[j].c_str(), "5")) 
			{
				type = "C1X" + nul + epoch_otype[j];
				double value = find_Dcb_value(type, PtrDcb);
				p2_ =p2+ value;
			}
			if (strstr(epoch_otype[j].c_str(), "7"))
			{
				type = "C1X" + nul + epoch_otype[j];
				double value = find_Dcb_value(type, PtrDcb);
				p3_ = p3 + value;
				type = "C1X  C5X";
				double value1 = find_Dcb_value(type, PtrDcb);
				p3 = p3 + (value + (E_B_12 / E_B_13) * (0 - value1));
			}
			if (strstr(epoch_otype[j].c_str(), "8"))
			{
				type = "C1X" + nul + epoch_otype[j];
				double value = find_Dcb_value(type, PtrDcb);
				p4_ = p4 + value;
				type = "C1X  C5X";
				double value1 = find_Dcb_value(type, PtrDcb);
				p4 = p4 + (value + (E_B_12 / E_B_14) * (0 - value1));
			}
			/*if (strstr(epoch_otype[j].c_str(), "6"))    //DCB 文件中没有  可以考虑每个卫星单独加参数
			{
				type = "C1X" + nul + epoch_otype[j];
				double value = find_Dcb_value(type, PtrDcb);
				p5_ = p5 + value;
				type = "C1X  C5X";
				double value1 = find_Dcb_value(type, PtrDcb);
				p5 = p5 + (value + (E_B_12 / E_B_15) * (0 - value1));
			}*/
		}
		
	}
#endif
	else if ((GetSystem(satno) == 'R')) 
	{
		p1_ = p1; p2_ = p2; p3_ = p3; p4_ = p4; p5_ = p5;
	}
	
	int a = 0;
	double key1 = 0, key2 = 0, key3 = 0, key4 = 0;
	key1 = fabs(p1_ - p2_);
	if (p3_)key2 = fabs(p1_ - p3_); if (p4_)key3 = fabs(p1_ - p4_); if (p5_)key4 = fabs(p1_ - p5_);
	if (key1) 
	{
		key1 < 30 ? a = 1 : a = 0;
	}
	if (key2)
	{
		key2 < 30 ? a = 1 : a = 0;
	}
	if (key3)
	{
		key3 < 30 ? a = 1 : a = 0;
	}
	if (key4)
	{
		key4 < 30 ? a = 1 : a = 0;
	}
	if (a) return 1;
	else return 0;
}

/// Check observation data
int PPPObsData::CheckObsData(ppp_option_t& popt, GenData& gendata, unsigned int satno, gnsstime& t, PreciseData *predata, ClockData *clkdata)
{
	char Sysflag = GetSystem(satno);

	// Check DCB data ----------------------------------------
	if ((Sysflag == 'G'))//|| (Sysflag == 'R')
	{
		DCB_Item *PtrDcb = gendata.m_Dcb.Find_DCB_Data(satno);
		if (PtrDcb == NULL)	return 0;
	}

	if (Sysflag == 'R')
	{
		if (Get_Nav_Orbitn(satno, predata->m_NavData) == -1)
		{
			//Log
			return 0;
		}
	}

	if (Sysflag == 'C')//北斗dcb改正；
	{
		DCB_Item* PtrDcb = gendata.m_Dcb.Find_mgex_DCB_Data(satno);
		if (PtrDcb == NULL)	return 0;//
	}
	if (Sysflag == 'E')
	{
		//
	}
	// Check DCB data end ------------------------------------

	// Check Ephemeris and Clock data ----------------------------------------------------------
	if (popt.process_time == MODE_POST && popt.process_mode==MODE_PPP)
	{
		static double clk_value[1000] = { 0 };

		if (popt.use_orbclk)
		{
			// Check precise clock --------------------------------------------------------------
			if (clkdata == NULL) return 0;
			Clk_Data_Satellite *cp = clkdata->Search_clock_Prn(satno);
			if (cp == NULL)	return 0;
			else
			{
				double clk_tmp = 0;
				if ((int)(clk_tmp = Clock_Interpolate_w(satno, t, clkdata, 10) + 0.5) == 999999)
				{
					//Log
					return 0;
				}
				else
				{
					if (clk_value[satno] == 0)
						clk_value[satno] = clk_tmp;
					else
					{
						if (fabs(clk_value[satno] * LIGHTSPEED - clk_tmp*LIGHTSPEED) > 100)
						{
							clk_value[satno] = clk_tmp;
							//Log
							return 0;
						}
					}
				}
			}
			// -------------------------------------------------------------------------------------

			// Check precise ephemeris -------------------------------------------------------------
			if (predata == NULL) return 0;
			PreDataSat *sp = predata->Precise_Search_Prn(satno);
			if (sp == NULL)
			{
				//Log
				return 0;
			}
			else
			{
				double satpos[3] = { 0.0 };
				if (predata->Precise_Interp_Pos_w(satno, t, satpos, 10) == 0)
				{
					//Log
					return 0;
				}
			}
			// -------------------------------------------------------------------------------------
		}
	}
	// Check Ephemeris and Clock data ----------------------------------------------------------

	return 1;
}

/* Get Observation Value ----------------------------------------------------
* flag<0:get p1(c1),p2,L1,L2
* flag=1:get p1(c1),p2
* flag=2:get L1,L2
* --------------------------------------------------------------------------- */
int PPPObsData::GetBlockData(ppp_option_t& popt, GenData *gendata, ObsEpochData& obsdata,
	PreciseData *predata, ClockData *clkdata, int flag)
{
	char sys = 0;
	unsigned int i = 0, k = 0, SatCount = 0;
	unsigned int FreqN = popt.Freqn;
	int mathmode = popt.math_model;
	SatCount = obsdata.m_SatCount;
	/*vector<int>ifb_num1;*/
	for (i = 0; i < SatCount; ++i)
	{
		unsigned int gf_flag = 0, mw_flag = 0;
		int sysflag = 0;

		switch (sys = GetSystem(obsdata.EpochData[i].m_Prn))
		{
			case 'G':
			{
				if (popt.IS_ISB == 1) { sysflag = 0; break; }
				else { sysflag = popt.System & PRO_SYS_GPS ? 0 : 1; break; }
			}
			case 'R':
			{
				if (popt.IS_ISB == 1) { sysflag = 0; break; }
				else { sysflag = popt.System & PRO_SYS_GLO ? 0 : 1; break; }
			}
			case 'C':
			{
				if (popt.IS_ISB == 1) { sysflag = 0; break; }
				else { sysflag = popt.System & PRO_SYS_BDS ? 0 : 1; break; }
			}
			case 'E':
			{
				if (popt.IS_ISB == 1) { sysflag = 0; break; }
				else { sysflag = popt.System & PRO_SYS_GAL ? 0 : 1; break; }
			}
			default:sysflag = 1; break;
		}

		if (sysflag) continue;
		if (CheckObsData(popt, *gendata, obsdata.EpochData[i].m_Prn, obsdata.m_EpochTime, predata, clkdata) == 0)//自己加上北斗DCB； 
		{
			//Log
			continue;
		}

		switch (sys)
		{
			case 'G':
			{
				int freq = 2;
				if (FreqN >= 3) { freq = 3; }
						
						if (GetGPSData(gendata->m_Dcb, obsdata, obsdata.EpochData[i], freq, flag,mathmode) == 0)
							continue;
						break;
			}
			case 'R':
			{
						if (GetGLOData(gendata->m_Dcb, obsdata, obsdata.EpochData[i], FreqN, flag,mathmode) == 0)
							continue;
						break;
			}
			case 'E':
			{
						if (GetGALData(gendata->m_Dcb, obsdata, obsdata.EpochData[i], FreqN, flag, mathmode) == 0)
							continue;
						break;
			}
			case 'C':
			{
						if (obsdata.EpochData[i].m_Prn <= 316)continue;
						if (GetBDSData(gendata->m_Dcb, obsdata, obsdata.EpochData[i], FreqN, flag, obsdata.EpochData[i].m_Prn, mathmode) == 0)
							continue;
						break;
			}	
		default: continue;
		}
#if 0
		if ((mathmode == 0 || mathmode == 4))
		{
			if (!(P_detec(gendata, popt, obsdata.EpochData[i].m_Prn))) 
			{ 
				
				Prn_Otype.pop_back();
				Prn_Otype_int.pop_back();
				Prn_Otype_channel.pop_back();
				epoch_otype.pop_back();
				continue;
			}
		}
#endif
		
		m_Prn[k] = obsdata.EpochData[i].m_Prn;
		m_Pr1[k] = p1; m_Pr2[k] = p2; m_Pr3[k] = p3; m_Pr4[k] = p4; m_Pr5[k] = p5;
		m_L1[k]  = l1; m_L2[k]  = l2; m_L3[k]  = l3; m_L4[k]  = l4; m_L5[k]  = l5;          
		m_D1[k] = d1; _num[k] = obstype_num;
		 cur_num += obstype_num; popt.CUR_mum = cur_num;
		_num_accu[k] = cur_num - obstype_num;
		// Carrier Phase Offset ---------------------------------------------------------
		unsigned int satno = m_Prn[k];
		if ((Carrier1ObsOffset[satno] == 0) && (l1 != 0 || l2 != 0))
		{			
			int orbn = 0;  // Get orbn for GLONASS satellites
			if (sys == 'R')	{
				orbn = Get_Nav_Orbitn(satno, predata->m_NavData);
				if (orbn == -1) continue;
			}
		}
		
		m_Valid_flag[k] = 1;
		k++;
	}

	
	m_SatCount = k;

	return k;
}

/// Get gps data

int PPPObsData::GetGPSData(CodeDCB &dcbdata, ObsEpochData& obsdata, SatelliteData &epochdata, unsigned int freqn, int flag,int mode)
{
	fl = 0, f2 = 0, f3 = 0, f4 = 0, f5 = 0, pflag = 0, pflag2 = 0, pflag3 = 0, pflag4 = 0, pflag5 = 0;
	c1 = 0, p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0;

	d1 = d2 = d3 = 0; obstype_num = 1; vector<string> obstype; vector<int> obstypeint; vector<string>channel; 
	obstype.push_back(to_string(epochdata.m_Prn));
	if (flag < 2)
	{
		/// Rinex 2.0 version
		if (obsdata.m_ObsType->m_Version < 3)
		{
			p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, P1);
			fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, P1);

			if ((int)(p1 + 0.5) == 0)
			{
				c1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1);
				fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1);
				if ((int)(c1 + 0.5) == 0) return 0;
				else
				{
					double dcb = dcbdata.Get_DCB_P1C1(epochdata.m_Prn);
					p1 = c1 + dcb*LIGHTSPEED*1.0e-9;
				}
			}

			if (fl == 1) return 0;

			/// pseudorange quality check
			c1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1);
			if (c1 != 0 && p1 != 0)
			{
				if (PseudoRange_C1P1_Detect(c1, p1)) return 0;
			}
		}
		else /// Rinex 3.0 version
		{
			 
			p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1P);
			p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C2W);
			
			fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1P);
			f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C2W);
			
			d1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, D1P);
			
			if ((int)(p1 + 0.5) == 0)
			{
				p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1Y);
				fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1Y);
				
			}
			if ((int)(p1 + 0.5) == 0)
			{
				p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1W);
				fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1W);	
			}
			if ((int)(p1 + 0.5) == 0)
			{
				p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1X);
				fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1X);
			}
			if ((int)(p1 + 0.5) == 0)
			{
				c1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1C);
				fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1C);
				if ((int)(c1 + 0.5) == 0) return 0;
				else
				{
					double dcb = dcbdata.Get_DCB_P1C1(epochdata.m_Prn);
					p1 = c1 + dcb*LIGHTSPEED*1.0e-9;
				}
			}
			if (fl == 1) return 0;

			/// pseudorange quality check
			c1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1C);
			if (c1 != 0 && p1 != 0)
			{
				if (PseudoRange_C1P1_Detect(c1, p1)) return 0;
			}
			if ((int)(p2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C2P);
				f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C2P);
			}
			if ((int)(p2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C2X);
				f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C2X);
			}
			if ((int)(p2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C2L);
				f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C2L);
			}
			if ((int)(p2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C2W);
				f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C2W);
			}
			if ((int)(p2 + 0.5) == 0)
			{
				c1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C2C);
				f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C2C);
				if ((int)(c1 + 0.5) == 0); //return 0;
				else
				{
					double dcb = dcbdata.Get_DCB_P2C2(epochdata.m_Prn);
					p2 = c1 + dcb * LIGHTSPEED * 1.0e-9;
				}
			}
			if (((int)(p2 + 0.5) == 0 && (mode == 0)) || (f2 == 1))return 0;//||mode == 4
			if ((int)(p3 + 0.5) == 0 && mode == 4 ) {
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C5Q);
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C5Q);
			}//&&freqn>=3
			if ((int)(p3 + 0.5) == 0 && mode == 4 )
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C5X);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C5X);
			}//&& freqn >= 3
			if ((int)(p3 + 0.5) == 0 && mode == 4 )
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C5I);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C5I);
			}//&& freqn >= 3
		}
	}

	/// Carrier phase data
	if (flag > 1 || flag < 0)
	{
		/// Rinex 2 version
		if (obsdata.m_ObsType->m_Version < 3)
		{
			l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1);
			pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1);
			if ((int)(l1 + 0.5) == 0 || pflag == 1)	return 0;
		}
		else /// Rinex 3 version
		{
			if ((int)(l1 + 0.5) == 0)
			{
				l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1P);
				pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1P);
				if ((int)(l1 + 0.5) != 0) channel.push_back("C1P");
			}
			if ((int)(l1 + 0.5) == 0)
			{
				l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1W);
				pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1W);
				if ((int)(l1 + 0.5) != 0) channel.push_back("C1W");
			}
			if ((int)(l1 + 0.5) == 0)
			{
				l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1C);
				pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1C);
				if ((int)(l1 + 0.5) != 0) channel.push_back("C1P");
			}
			
			if ((int)(l1 + 0.5) == 0)
			{
				l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1X);
				pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1X);
				if ((int)(l1 + 0.5) != 0) channel.push_back("C1X");
			}
			if ((int)(l1 + 0.5) == 0 || pflag == 1)	return 0;
			else {
				obstype.push_back("L1"); obstypeint.push_back(1);
			}
			if ((int)(l2 + 0.5) == 0)
			{
				l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L2P);
				pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L2P);
				if ((int)(l2 + 0.5) != 0) channel.push_back("C2P");
			}
			if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L2X);
				pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L2X);
				if ((int)(l2 + 0.5) != 0) channel.push_back("C2X");
			}
			if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L2L);
				pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L2L);
				if ((int)(l2 + 0.5) != 0) channel.push_back("C2L");
			}
			if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L2W);
				pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L2W);
				if ((int)(l2 + 0.5) != 0) channel.push_back("C2W");
			}
			if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L2C);
				pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L2C);
				if ((int)(l2 + 0.5) != 0) channel.push_back("C2P");
			}
			if ((int)(l2 + 0.5) == 0 && (mode == 0) || pflag2 == 1)return 0;//||mode==4
			else if ((int)(l2 + 0.5) == 0) {}
			else{ obstype_num +=1; obstype.push_back("L2"); obstypeint.push_back(2);
			}
			if ((int)(l3 + 0.5) == 0 && mode == 4 )
			{
				l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L5Q);
				pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L5Q);
				if ((int)(l3 + 0.5) != 0) channel.push_back("C5Q");
			}//&&freqn>=3
			if ((int)(l3 + 0.5) == 0 && mode == 4 )
			{
				l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L5X);
				pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L5X);
				if ((int)(l3 + 0.5) != 0) channel.push_back("C5X");
			}//&&freqn>=3
			if ((int)(l3 + 0.5) == 0 && mode == 4 )
			{
				l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L5I);
				pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L5I);
				if ((int)(l3 + 0.5) != 0) channel.push_back("C5I");
			}//&&freqn>=3
			if ((int)(l3 + 0.5) != 0 && pflag3 != 1 && (int)(p3 + 0.5) != 0 && f3 != 1) {
				obstype_num += 1; obstype.push_back("L5"); obstypeint.push_back(3); 
			}
			epoch_otype = channel;
		}
	}
	if (channel.size() != freqn) { return 0; }
	else
	{
		Prn_Otype.push_back(obstype); Prn_Otype_int.push_back(obstypeint); Prn_Otype_channel.push_back(channel);
		return 1;
	}
}

/// Get glonass data
int PPPObsData::GetGLOData(CodeDCB &dcbdata, ObsEpochData& obsdata, SatelliteData &epochdata, unsigned int freqn, int flag, int mode)
{
	fl = 0, f2 = 0, f3 = 0, f4 = 0, f5 = 0, pflag = 0, pflag2 = 0, pflag3 = 0, pflag4 = 0, pflag5 = 0;
	c1 = 0, p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0;

	d1 = d2 = d3 = 0; obstype_num = 1; vector<string> obstype; vector<int>obstypeint; vector<string> channel;
	obstype.push_back(to_string(epochdata.m_Prn));
	if (flag < 2)
	{
		/// Rinex 2.0 version
		if (obsdata.m_ObsType->m_Version < 3)
		{
			p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, P1);
			fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, P1);

			if ((int)(p1 + 0.5) == 0)
			{
				c1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1);
				fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1);
				if ((int)(c1 + 0.5) == 0) return 0;
				else
				{
					double dcb = dcbdata.Get_DCB_P1C1(epochdata.m_Prn);
					p1 = c1 + dcb*LIGHTSPEED*1.0e-9;
				}
			}

			if (fl == 1) return 0;

			/// pseudorange quality check
			c1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1);
			if (c1 != 0 && p1 != 0)
			{
				if (PseudoRange_C1P1_Detect(c1, p1)) return 0;
			}
		}
		else /// Rinex 3.0 version
		{
			p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1P);//G1
			fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1P);
			p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C2P);//G2
			f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C2P);

			
			
			if ((int)(p1 + 0.5) == 0)
			{
				c1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1C);
				fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1C);

				if ((int)(c1 + 0.5) == 0) return 0;
				else
				{
					double dcb = dcbdata.Get_DCB_P1C1(epochdata.m_Prn);
					p1 = c1 + dcb * LIGHTSPEED * 1.0e-9;
				}
			}
			if (fl == 1) return 0;
			/// pseudorange quality check
			c1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1C);
			if (c1 != 0 && p1 != 0)
			{
				if (PseudoRange_C1P1_Detect(c1, p1)) return 0;
			}
			if ((int)(p2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C2C);
				f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C2C);
			}
			if ((p2 == 0 || f2 == 1) && (mode == 0 || mode == 4)) return 0;
			if ((int)(p3 + 0.5) == 0 && mode == 4&&freqn>=3) {
			p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C4A);//G1a
			f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C4A);
			}
			if ((int)(p3 + 0.5) == 0 && mode == 4 && freqn >= 3)
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C4B);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C4B);
			}
			if ((int)(p3 + 0.5) == 0 && mode == 4 && freqn >= 3)
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C4X);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C4X);
			}
			if ((int)(p4 + 0.5) == 0 && mode == 4 && freqn >= 4) {
				p4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C6A);//G2a
				f4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C6A);
			}
			if ((int)(p4 + 0.5) == 0 && mode == 4 && freqn >= 4)
			{
				p4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C6B);
				f4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C6B);
			}
			if ((int)(p4 + 0.5) == 0 && mode == 4 && freqn >= 4)
			{
				p4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C6X);
				f4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C6X);
			}
			if ((int)(p5 + 0.5) == 0 && mode == 4 && freqn == 5) {
				p5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C3I);//G3
				f5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C3I);
			}
			if ((int)(p5 + 0.5) == 0 && mode == 4 && freqn == 5)
			{
				p5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C3Q);
				f5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C3Q);
			}
			if ((int)(p5 + 0.5) == 0 && mode == 4 && freqn == 5)
			{
				p5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C3X);
				f5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C3X);
			}

		}
	}
	/// Carrier phase data
	if (flag > 1 || flag < 0)
	{
		/// Rinex 2 version
		if (obsdata.m_ObsType->m_Version < 3)
		{
			l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1);
			pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1);
			if ((int)(l1 + 0.5) == 0 || pflag == 1)	return 0;
		}
		else /// Rinex 3 version
		{
			
			if ((int)(l1 + 0.5) == 0) {
				l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1P);
				pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1P);
				if ((int)(l1 + 0.5) != 0) channel.push_back("C1P");
			}
			if ((int)(l1 + 0.5) == 0)
			{
				l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1C);
				pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1C);
				if ((int)(l1 + 0.5) != 0) channel.push_back("C1C");
			}
			if ((int)(l1 + 0.5) == 0 || pflag == 1)	return 0;
			else
			{
				obstype.push_back("G1"); obstypeint.push_back(1);;
			}
			if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4)) {
				l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L2P);
				pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L2P);
				if ((int)(l2 + 0.5) != 0) channel.push_back("C2P");
			}
			if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L2C);
				pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L2C);
				if ((int)(l2 + 0.5) != 0) channel.push_back("C2C");
			}
			if ((int)(l2 + 0.5) == 0 || pflag2 == 1&&(mode==0 || mode==4)) return 0;
			else { obstype_num += 1; obstype.push_back("G2"); obstypeint.push_back(2);
			}
			if ((int)(l3 + 0.5) == 0 && mode == 4&&freqn>=3) {
				l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L4A);
				pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L4A);
				if ((int)(l3 + 0.5) != 0) channel.push_back("C4A");
			}
			if ((int)(l3 + 0.5) == 0 && mode == 4 &&freqn>= 3)
			{
				l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L4B);
				pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L4B);
				if ((int)(l3 + 0.5) != 0) channel.push_back("C4B");
			}
			if ((int)(l3 + 0.5) == 0 && mode == 4 && freqn >= 3)
			{
				l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L4X);
				pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L4X);
				if ((int)(l3 + 0.5) != 0) channel.push_back("C4X");
			}
			if (((int)(l3 + 0.5) == 0 || (int)(p3 + 0.5) == 0 || pflag3 == 1 || f3 == 1) && mode == 4);
			else { obstype_num += 1; obstype.push_back("G1a"); obstypeint.push_back(3);
			}
			if ((int)(l4 + 0.5) == 0 && mode == 4&&freqn >= 4) {
				l4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L6A);
				pflag4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L6A);
				if ((int)(l4 + 0.5) != 0) channel.push_back("C6A");
			}
			if ((int)(l4 + 0.5) == 0 && mode == 4 && freqn >= 4)
			{
				l4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L6B);
				pflag4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L6B);
				if ((int)(l4 + 0.5) != 0) channel.push_back("C6B");
			}
			if ((int)(l4 + 0.5) == 0 && mode == 4 && freqn >= 4)
			{
				l4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L6X);
				pflag4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L6X);
				if ((int)(l4 + 0.5) != 0) channel.push_back("C6X");
			}
			if (((int)(l4 + 0.5) == 0 || (int)(p4 + 0.5) == 0 || pflag4 == 1 || f4 == 1) && mode == 4);
			else { obstype_num += 1; obstype.push_back("G2a"); obstypeint.push_back(4); 
			}
			if ((int)(l5 + 0.5) == 0 && mode == 4&&freqn==5) {
				l5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L3I);
				pflag5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L3I);
				if ((int)(l5 + 0.5) != 0) channel.push_back("C3I");
			}
			if ((int)(l5 + 0.5) == 0 && mode == 4 && freqn == 5)
			{
				l5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L3Q);
				pflag5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L3Q);
				if ((int)(l5 + 0.5) != 0) channel.push_back("C3Q");
			}
			if ((int)(l5 + 0.5) == 0 && mode == 4 && freqn == 5)
			{
				l5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L3X);
				pflag5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L3X);
				if ((int)(l5 + 0.5) != 0) channel.push_back("C3X");
			}
			if (((int)(l5 + 0.5) == 0 || (int)(p5 + 0.5) == 0 || pflag5 == 1 ||f5==1)&& mode == 4);
			else { obstype_num += 1; obstype.push_back("G3"); obstypeint.push_back(5);
			}
			epoch_otype = channel; 
		}
	}
	if (channel.size() != freqn) { return 0; }
	else
	{
		Prn_Otype.push_back(obstype); Prn_Otype_int.push_back(obstypeint); Prn_Otype_channel.push_back(channel);
		return 1;
	}
}

/// Get galileo data
int PPPObsData::GetGALData(CodeDCB &dcbdata, ObsEpochData& obsdata, SatelliteData &epochdata, unsigned int freqn, int flag, int mode)
{
	fl = 0, f2 = 0, f3 = 0, f4 = 0, f5 = 0, pflag = 0, pflag2 = 0, pflag3 = 0, pflag4 = 0, pflag5 = 0;
	c1 = 0, p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0;

	d1 = d2 = d3 = 0; obstype_num = 1; vector<string> obstype; vector<int>obstypeint; vector<string> channel;
	obstype.push_back(to_string(epochdata.m_Prn));
	if (flag < 2)
	{
		/// Rinex 2.0 version
		if (obsdata.m_ObsType->m_Version < 3)
		{
			return 0;
		}
		else /// Rinex 3.0 version
		{
			p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1X);
			fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1X);
			p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C5I);
			f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C5I);
			p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C7I);
			f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C7I);
			p4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C8I);
			f4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C8I);
			p5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C6X);
			f5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C6X);
			
			if ((int)(p1 + 0.5) == 0)
			{
				p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1C);
				fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1C);
			}
			if (fl == 1|| (int)(p1 + 0.5) == 0) return 0;
			if ((int)(p2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C5Q);
				f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C5Q);
			}
			if ((int)(p2 + 0.5) == 0 && (mode == 0 || mode == 4))
			{
				p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C5X);
				f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C5X);
			}
			if (((int)(p2 + 0.5) == 0 || f2 == 1) && ( mode == 0||mode ==4))return 0;
			if ((int)(p3 + 0.5) == 0 && mode == 4)
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C7Q);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C7Q);
			}
			if ((int)(p3 + 0.5) == 0 && mode == 4)
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C7X);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C7X);
			}
			if ((int)(p4 + 0.5) == 0 && mode == 4)
			{
				p4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C8Q);
				f4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C8Q);
			}
			if ((int)(p4 + 0.5) == 0 && mode == 4)
			{
				p4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C8X);
				f4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C8X);
			}
			if ((int)(p5 + 0.5) == 0 && mode == 4)
			{
				p5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C6C);
				f5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C6C);
			}
			if ((int)(p5 + 0.5) == 0 && mode == 4)
			{
				p5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C6B);
				f5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C6B);
			}
		}
	}

	if (flag > 1 || flag < 0)
	{
		if ((int)(l1 + 0.5) == 0) {
			l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1X);
			pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1X);
			if ((int)(l1 + 0.5) != 0) channel.push_back("C1X");
		}
		if ((int)(l1 + 0.5) == 0)
		{
			l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1C);
			pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1C);
			if ((int)(l1 + 0.5) != 0) channel.push_back("C1C");
		}
		if (((int)(l1 + 0.5) == 0 || pflag == 1)) return 0;
		else { obstype.push_back("E1"); obstypeint.push_back(1); }
		if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4)) {
			l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L5I);
			pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L5I);
			if ((int)(l2 + 0.5) != 0) channel.push_back("C5I");
		}
		if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4))
		{
			l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L5Q);
			pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L5Q);
			if ((int)(l2 + 0.5) != 0) channel.push_back("C5Q");
		}
		if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4))
		{
			l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L5X);
			pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L5X);
			if ((int)(l2 + 0.5) != 0) channel.push_back("C5X");
		}
		if (((int)(l2 + 0.5) == 0 || pflag2 == 1) && (mode == 0 || mode == 4)) return 0;
		else if ((int)(l2 + 0.5) == 0) {}
		else{ obstype_num += 1; obstype.push_back("E5a"); obstypeint.push_back(2); }
		if ((int)(l3 + 0.5) == 0 && mode == 4&&freqn>=3) {
		l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L7I);
		pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L7I);
		if ((int)(l3 + 0.5) != 0) channel.push_back("C7I");
		}
		if ((int)(l3 + 0.5) == 0 && mode == 4 && freqn >= 3)
		{
			l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L7Q);
			pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L7Q);
			if ((int)(l3 + 0.5) != 0) channel.push_back("C7Q");
		}
		if ((int)(l3 + 0.5) == 0 && mode == 4 && freqn >= 3)
		{
			l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L7X);
			pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L7X);
			if ((int)(l3 + 0.5) != 0) channel.push_back("C7X");
		}
		if (((int)(l3 + 0.5) == 0 || pflag3 == 1|| (int)(p3 + 0.5) == 0 || f3 == 1) && (mode == 0 || mode == 4)) ;
		else { obstype_num += 1; obstype.push_back("E5b"); obstypeint.push_back(3);
		}
		if ((int)(l4 + 0.5) == 0 && mode == 4 && freqn >= 4) {
			l4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L8I);
			pflag4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L8I);
			if ((int)(l4 + 0.5) != 0) channel.push_back("C8I");
		}
		if ((int)(l4 + 0.5) == 0 && mode == 4 && freqn >= 4)
		{
			l4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L8Q);
			pflag4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L8Q);
			if ((int)(l4 + 0.5) != 0) channel.push_back("C8Q");
		}
		if ((int)(l4 + 0.5) == 0 && mode == 4 && freqn >= 4)
		{
			l4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L8X);
			pflag4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L8X);
			if ((int)(l4 + 0.5) != 0) channel.push_back("C8X");
		}
		if (((int)(l4 + 0.5) == 0 || pflag4 == 1 || (int)(p4 + 0.5) == 0 || f4 == 1) && (mode == 0 || mode == 4));
		else { obstype_num += 1; obstype.push_back("E5(ab)"); obstypeint.push_back(4); 
		}
		if ((int)(l5 + 0.5) == 0 && mode == 4 && freqn >= 5) {
			l5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L6X);
			pflag5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L6X);
			if ((int)(l5 + 0.5) != 0) channel.push_back("C6X");
		}
		if ((int)(l5 + 0.5) == 0 && mode == 4 && freqn >= 5)
		{
			l5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L6C);
			pflag5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L6C);
			if ((int)(l5 + 0.5) != 0) channel.push_back("C6C");
		}
		if ((int)(l5 + 0.5) == 0 && mode == 4 && freqn >= 5)
		{
			l5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L6B);
			pflag5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L6B);
			if ((int)(l5 + 0.5) != 0) channel.push_back("C6B");
		}
		if (((int)(l5 + 0.5) == 0 || pflag5 == 1 || (int)(p5 + 0.5) == 0 || f5 == 1) && (mode == 0 || mode == 4));
		else { 
			obstype_num += 1; obstype.push_back("E6"); obstypeint.push_back(5); 
		}
		epoch_otype = channel; 
	}
	if (channel.size() != freqn) { return 0; }
	else
	{
		Prn_Otype.push_back(obstype); Prn_Otype_int.push_back(obstypeint); Prn_Otype_channel.push_back(channel);
		return 1;
	}
}

/// Get bds data
int PPPObsData::GetBDSData(CodeDCB& dcbdata, ObsEpochData& obsdata, SatelliteData& epochdata, unsigned int freqn, int flag, int prn, int mode)
{
	fl = 0, f2 = 0, f3 = 0, f4 = 0, f5 = 0, pflag = 0, pflag2 = 0, pflag3 = 0, pflag4 = 0, pflag5 = 0;
	c1 = 0, p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0;

	d1 = d2 = d3 = 0; obstype_num = 1; vector<string> obstype; vector<int>obstypeint; vector<string> channel;
	obstype.push_back(to_string(epochdata.m_Prn)); int a;

	if (prn > 360)
	{
		return 0;
	}
	if (flag < 2)
	{
		/// Rinex 2.0 version
		if (obsdata.m_ObsType->m_Version < 3)
		{
			//
			return 0;
		}
		else /// Rinex 3.0 version
		{


			if ((int)(p1 + 0.5) == 0)
			{
				p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C2I);
				fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C2I);
			}
			/*if ((int)(p1 + 0.5) == 0)
			{
				p1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C2X);
				fl = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C2X);
			}*/
			if ((int)(p1 + 0.5) == 0 || fl == 1)return 0;
			if ((int)(p2 + 0.5) == 0)
			{
				p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C6I);
				f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C6I);
			}
			/*if ((int)(p2 + 0.5) == 0 )
			{
				p2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C6X);
				f2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C6X);
			}*/
			if (((int)(p2 + 0.5) == 0 || f2 == 1) && (mode == 0 || mode == 4));//return 0;
			if ((int)(p3 + 0.5) == 0 && mode == 4)
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C7I);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C7I);
			}
			if ((int)(p3 + 0.5) == 0 && mode == 4)
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C7Q);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C7Q);
			}
			if ((int)(p3 + 0.5) == 0 && mode == 4)
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C7X);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C7X);
			}
			if ((int)(p3 + 0.5) == 0 && mode == 4)
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C7D);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C7D);
			}
			if ((int)(p3 + 0.5) == 0 && mode == 4)
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C7P);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C7P);
			}
			if ((int)(p3 + 0.5) == 0 && mode == 4)
			{
				p3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C7Z);
				f3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C7Z);
			}
			if ((int)(p4 + 0.5) == 0 && mode == 4)
			{
				p4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1D);
				f4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1D);
			}
			if ((int)(p4 + 0.5) == 0 && mode == 4)
			{
				p4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1P);
				f4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1P);
			}
			if ((int)(p4 + 0.5) == 0 && mode == 4)
			{
				p4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C1X);
				f4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C1X);
			}
			if ((int)(p5 + 0.5) == 0 && mode == 4)
			{
				p5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C5D);
				f5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C5D);
			}
			if ((int)(p5 + 0.5) == 0 && mode == 4)
			{
				p5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C5P);
				f5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C5P);
			}
			if ((int)(p5 + 0.5) == 0 && mode == 4)
			{
				p5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, C5X);
				f5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, C5X);
			}
		}
	}

	if (flag > 1 || flag < 0)
	{

		if ((int)(l1 + 0.5) == 0) {
			l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L2I);
			pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L2I);
			if ((int)(l1 + 0.5) != 0) channel.push_back("C2I");
		}
		/*if ((int)(l1 + 0.5) == 0)
		{
			l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L2Q);
			pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L2Q);
			if ((int)(l1 + 0.5) != 0) channel.push_back("C2Q");
		}
		if ((int)(l1 + 0.5) == 0)
		{
			l1 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L2X);
			pflag = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L2X);
			if ((int)(l1 + 0.5) != 0) channel.push_back("C2X");
		}*/
		if ((int)(l1 + 0.5) == 0 || pflag == 1)return 0;
		else { obstype.push_back("B1"); obstypeint.push_back(1); }
		if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4))
		{
			l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L6I);
			pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L6I);
			if ((int)(l2 + 0.5) != 0) channel.push_back("C6I");
		}
		/*if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4))
		{
			l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L6Q);
			pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L6Q);
			if ((int)(l2 + 0.5) != 0) channel.push_back("C6Q");
		}
		if ((int)(l2 + 0.5) == 0 && (mode == 0 || mode == 4))
		{
			l2 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L6X);
			pflag2 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L6X);
			if ((int)(l2 + 0.5) != 0) channel.push_back("C6X");
		}*/
		if ((l2 == 0 || pflag == 1) && (mode == 0 || mode == 4));//return 0;
		else if (l2 == 0) {}
		else{
			obstype_num = 2; obstype.push_back("B3"); obstypeint.push_back(2);
		}
		if ((int)(l3 + 0.5) == 0 && mode == 4 && freqn >= 3)
		{
			l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L7I);
			pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L7I);
			if ((int)(l3 + 0.5) != 0) channel.push_back("C7I");
		}
		if ((int)(l3 + 0.5) == 0 && mode == 4 && freqn >= 3)
		{
			l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L7Q);
			pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L7Q);
			if ((int)(l3 + 0.5) != 0) channel.push_back("C7Q");
		}
		if ((int)(l3 + 0.5) == 0 && mode == 4 && freqn >= 3)
		{
			l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L7X);
			pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L7X);
			if ((int)(l3 + 0.5) != 0) channel.push_back("C7X");
		}
		if ((int)(l3 + 0.5) == 0 && mode == 4 && freqn >= 3)
		{
			l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L7D);
			pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L7D);
			if ((int)(l3 + 0.5) != 0) channel.push_back("C9D");
		}
		if ((int)(l3 + 0.5) == 0 && mode == 4 && freqn >= 3)
		{
			l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L7P);
			pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L7P);
			if ((int)(l3 + 0.5) != 0) channel.push_back("C9P");
		}
		if ((int)(l3 + 0.5) == 0 && mode == 4 && freqn >= 3)
		{
			l3 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L7Z);
			pflag3 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L7Z);
			if ((int)(l3 + 0.5) != 0) channel.push_back("C9Z");
		}
		if ((l3 == 0 || p3 == 0 || f3 == 1 || pflag3 == 1) && mode == 4);
		else {
			obstype_num += 1; obstype.push_back("B2"); obstypeint.push_back(3); //包含  B3:b2b
		}
		if ((int)(l4 + 0.5) == 0 && mode == 4 && freqn >= 4)
		{
			l4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1P);
			pflag4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1P);
			if ((int)(l4 + 0.5) != 0) channel.push_back("C1P");
		}
		if ((int)(l4 + 0.5) == 0 && mode == 4 && freqn >= 4)
		{
			l4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1D);
			pflag4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1D);
			if ((int)(l4 + 0.5) != 0) channel.push_back("C1D");
		}
		if ((int)(l4 + 0.5) == 0 && mode == 4 && freqn >= 4)
		{
			l4 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L1X);
			pflag4 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L1X);
			if ((int)(l4 + 0.5) != 0) channel.push_back("C1X");
		}
		if ((l4 == 0 || p4 == 0 || f4 == 1 || pflag4 == 1) && mode == 4);
		else {
			obstype_num += 1; obstype.push_back("B1C"); obstypeint.push_back(4);
		}
		if ((int)(l5 + 0.5) == 0 && mode == 4 && freqn >= 5)
		{
			l5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L5P);
			pflag5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L5P);
			if ((int)(l5 + 0.5) != 0) channel.push_back("C5P");
		}
		if ((int)(l5 + 0.5) == 0 && mode == 4 && freqn >= 5)
		{
			l5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L5D);
			pflag5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L5D);
			if ((int)(l5 + 0.5) != 0) channel.push_back("C5D");
		}
		if ((int)(l5 + 0.5) == 0 && mode == 4 && freqn >= 5)
		{
			l5 = GET_OBS_DATA(epochdata, obsdata.m_ObsType, L5X);
			pflag5 = GET_OBS_FLAG(epochdata, obsdata.m_ObsType, L5X);
			if ((int)(l5 + 0.5) != 0) channel.push_back("C5X");
		}
		if ((l5 == 0 || p5 == 0 || f5 == 1 || pflag5 == 1) && mode == 4);
		else {
			obstype_num += 1; obstype.push_back("B2a"); obstypeint.push_back(5);
		}
		epoch_otype = channel;
	}
	if (channel.size() != freqn) { return 0; }
	else 
	{ 
		Prn_Otype.push_back(obstype); Prn_Otype_int.push_back(obstypeint); Prn_Otype_channel.push_back(channel);
		return 1; 
	}

}

/// PPP Rank Detection  0:rank defect  1:full rank
int RankDetect(ppp_option_t& popt, unsigned int satn)
{
	unsigned MathModel = 0, FreqNum = 0, TropNum = 0;
	unsigned int para_n = 0, eq_n = 0;   // Parameter number, Equation number

	MathModel = popt.math_model;
	FreqNum = popt.freqn;

	if (popt.process_mode == MODE_SPP)
	{
		para_n = 3 * popt.estpara.xyz + popt.estpara.rclk;

		eq_n = satn;

		if (eq_n < para_n)
			return 1;
		else
			return 0;
	}

	// troposphere parameter numbers in different trop model
	switch (popt.trop_model)
	{
		case TROP_ESTIMATE:{TropNum = 1; break; }
		case TROP_GRADIENT:{TropNum = 3; break; }
		case TROP_GRID:{TropNum = 0; break; }
		default:{TropNum = 1; break; }
	}

	if (MathModel == MODE_UNCOMBINE)
	{
		if (FreqNum == 1)
		{
			para_n = 5;
			eq_n = 2 * satn;
		}
	}

	if (satn < 5)
		return 1;

	return 0;
}

void DeleteSatData(PPPObsData& input, unsigned int index)
{
	if (input.m_SatCount == 0) return;
	
	input.m_Valid_flag[index] = input.m_Valid_flag[input.m_SatCount - 1];
	input.m_Prn[index] = input.m_Prn[input.m_SatCount - 1];
	input.m_Pr1[index] = input.m_Pr1[input.m_SatCount - 1];
	input.m_Pr2[index] = input.m_Pr2[input.m_SatCount - 1];
	input.m_Pr3[index] = input.m_Pr3[input.m_SatCount - 1];
	
	input.m_L1[index] = input.m_L1[input.m_SatCount - 1];
	input.m_L2[index] = input.m_L2[input.m_SatCount - 1];
	input.m_L3[index] = input.m_L3[input.m_SatCount - 1];
	
	input.m_Prn[input.m_SatCount - 1] = 0;
	
	//input.Prn_Otype.erase(input.Prn_Otype.begin()+ index);
	--input.m_SatCount;
}
