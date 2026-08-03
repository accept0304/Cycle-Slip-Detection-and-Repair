#include <Windows.h>
#include <string.h>

#include "SUPREME_PPP.h"
#include "SUPREME_ObsRNX.h"
#include "SUPREME_NavRNX.h"
#include "SUPREME_ClkRNX.h"
#include "SUPREME_GenData.h"
#include "SUPREME_ParaEst.h"
#include "SUPREME_CommonFunction.h"
#include "SUPREME_QC.h"
#include "SUPREME_Output.h"
#include "SUPREME_GraphPlot.h"
#include "SUPREME_STACRD.h"
#include "SUPREME_DownLoad.h"
#include "SUPREME_Ionosphere.h"
extern vector<double> SimulatedCbias;

/// Load gnss data: Observation, Navigation, Precise Ephemeris, Precise Clock, General data...
int LoadGnssData(ppp_option_t &popt, ObsData *obsdata, PreciseData *predata, ClockData *clkdata, GenData &gendata)
{
	// Get precise ephemeris
	if (predata)
		predata->Read_PreciseFile(popt.IOfile.rinex_f.prephf, 1);

	// Get precise clock
	if (clkdata)
		clkdata->Read_ClockFile(popt.IOfile.rinex_f.preclkf);

	// Get Navigation
	if ((!popt.use_orbclk) || (popt.System&PRO_SYS_GLO))
		predata->Read_PreciseFile(popt.IOfile.rinex_f.navf, 0);

	// Get gnss general data
	gendata.LoadGenData(popt);

	return  1;
}

/// Post PPP by least square filter
int SUPREME_PPP_LSF_Post(ppp_config &cfg)
{
	ppp_option_t *popt;
	popt = cfg.popt;

	// Get station number by station list
	unsigned int StaN = 0;
	char StationList[MAXSTATION][5] = { 0 };
	if (cfg.batchmode) StaN = GetStationList(popt->stationlist, StationList);
	else StaN = 1;

	// Station-wised process
	for (unsigned int i_sta = 0; i_sta < StaN; ++i_sta)            
	{
		// Process time
		gnsstime process_time(popt->starttime);
		popt->CurrentEpoch = process_time;

		char ACcentre[5] = { 0 };  // Analysis Center name
		char CENTER[5] = { 0 };
		sprintf(ACcentre, "%s", popt->CentreName);
		if (strstr(ACcentre,"com")!=NULL)
		{
			sprintf(CENTER, "%s", "cod");
		}
		else if (strstr(ACcentre, "gbm") != NULL)
		{
			sprintf(CENTER, "%s", "gfz");
		}
		char StaName[5] = { 0 };   // Station name
		sprintf(StaName, "%s", popt->StaName);

		unsigned int DayN = 0;     // Get process day number
		if (cfg.batchmode)
			DayN = GetDayNumber(popt->starttime, popt->endtime);
		else DayN = 1;

		LeastSquarePPP LsqPPP(*popt);

		ObsEpochData *ptr = NULL;
		ObsEpochData *zhy = NULL;
		if (cfg.batchmode)
		{
			memset(StaName, 0, sizeof(char)* 5);
			strcpy(StaName, StationList[i_sta]);
			memset(popt->StaName, 0, sizeof(char)* 50);
			sprintf(popt->StaName, "%s", StaName);
		}

		//% Get station coordinates ------------------------
		StationCrd sta;char crdpath[PATH_LENGTH] = { 0 };
		sprintf(crdpath, "%s%s", popt->ObsPath,"STATION.CRD");
		Lower2Cap(StaName, sta.m_StaName);
		if (sta.GetStaSnx(popt->IOfile.gen_f.crd_fl) == 0) 
		{
			sta.GetStaCrd("F:\\STATION.CRD");
			sta.Display();
		}
		popt->StaPos_CRD[0] = sta.StaCoord._X;
		popt->StaPos_CRD[1] = sta.StaCoord._Y;
		popt->StaPos_CRD[2] = sta.StaCoord._Z;

		if (Norm(popt->StaPos_CRD, 3) == 0 && popt->estpara.xyz == 0)
		{
			printf("Warning: No Station Coordinates to Fix...\n");
			continue;
		}

		printf("******* PROCESS START *******\n");

		//随机历元
		double percentvalue = 0.0;//周跳占比
		if (popt->rateorppp == 1)
		{
			percentvalue = popt->pppsrate;//周跳占比
		}
		int cycletime = 0, cyc_ount = 0;//
		int cycsource = 6000;
		cyc_ount = percentvalue * cycsource;
		int tm = 0;
		set<int> cyctime_set; popt->simnumc = cyc_ount;
		do {
			cycletime = rand() % cycsource + 1;
			cyctime_set.insert(cycletime);
		} while (cyctime_set.size() < cyc_ount);
		if (cyc_ount != 0)
		{
			for (set<int>::iterator it = cyctime_set.begin(); it != cyctime_set.end(); it++)
			{
				popt->simulcycletime[tm] = (*it);
				tm++;
			}
		}

		string stationname1;
		stationname1.assign(StationList[i_sta]);
		// day-wised process
		for (unsigned int i_day = 0; i_day < DayN; ++i_day)
		{
			int year = 0, doy = 0;
			year = process_time.m_Year;
			doy = process_time.m_Doy;
			string stationname2, stationname3, stationname4;
			stationname2 = '.'; stationname3 = to_string(doy) + '0'; stationname4 = to_string(year - 2000) + 'o';
			stationname1 += stationname3;
			stationname2 += stationname4;
			stationname1 += stationname2;

			// Observation path
			char obs_path[PATH_LENGTH] = { 0 };
			char obs_dir[PATH_LENGTH] = { 0 };
			char obs_path_tmp[PATH_LENGTH] = { 0 };
			char obs_BRDM[PATH_LENGTH] = { 0 };
 			if (cfg.batchmode)
			{
				sprintf(obs_dir, "%s%4d/%03d/", popt->ObsPath, year, doy);
				sprintf(obs_path, "%s%s%03d0.%02do", popt->ObsPath, StaName, doy, year % 100);
				sprintf(obs_BRDM, "%s%s%03d0.%02dp", popt->ObsPath, "brdm", doy, year % 100);
				memset(popt->IOfile.rinex_f.navf, 0, sizeof(char)*PATH_LENGTH);
				strcpy(popt->IOfile.rinex_f.obsf, obs_path);
				strcpy(popt->IOfile.rinex_f.navf, obs_BRDM);
			}
			else
				sprintf(obs_path, "%s", popt->IOfile.rinex_f.obsf);

			// Precise orbit and clock path (three days)
			char eph_path1[PATH_LENGTH] = { 0 }, eph_path2[PATH_LENGTH] = { 0 }, eph_path3[PATH_LENGTH] = { 0 };
			char clk_path1[PATH_LENGTH] = { 0 }, clk_path2[PATH_LENGTH] = { 0 }, clk_path3[PATH_LENGTH] = { 0 };
			char eop_path[PATH_LENGTH] = { 0 };
			gnsstime t1, t2;
			t1 = process_time + (-86400);
			t2 = process_time + 86400;

			if (cfg.batchmode)    
			{
				sprintf(eop_path, "%s%s%4d7.erp", cfg.popt->PrePath, CENTER, t1.m_Week);
				sprintf(eph_path1, "%s%s%4d%d.sp3", cfg.popt->PrePath, ACcentre, t1.m_Week, t1.m_Dow);
				sprintf(eph_path2, "%s%s%4d%d.sp3", cfg.popt->PrePath, ACcentre, process_time.m_Week, process_time.m_Dow);
				sprintf(eph_path3, "%s%s%4d%d.sp3", cfg.popt->PrePath, ACcentre, t2.m_Week, t2.m_Dow);
				if (!FileExist(eph_path1)) Download_EphClk(cfg.popt->PrePath, ACcentre, t1, EPHCLK);
				if (!FileExist(eph_path2)) Download_EphClk(cfg.popt->PrePath, ACcentre, process_time, EPHCLK);
				if (!FileExist(eph_path3)) Download_EphClk(cfg.popt->PrePath, ACcentre, t2, EPHCLK);

				sprintf(clk_path1, "%s%s%4d%d.clk", cfg.popt->PrePath, ACcentre, t1.m_Week, t1.m_Dow);
				sprintf(clk_path2, "%s%s%4d%d.clk", cfg.popt->PrePath, ACcentre, process_time.m_Week, process_time.m_Dow);
				sprintf(clk_path3, "%s%s%4d%d.clk", cfg.popt->PrePath, ACcentre, t2.m_Week, t2.m_Dow);
				if (!FileExist(clk_path1)) Download_EphClk(cfg.popt->PrePath, ACcentre, t1, CLK);
				if (!FileExist(clk_path2)) Download_EphClk(cfg.popt->PrePath, ACcentre, process_time, CLK);
				if (!FileExist(clk_path3)) Download_EphClk(cfg.popt->PrePath, ACcentre, t2, CLK);
			}
			else
			{
				strcpy(eph_path2, popt->IOfile.rinex_f.prephf);
				strcpy(clk_path2, popt->IOfile.rinex_f.preclkf);
			}

			if (cfg.batchmode)
			{
				memset(popt->IOfile.rinex_f.prephf, 0, sizeof(char)*PATH_LENGTH);
				memset(popt->IOfile.rinex_f.preclkf, 0, sizeof(char)*PATH_LENGTH);
				strcpy(popt->IOfile.rinex_f.prephf, eph_path2);
				strcpy(popt->IOfile.rinex_f.preclkf, clk_path2);
				strcpy(popt->IOfile.gen_f.erp_fl, eop_path);
			}

			if (!FileExist(obs_path))
			{
				//ObsNav_download(StaName, obs_path2, pro_time, MGEX, 'o');
				printf("Warrning: Observation %s is not exsiting...\nContinue...\n", obs_path);
				continue;
			}

			if (!FileExist(eph_path2) || !FileExist(clk_path2))
			{
				printf("Warrning: Precise orbit or clock product are not exsiting...\n");
				continue;
			}

			if (!cfg.popt->ContinueMode)	LsqPPP.Initialize();

			PreciseData  predata;                     // Precise orbit product
			ClockData    clkdata;                     // Precise clock product

			predata.Read_PreciseFile(eph_path2, 1);
			clkdata.Read_ClockFile(clk_path2);

			if ((popt->use_orbclk) || (popt->System&PRO_SYS_GLO)) // Read navigation
				predata.Read_PreciseFile(popt->IOfile.rinex_f.navf, 0);

			if (cfg.batchmode)
			{
				predata.Read_PreciseFile(eph_path3, 1);
				clkdata.Read_ClockFile(clk_path3);
			}

			ObsData obsdata;                          // Observation data

			// Get observation header information
			obsdata.m_Header.ReadObsFileHeader(obs_path);
			if (obsdata.m_Header.m_AntennaType[0] != ' ')   // Antenna type
			{
				memset(popt->antenna_type, 0, sizeof(char)* 21);
				strcpy(popt->antenna_type, obsdata.m_Header.m_AntennaType);
			}
			if (obsdata.m_Header.m_ReceiverType[0] != ' ')  // Receiver type
			{
				memset(popt->receiver_type, 0, sizeof(char)* 21);
				strcpy(popt->receiver_type, obsdata.m_Header.m_ReceiverType);
			}
			if (((i_day == 0) && cfg.popt->ContinueMode) || (!cfg.popt->ContinueMode))
			{
				popt->ProXYZ[0] = obsdata.m_Header.m_StaPos_x;
				popt->ProXYZ[1] = obsdata.m_Header.m_StaPos_y;
				popt->ProXYZ[2] = obsdata.m_Header.m_StaPos_z;
				popt->antenna_offset[0] = obsdata.m_Header.m_AntPos_N;
				popt->antenna_offset[1] = obsdata.m_Header.m_AntPos_E;
				popt->antenna_offset[2] = obsdata.m_Header.m_AntPos_H;
			}

			popt->CurrentEpoch = obsdata.m_Header.m_StartEpoch;

			GenData gendata;   // General information
			
			

			sprintf(popt->IOfile.gen_f.dcb_mgex, "%sCAS0MGXRAP_%d%03d0000_01D_01D_DCB.BSX", cfg.popt->PrePath, process_time.m_Year,process_time.m_Doy);
			sprintf(popt->IOfile.gen_f.ion_fl, "%s%s%03d0.%02di", cfg.popt->PrePath,"codg", doy, year % 100);

			sprintf(popt->IOfile.gen_f.ifcb_f, "%s%03d.ifcb", cfg.popt->PrePath, doy);

			if(!(gendata.LoadGenData(*popt)))return 0;  //读取  DCB   PCO and PCV  ERP  BLQ  *GIM*


			obsdata.ReadObsRinex(obs_path);

			LsqPPP._OutputFile = 1;
			if (popt->mode_filter == MODE_FILTER_FORWARD)         // Forward Filting
				ptr = obsdata.m_FirstEpoch;
			else if (popt->mode_filter == MODE_FILTER_BACKWARD|| popt->mode_filter == MODE_FILTER_COMBINE)   // Backward Filting
				ptr = obsdata.m_LastEpoch;
			else
			{
				printf("Warrning: Unrecognized filting mode.\n");
				return 0;
			}

			vector<zhydata>prepppdata;
			vector<zhyinfo>presatinfo;
			IonexModel iondata;
			deque<vector<int>> compareprn;
			deque<vector<double>> iongim123;
		
			ofstream rinexo, rinexopre, rineosimul;
			if (percentvalue != 0.0)
			{
				string rinexocom, rinexoprecom, rineosimulcom;
				string sbase, s1, s2, s3; sbase = "C:\\Users\\15149\\Desktop\\shiy";
				rinexocom = "\\stations\\"; rinexocom = sbase + rinexocom;
				rinexoprecom = "\\stations\\PRE_"; rinexoprecom = sbase + rinexoprecom;
				rineosimulcom = "\\stations\\SIM_"; rineosimulcom = sbase + rineosimulcom;
				rinexocom += stationname1; rinexoprecom += stationname1; rineosimulcom += stationname1;
				rinexo.open(rinexocom, ios::out); rinexopre.open(rinexoprecom, ios::out); rineosimul.open(rineosimulcom, ios::out);
			}

			int kk = 0; int zhy_flag=0;
			while (ptr)
			{
				gnsstime obst = ptr->m_EpochTime;
				
				// dt of current and previous epoch
				if (popt->dt == 0)
					popt->dt = popt->interval;
				else
					popt->dt = obst - popt->CurrentEpoch;
				popt->CurrentEpoch = obst;

				// Generate output file path
				Create_SUPREME_PPP_Solutions(popt, ptr->m_SatCount, popt->trace);

				if (LsqPPP._OutputFile)
					printf("%s %4d/%02d/%02d %02d:%02d:%02d ", popt->StaName,
					obst.m_Year, obst.m_Month, obst.m_Day, obst.m_Hour, obst.m_Min, (int)obst.m_Sec);

				int ctrl = 0;
				int templsq = LsqPPP.Lsq_Filter(*popt, obst, gendata, *ptr, &predata, &clkdata, prepppdata, presatinfo, zhy_flag, compareprn, iongim123);
				switch (templsq)
				{
					case -1:{
							   printf("Warrning:Observation break up...\n");
								break;
					}
					case 0:{
							   printf("Least Square Filting Failed..\n");
							   break;
					}
					case 1:{
							if (LsqPPP.flag_epoch_count > 1)
							{
								if (percentvalue != 0.0) {
									if (LsqPPP.flag_epoch_count == 2)
									{
										Rinexouthead(*popt, rinexo); Rinexouthead(*popt, rinexopre); Rinexouthead(*popt, rineosimul);
									}
									Rinexoutbodypre(rinexopre, gendata.m_Dcb, *popt, prepppdata);
									Rinexoutbody(rinexo, gendata.m_Dcb, *popt, prepppdata);
									Rinexoutbodypsim(rineosimul, gendata.m_Dcb, *popt, prepppdata);
								}
							}
							   break;
					}
					case 2:{
							   printf("Single Frequency PPP Initialize...\n");
							   break;
					}
					default:{
								break;
					}
				}
				
				if (popt->mode_filter == MODE_FILTER_FORWARD)
					ptr = ptr->m_Next;
					
				else if (popt->mode_filter == MODE_FILTER_BACKWARD)
					ptr = ptr->m_Previous;
				else
				{
					
					if (kk==0)	ptr = ptr->m_Previous;
					if (ptr->m_Previous == NULL) { 
						kk = 1; }
					if (kk == 1)
					{
						ptr = obsdata.m_FirstEpoch; kk += 1,zhy_flag=1;
					}
					if (kk == 2) 
					{
						ptr = ptr->m_Next;
					}
						
				}

				if (prepppdata.size()>1)
				{
					vector<zhydata>::iterator k = prepppdata.begin();
					prepppdata.erase(k);
				}
				if (presatinfo.size() > 1)
				{
					vector<zhyinfo>::iterator k = presatinfo.begin();
					presatinfo.erase(k);
				}

				if (ctrl) break;
			}

			printf("******* PROCESS END *******\n");

			Plot_SUPREME_Post_Results(popt); // Plot results

			process_time = process_time + 86400.0;   // Day++
		}// end loop days   
	}// end loop stations   

	return 1;
}

/// SUPREME SFPPP Processing: 1. Post mode; 2. Realtime mode
int SUPREME_PPP_Processing(char *cfg_path)
{
	printf(">>>>>>Welcome to use SUPREME (Version%.1f).\n\n", VERSION);

	ppp_config cfg;
	if (cfg_path)	cfg.ReadConfig(cfg_path);          // Read configuration file
	else cfg.ReadConfig("E:\\2020033\\PPP_Config.dat");        // Default configuratin path

		
	int fastS_D_T = 0;					
	int cycset = 0;

	if (cfg.popt->process_time == MODE_POST)           // Post poscessing
	{
		int jjt = 0, iit = 0;
		if (cfg.popt->rateorppp == 1) { jjt = 1; iit = 1; }
		else if (cfg.popt->rateorppp == 2) { jjt = 1; iit = 3; }
		else { jjt = 1; iit = 3; }	//5,3
		for (int j = 0; j < jjt; j++)
		{
			cfg.popt->schemeset = 0;//j			
			cfg.popt->cccount = 0;
			for (int i = 0; i < iit; i++)
			{
				cfg.popt->frennn = 1 + i;//1+i				
				cfg.popt->Zcount = 0;
				if (cfg.popt->rateorppp == 1) { cfg.popt->frennn == 3; fastS_D_T = 1; }
				switch (cfg.popt->frennn)
				{
				case 1:cycset = fastS_D_T > 0 ? 1 : 5; break;	//5
				case 2:cycset = fastS_D_T > 0 ? 1 : 6; break;	//6
				case 3:cycset = fastS_D_T > 0 ? 1 : 10; break;	//10
				default:break;
				}

				do {
					cfg.popt->Zcount++;
					cfg.popt->cccount++;

					SUPREME_PPP_LSF_Post(cfg);

				} while ((cfg.popt->Zcount) < cycset);//(cfg.popt->Zcount) < cycset
			}
		}
	}

	return 1;
}

