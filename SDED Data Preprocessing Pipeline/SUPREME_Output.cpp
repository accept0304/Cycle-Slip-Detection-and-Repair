#include "SUPREME_Output.h"
#include "SUPREME_Coordinate.h"
#include "SUPREME_Ionosphere.h"
#include "SUPREME_Troposphere.h"
#include "SUPREME_GraphPlot.h"

/* Create PPP Solutions Output File Path -----------------------------------------
* Parameter:
*            ppp_opt        *popt           IO           ppp option
*            unsigned int    satn           I            satellite number
*            int             Stano          I            station id
*            unsigned int    trace          I            trace mode
* Note: create output file path and file name
* ------------------------------------------------------------------------------- */
int Create_SUPREME_PPP_Solutions(ppp_option_t *popt, unsigned int satn, unsigned int trace)
{
	char freq_flag = 0;                          // frequency flag
	char mode_flag[10] = { 0 };                  // mode flag: S->static  K->kinematic
	char filter_mode[10] = { 0 };
	char sta_name[6] = { 0 };                    // station name
	char time_str[10] = { 0 };                   // time string
	char sys_str[10] = { 0 };                    // system string: G, R, E, C, GREC...
	char SUPREME_Solutions[PATH_LENGTH] = { 0 }; // project name
	char OutName[OUTPATH_LENGTH] = { 0 };        // outpath name

	unsigned int FreqNum = popt->freqn;          // frequency number

	gnsstime obst = popt->CurrentEpoch;          // current time

	memcpy(sta_name, popt->StaName, sizeof(char)* 4);

	sprintf(SUPREME_Solutions, "%s%s/", popt->IOfile.out_f.out_fp, popt->ProjName);
	memset(popt->IOfile.out_f.out_IFB, 0, sizeof(char) * 20);
	
	sprintf(time_str, "%04d%02d%02d", obst.m_Year, obst.m_Month, obst.m_Day);
	sprintf(popt->IOfile.out_f.out_IFB, "%s%s%s%s", SUPREME_Solutions, popt->StaName, time_str, "IFB.txt");
	memset(popt->IOfile.out_f.time_str, 0, sizeof(char)* 20);
	strcpy(popt->IOfile.out_f.time_str, time_str);

	switch (popt->freqn)
	{
		case 1:{freq_flag = 'S'; break; }
		case 2:{freq_flag = 'D'; break; }
		case 3:{freq_flag = 'T'; break; }
		default:{printf("Error:Frequency nubmer setting error...\n"); return 0; }
	}

	switch (popt->mode_ppp)
	{
		case MODE_STATIC:{sprintf(mode_flag,"%s","Static"); break; }
		case MODE_KINEMATIC:{sprintf(mode_flag, "%s", "Kinematic"); break; }
		default:{printf("Error:Mode setting error...\n"); return 0; }
	}

	if (popt->process_time == MODE_POST){
		switch (popt->mode_filter)
		{
			case MODE_FILTER_FORWARD:{sprintf(filter_mode, "%s", "Forward"); break; }
			case MODE_FILTER_BACKWARD:{sprintf(filter_mode, "%s", "Backward"); break; }
			case MODE_FILTER_COMBINE: {sprintf(filter_mode, "%s", "Combine"); break; }
			default:return 0;
		}
	}

	if (popt->System&PRO_SYS_GPS)
		strcat(sys_str, "G");
	if (popt->System&PRO_SYS_GLO)
		strcat(sys_str, "R");
	if (popt->System&PRO_SYS_GAL)
		strcat(sys_str, "E");
	if (popt->System&PRO_SYS_BDS)
		strcat(sys_str, "C");

	sprintf(OutName, "%s_%s_%cFPPP_%s_%s_%s", sta_name, time_str, freq_flag, sys_str, mode_flag, filter_mode);
	memset(popt->IOfile.out_f.OutName, 0, sizeof(char)*OUTPATH_LENGTH);
	strcpy(popt->IOfile.out_f.OutName, OutName);

	// PPP CRD Output File Path -----------------------------------------------
	char crdout_path[PATH_LENGTH] = { 0 };
	sprintf(crdout_path, "%s%s/", SUPREME_Solutions, time_str);
	if (!FileExist(crdout_path)) CreatePath(crdout_path);
	memset(popt->IOfile.out_f.out_fp_fpath, 0, sizeof(char)*PATH_LENGTH);
	sprintf(popt->IOfile.out_f.out_fp_fpath, "%s%s.OUT", crdout_path, OutName);

	if (trace)
	{
		// Ionosphere slant delay -------------------------------------------------
		if (popt->math_model == MODE_UNCOMBINE)
		{
			char ion_out[200] = { 0 };
			sprintf(ion_out, "%s%s/%s_STEC/", popt->IOfile.out_f.out_fp, popt->ProjName, popt->StaName);
			strcpy(popt->IOfile.out_f.out_ion, ion_out);

			if (!FileExist(ion_out)) CreatePath(ion_out);
		}

		// Trace file path --------------------------------------------------------
		char trace_path[200] = { 0 };
		sprintf(trace_path, "%s%s", SUPREME_Solutions, popt->IOfile.out_f.out_trace);
		if (!FileExist(trace_path))	CreatePath(trace_path);
		memset(popt->IOfile.out_f.out_trace_fpath, 0, sizeof(char)* PATH_LENGTH);
		sprintf(popt->IOfile.out_f.out_trace_fpath, "%s%s.trace", trace_path, OutName);

		/*FILE *tracefp = fopen(popt->IOfile.out_f.out_trace_fpath, "a+");
		if (tracefp)
			fprintf(tracefp, ">%04d-%02d-%02d %02d:%02d:%02d  Satn:%2d\n",
			popt->CurrentEpoch.m_Year, popt->CurrentEpoch.m_Month, popt->CurrentEpoch.m_Day,
			popt->CurrentEpoch.m_Hour, popt->CurrentEpoch.m_Min, (int)popt->CurrentEpoch.m_Sec, satn);
		if (tracefp) fclose(tracefp);*/

		// Log file path -----------------------------------------------------------
		char log_path[200] = { 0 };
		sprintf(log_path, "%s%s", SUPREME_Solutions, popt->IOfile.out_f.out_log);
		if (!FileExist(log_path)) CreatePath(log_path);
		memset(popt->IOfile.out_f.out_log_fpath, 0, sizeof(char)* PATH_LENGTH);
		sprintf(popt->IOfile.out_f.out_log_fpath, "%s%s.log", log_path, OutName);

		FILE *logfp = fopen(popt->IOfile.out_f.out_log_fpath, "a+");
		if (logfp)
			fprintf(logfp, ">%04d-%02d-%02d %02d:%02d:%02d  \n",
			popt->CurrentEpoch.m_Year, popt->CurrentEpoch.m_Month, popt->CurrentEpoch.m_Day,
			popt->CurrentEpoch.m_Hour, popt->CurrentEpoch.m_Min, (int)popt->CurrentEpoch.m_Sec);
		if (logfp) fclose(logfp);

		// Residual file path ------------------------------------------------------
		char residual_path[200] = { 0 };
		sprintf(residual_path, "%s%s", SUPREME_Solutions, popt->IOfile.out_f.out_res);
		if (!FileExist(residual_path))	CreatePath(residual_path);
		memset(popt->IOfile.out_f.out_res_fpath, 0, sizeof(char)* PATH_LENGTH);
		sprintf(popt->IOfile.out_f.out_res_fpath, "%s%s.resl", residual_path, OutName);

		// Variables file path -----------------------------------------------------
		char variable_paht[200] = { 0 };
		sprintf(variable_paht, "%s%s", SUPREME_Solutions, popt->IOfile.out_f.out_vbl);
		if (!FileExist(variable_paht))	CreatePath(variable_paht);
		memset(popt->IOfile.out_f.out_vbl_fpath, 0, sizeof(char)* PATH_LENGTH);
		sprintf(popt->IOfile.out_f.out_vbl_fpath, "%s%s.vbl", variable_paht, OutName);
	}

	return 1;
}

/// Write SPP CRD Output File Header
int LeastSquarePPP::LsqSPP_Output(ppp_option_t& popt)
{
	gnsstime t = popt.CurrentEpoch;

	Cart_Crd crdpos(popt.StaPos_CRD[0], popt.StaPos_CRD[1], popt.StaPos_CRD[2]);
	Coordinate neu;
	neu.XYZ.setXYZ(popt.StaPos_SPP[0], popt.StaPos_SPP[1], popt.StaPos_SPP[2]);

	neu._xyz2blh();

	/*if (!FileExist(popt.IOfile.out_f.out_fp_fpath))
	Write_PPP_CRD_Header(popt, popt.IOfile.out_f.out_fp_fpath);*/

	FILE *fp = fopen(popt.IOfile.out_f.out_fp_fpath, "a+");

	if (Norm(popt.StaPos_CRD, 3) == 0)
	{
		printf("(%s)-%2d %15.4f  %15.4f  %15.4f  %10.4f %10.4f %10.4f\n",
			popt.StaName, m_SatValidN, popt.StaPos_SPP[0], popt.StaPos_SPP[1], popt.StaPos_SPP[2],
			neu.BLH._B*180.0 / PI, neu.BLH._L*180.0 / PI, neu.BLH._H);

		if (fp)
		{
			fprintf(fp, "%04d:%02d:%02d:%02d:%02d:%02d  %2d  999.9  999.9  999.9 %17.4f%17.4f%17.4f %8.4f",
				t.m_Year, t.m_Month, t.m_Day, t.m_Hour, t.m_Min, (int)t.m_Sec, m_SatValidN,
				popt.StaPos_SPP[0], popt.StaPos_SPP[1], popt.StaPos_SPP[2], m_TropZenith);

			//% Satellite number used in ppp
			fprintf(fp, " %5d %5d %5d %5d\n", _GpsSatNum, _GloSatNum, _BdsSatNum, _GalSatNum);
		}
	}
	else
	{
		//% Output on screen ----------------------------------------------------------------
		printf("(%s):%2d ", popt.StaName, popt.SatN_SPP);

		//if (est_xyz)
		{
			neu._xyz2neu(crdpos);
			printf("%8.4f  %8.4f  %8.4f  ", neu.NEU._N, neu.NEU._E, neu.NEU._U);
		}
		printf("\n");
		//% ---------------------------------------------------------------------------------

		//% Output in file ------------------------------------------------------------------
		if (fp)
		{
			fprintf(fp, "%04d:%02d:%02d:%02d:%02d:%02d %4d",
				t.m_Year, t.m_Month, t.m_Day, t.m_Hour, t.m_Min, (int)t.m_Sec, popt.SatN_SPP);

			//if (est_xyz)
			{
				fprintf(fp, "%8.4f %8.4f %8.4f  %17.4f %17.4f %17.4f",
					neu.NEU._N, neu.NEU._E, neu.NEU._U,
					popt.StaPos_SPP[0], popt.StaPos_SPP[1], popt.StaPos_SPP[2]);
			}
			//if (est_rclk)
			{
				fprintf(fp, "%15.3f %5.1f", popt.Rclk / LIGHTSPEED*1e9,popt.GDOP);
			}
			
			//% Satellite number used in ppp
			fprintf(fp, " %5d %5d %5d %5d", _GpsSatNum, _GloSatNum, _BdsSatNum, _GalSatNum);

			fprintf(fp, "\n");
		}
	}
	if (fp) fclose(fp);

	return 1;
}

int STEC_Header(ppp_option_t &popt, ObsFileHeader &header, char *ionf)
{
	FILE *fp = fopen(ionf, "w");
	if (fp)
	{
		fprintf(fp, " Station(name,number):                %s\n", popt.StaName);
		fprintf(fp, " Antenna(type,number):%s\n", header.m_AntennaType);
		fprintf(fp, "Receiver(type,number):%s\n", header.m_ReceiverType);
		fprintf(fp, " Receiver coordinates:%20.3f%20.3f%20.3f\n", popt.StaPos_CRD[0], popt.StaPos_CRD[1], popt.StaPos_CRD[2]);
		fprintf(fp, "       Date(UTC)          PRN SAT_ID   STEC[TECu] RMS[TECu]  SatDcb RecDcb Sat_Elev  Sat_Azim      Sat_X[m]      Sat_Y[m]       Sat_Z[m]   Ipp_Elev  Ipp_Lat   Ipp_Lon  Coef_FRE1    Coef_FRE2  Coef_FRE3    ObsComb  Flag\n");
		fprintf(fp, "--------------------------------------------------------------------------------------------------------------------------------------------\n");
	}
	if (fp) fclose(fp);

	return 1;
}

int LeastSquarePPP::PPP_STEC_Output(ppp_option_t &popt, PPPObsData &inputdata)
{
	FILE *fp = NULL;
	char sf_ion_outf[200] = { 0 };

	gnsstime obst = popt.CurrentEpoch;

	double stablh[3] = { 0 };
	ecef2pos(m_Sta, stablh);

	sprintf(sf_ion_outf, "%s%s/%s%03d0.%02di",
		popt.IOfile.out_f.out_fp, popt.ProjName, popt.StaName, obst.m_Doy, obst.m_Year % 100);

	if (_FreqNum == 1 && flag_epoch_count == 1)
		STEC_Header(popt, _ObsHeader, sf_ion_outf);
	else if (flag_epoch_count == 0)
		STEC_Header(popt, _ObsHeader, sf_ion_outf);

	fp = fopen(sf_ion_outf, "a+");
	if (fp)
	{
		for (int i = 0; i < m_SatValidN; ++i)
		{
			fprintf(fp, "%4d/%02d/%02d %02d:%02d:%05.2f    ",
				obst.m_Year, obst.m_Month, obst.m_Day, obst.m_Hour, obst.m_Min, obst.m_Sec);

			double L1 = 0.0, L2 = 0.0, L_gf = 0.0;
			double wavelength1 = 0.0, wavelength2 = 0.0;
			wavelength1 = Get_WaveLength(m_Prn[i], 1, 0);
			wavelength2 = Get_WaveLength(m_Prn[i], 2, 0);
			L1 = inputdata.m_L1[i] * wavelength1; L2 = inputdata.m_L2[i] * wavelength2;
			L_gf = L1 - L2;

			unsigned int prn = 0;
			int slot = 0;
			char sys = 0;
			sys = GetSysPrn(m_Prn[i], prn);

			double ion_value = X(_ParaN_const + _PsuedoObsN + i + 1, 1);
			double var_ion = 0.0;
			var_ion = Nbb(_ParaN_const + _PsuedoObsN + i + 1, _ParaN_const + _PsuedoObsN + i + 1);
			var_ion = sqrt(var_ion);

			if (flag_epoch_count == 1 && _FreqNum == 1 && _IonModel == ION_MODEL_EST_WHITE)
				ion_value = X(_ParaN_const + _PsuedoObsN + m_SatValidN + i + 1, 1);

			//% convert meter to tecu
			if (sys == 'G')
			{
				ion_value = Meter2TECu(ion_value, GPS_FREQU_L1);
				var_ion = Meter2TECu(var_ion, GPS_FREQU_L1);
			}
			else if (sys == 'C')
			{
				ion_value = Meter2TECu(ion_value, BDS_FREQU_B1);
				var_ion = Meter2TECu(var_ion, BDS_FREQU_B1);
			}
			else if (sys == 'R')
			{
				slot = m_Slot[i];

				ion_value = Meter2TECu(ion_value, GLO_FREQU_L1(slot));
				var_ion = Meter2TECu(var_ion, GLO_FREQU_L1(slot));
			}
			else if (sys == 'E')
			{
				ion_value = Meter2TECu(ion_value, GAL_FREQU_E1);
				var_ion = Meter2TECu(var_ion, GAL_FREQU_E1);
			}

			double ipp_pos[3] = { 0 }, azel[2] = { 0.0 }, xs[3] = { 0 };
			azel[0] = m_SatInfo.m_Azi[i];
			azel[1] = m_SatInfo.m_Ele[i];
			xs[0] = m_SatInfo.m_Xs[3 * i + 0];
			xs[1] = m_SatInfo.m_Xs[3 * i + 1];
			xs[2] = m_SatInfo.m_Xs[3 * i + 2];

			IonPPP(stablh, azel, WGS84_A / 1000, 450, ipp_pos);

			fprintf(fp, "%c%02d    %02d %12.4f ", sys, prn, prn, ion_value);

			fprintf(fp, "%9.4f %6d %6d ", var_ion, 0, 0);

			fprintf(fp, "%9.2f %9.2f %14.3f %14.3f %14.3f ", R2D*azel[1], R2D*azel[0], xs[0], xs[1], xs[2]);

			fprintf(fp, "%9.2f %9.2f %9.2f ", 0.0, ipp_pos[0] * R2D, ipp_pos[1] * R2D);

			if (sys == 'G')
			{
				double gps_cof = 0;
				gps_cof = (GPS_FREQU_L1*GPS_FREQU_L1)*(GPS_FREQU_L2*GPS_FREQU_L2) / ((GPS_FREQU_L1*GPS_FREQU_L1 - GPS_FREQU_L2*GPS_FREQU_L2)*ION_FACTEC);

				fprintf(fp, "%11.7f %11.7f %11.7f    ", -9.5243700, 9.5243700, 0.0);
				fprintf(fp, "C1W_C2W    %2d %15.3f\n", popt.QC_Prn[prn - 1], L_gf);
			}
			else if (sys == 'C')
			{
				double bds_cof = 0;

				bds_cof = (BDS_FREQU_B1*BDS_FREQU_B1)*(BDS_FREQU_B2*BDS_FREQU_B2) / ((BDS_FREQU_B1*BDS_FREQU_B1 - BDS_FREQU_B2*BDS_FREQU_B2)*ION_FACTEC);

				fprintf(fp, "%11.7f %11.7f %11.7f    ", -bds_cof, bds_cof, 0);

				fprintf(fp, "C2I_C7I     1\n");
			}
			else if (sys == 'R')
			{
				double glo_cof = 0;
				slot = m_Slot[i];

				glo_cof = (GLO_FREQU_L1(slot)*GLO_FREQU_L1(slot))*(GLO_FREQU_L2(slot)*GLO_FREQU_L2(slot)) / ((GLO_FREQU_L1(slot)*GLO_FREQU_L1(slot) - GLO_FREQU_L2(slot)*GLO_FREQU_L2(slot))*ION_FACTEC);

				fprintf(fp, "%11.7f %11.7f %11.7f    ", -glo_cof, glo_cof, 0);

				fprintf(fp, "C1P_C2P     1\n");
			}
			else if (sys == 'E')
			{
				double gal_cof = 0;

				gal_cof = (GAL_FREQU_E1*GAL_FREQU_E1)*(GAL_FREQU_E5a*GAL_FREQU_E5a) / ((GAL_FREQU_E1*GAL_FREQU_E1 - GAL_FREQU_E5a*GAL_FREQU_E5a)*ION_FACTEC);

				fprintf(fp, "%11.7f %11.7f %11.7f    ", -gal_cof, gal_cof, 0);

				fprintf(fp, "C1X_C5X     1\n");
			}
		}
	}
	if (fp) fclose(fp);

	return 1;
}

/// Write PPP CRD Output File Header
int Write_PPP_CRD_Header(ppp_option_t& popt, char* outpath)
{
	FILE *fp = NULL;
	fp = fopen(outpath, "a+");
	if (fp)
	{
		fprintf(fp, "# Station Name:%s\n", popt.StaName);
		fprintf(fp, "# Antenna Type:%s\n", popt.antenna_type);
		fprintf(fp, "# ReceiverType:%s\n", popt.receiver_type);

		if (popt.mode_ppp == MODE_STATIC)
			fprintf(fp, "# Process Mode:Static\n");
		else if (popt.mode_ppp == MODE_KINEMATIC)
			fprintf(fp, "# Process Mode:Kinematic\n");

		if (popt.process_time == MODE_REALTIME)
		{
			fprintf(fp, "# SSR Source  :%s\n", popt.ssr_mountpoint);
			fprintf(fp, "# NAV Source  :%s\n", popt.nav_mountpoint);
		}
		fprintf(fp, "# Coarse Coord:%15.4f %15.4f %15.4f\n", popt.StaPos_SPP[0], popt.StaPos_SPP[1], popt.StaPos_SPP[2]);
		fprintf(fp, "# Created   by:C.Zhao\n");
		fprintf(fp, "#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
		fprintf(fp, "#    Time(GPST)     SATN    N(m)       E(m)       U(m)              X(m)      Cov_X(m)          Y(m)      Cov_Y(m)          Z(m)      Cov_Z(m)         RClk(s)       Cov_Rclk(m)     ZTD(m)  Cov_ZTD(m)   G     R     C     E    ISB   \n");
		fprintf(fp, "#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	}

	if (fp) fclose(fp);

	return 1;
}

/// Write PPP CRD Output File Header
int LeastSquarePPP::PPP_CRD_Trop_Output(ppp_option_t& popt)
{
	gnsstime t = popt.CurrentEpoch;

	Cart_Crd crdpos(popt.StaPos_CRD[0], popt.StaPos_CRD[1], popt.StaPos_CRD[2]);//有问题

	Coordinate neu;
	neu.XYZ.setXYZ(m_Sta[0], m_Sta[1], m_Sta[2]);

	neu._xyz2blh();

	if (!FileExist(popt.IOfile.out_f.out_fp_fpath))
		Write_PPP_CRD_Header(popt, popt.IOfile.out_f.out_fp_fpath);
	
	FILE *fp = fopen(popt.IOfile.out_f.out_fp_fpath, "a+");
	FILE *fp2 =NULL;//fopen("out/TEST2/IGGUBX1.pos", "a+"); //
	FILE *fp3 =  fopen(popt.IOfile.out_f.out_IFB, "a+");
	//% Output on screen ----------------------------------------------------------------
	printf("(%s):%2d ", popt.StaName, m_SatValidN);

	int num = popt.IFB_num;
	if (fp3) 
	{
		for (int i = 0; i < num; i++)
		{

			double ifb = X_Trans(_ParaN_const + m_SatValidN * _IonNum + popt.CUR_mum + i + 1, 1);
			fprintf(fp3, "%8.4f", ifb);

		}
		fprintf(fp3, "\n");
	}
	


	if (Norm(popt.StaPos_CRD, 3) == 0)
	{
		printf(" %15.4f  %15.4f  %15.4f  %10.4f %10.4f %10.4f",
			m_Sta[0], m_Sta[1], m_Sta[2], neu.BLH._B*180.0 / PI, neu.BLH._L*180.0 / PI, neu.BLH._H);

		neu.NEU._N = 999.999; neu.NEU._E = 999.999; neu.NEU._U = 999.999;
	}
	else
	{
		neu._xyz2neu(crdpos);
		printf("%8.4f  %8.4f  %8.4f  ", neu.NEU._N, neu.NEU._E, neu.NEU._U);
		printf("%10.4f", m_TropZenith);
		printf("%20.4f", (X_Trans(3*(est_xyz)+1, 1) ));
		printf("%20.4f", (X_Trans(3 * (est_xyz)+2, 1)));
		
	}
	printf("\n");
	//% ---------------------------------------------------------------------------------
	
	//% Output in file ------------------------------------------------------------------
	if (fp)
	{
		fprintf(fp, "%04d:%02d:%02d:%02d:%02d:%02d %3d",
			t.m_Year, t.m_Month, t.m_Day, t.m_Hour, t.m_Min, (int)t.m_Sec, m_SatValidN);
	}

	if (fp2)
	{
		fprintf(fp2, "%04d:%02d:%02d:%02d:%02d:%02d %3d",
			t.m_Year, t.m_Month, t.m_Day, t.m_Hour, t.m_Min, (int)t.m_Sec, m_SatValidN);
	}

	double var_x = 0.0, var_y = 0.0, var_z = 0.0, var_rclk = 0.0;
	double var_trop = 0.0, var_Ntrop = 0.0, var_Etrop = 0.0;
	/*if (popt.estpara.xyz)
	{
		var_x = sqrt(Nbb(1, 1));
		var_y = sqrt(Nbb(2, 2));
		var_z = sqrt(Nbb(3, 3));
	}*/
	if (popt.estpara.xyz)
	{
		var_x = sqrt(P_Trans(1, 1));
		var_y = sqrt(P_Trans(2, 2));
		var_z = sqrt(P_Trans(3, 3));
	}
	//var_rclk = sqrt(Nbb(3 * est_xyz + _CbiasNum + est_rclk, 3 * est_xyz + _CbiasNum + est_rclk));

	var_rclk = sqrt(P_Trans(3 * est_xyz + _CbiasNum + est_rclk, 3 * est_xyz + _CbiasNum + est_rclk));

	//if (est_xyz)
	{
		fprintf(fp, "%11.4f%11.4f%11.4f%18.4f%10.4f%18.4f%10.4f%18.4f%10.4f",
			neu.NEU._N, neu.NEU._E, neu.NEU._U,
			m_Sta[0], var_x, m_Sta[1], var_y, m_Sta[2], var_z);

		if (fp2)
		fprintf(fp2, "%15.10f %15.10f %12.4f\n", neu.BLH._B*180.0 / PI, neu.BLH._L*180.0 / PI, neu.BLH._H);
	}

	// Output Receiver Clock
	fprintf(fp, "%10.4f %10.4f", X_Trans(3 * est_xyz + _CbiasNum + est_rclk, 1) , var_rclk);//%25.12e

	// Output Trop Delay and Gradient ------------------------------------------
	if (est_trop)
		var_trop = sqrt(P_Trans(3 * est_xyz + _CbiasNum + est_rclk*_SystemN + est_trop,
		3 * est_xyz + _CbiasNum + est_rclk*_SystemN + est_trop));
	fprintf(fp, "%11.4f%10.4f", m_TropZenith, var_trop);

	if (_TropModel == TROP_GRADIENT)
	{
		var_Ntrop = sqrt(P_Trans(3 * est_xyz + _CbiasNum + est_rclk*_SystemN + est_trop + 1,
			3 * est_xyz + _CbiasNum + est_rclk*_SystemN + est_trop + 1));
		var_Etrop = sqrt(P_Trans(3 * est_xyz + _CbiasNum + est_rclk*_SystemN + est_trop + 2,
			3 * est_xyz + _CbiasNum + est_rclk*_SystemN + est_trop + 2));
		fprintf(fp, "%11.3f%10.3f", m_TropGradientN * 1000, var_Ntrop * 1000);
		fprintf(fp, "%11.3f%10.3f", m_TropGradientE * 1000, var_Etrop * 1000);
	}
	// Output Trop Delay and Gradient end ---------------------------------------

	//% Satellite number of each system used in ppp
	fprintf(fp, " %5d %5d %5d %5d", _GpsSatNum, _GloSatNum, _BdsSatNum, _GalSatNum);

	//输出北斗ISB;;;;;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	/*double isb = X_Trans(3 * est_xyz + 2, 1);
	double isb_ns = isb / LIGHTSPEED * 1e9;
	fprintf(fp, " %10.4f ", isb);
	fprintf(fp, "%10.4f", isb_ns);*/
	
	double B3 = X_Trans(3 * (est_xyz)+2, 1);
	
	fprintf(fp, "%20.4f", (B3 ));//b3
	

	fprintf(fp, "\n");

	if (fp) fclose(fp);
	if (fp2) fclose(fp2);
	if (fp3)fclose(fp3);
	return 1;
}

int LeastSquarePPP::PPP_Res_Output(ppp_option_t &popt)
{
	FILE *resfp = fopen(popt.IOfile.out_f.out_res_fpath, "a+");
	if (resfp)
		fprintf(resfp, "#%04d %02d %02d %02d %02d %02d  satn: %2d (",
		popt.CurrentEpoch.m_Year, popt.CurrentEpoch.m_Month, popt.CurrentEpoch.m_Day,
		popt.CurrentEpoch.m_Hour, popt.CurrentEpoch.m_Min, (int)popt.CurrentEpoch.m_Sec,
		m_SatValidN);
	for (int i = 0; i < m_SatValidN; ++i)
	{
		char sys = 0;
		int prn = 0;
		switch (sys = GetSystem(m_Prn[i]))
		{
			case 'G':{prn = m_Prn[i] - MIN_GPS_PRN + 1; break; }
			case 'R':{prn = m_Prn[i] - MIN_GLO_PRN + 1; break; }
			case 'E':{prn = m_Prn[i] - MIN_GAL_PRN + 1; break; }
			case 'C':{prn = m_Prn[i] - MIN_BDS_PRN + 1; break; }
			default:break;
		}
		if (resfp) fprintf(resfp, "%c%02d", sys, prn);
	}
	if (resfp) fprintf(resfp, ")\n");
	for (int i = 0; i < m_SatValidN; ++i)
	{
		char sys = 0;
		unsigned int prn = 0;
		sys = GetSysPrn(m_Prn[i], prn);

		if (_MathModel == MODE_IF_COMBINE)
		{
			if (_FreqNum == 2)
			{
				if (resfp)
					fprintf(resfp, "%c%02d %8.4f  %8.4f\n",
					sys, prn, V_zhy(_PsuedoObsNum + _PsuedoObsN + 2 * i + 1, 1), V_zhy(_PsuedoObsNum + _PsuedoObsN + 2 * i + 2, 1));
			}
			else if (_FreqNum == 1)
			{
				if (resfp)
					fprintf(resfp, "%c%02d %8.4f\n",
					sys, prn, V_zhy(_PsuedoObsNum + _PsuedoObsN + i + 1, 1));
			}
		}
		else if (_MathModel == MODE_UNCOMBINE||_MathModel==MODE_UNCOMBINE_ION_CONT)
		{
			if (_FreqNum == 2)
			{
				if (resfp)
					fprintf(resfp, "%c%02d %8.4f  %8.4f  %8.4f  %8.4f\n", sys, prn,
					V_zhy(_PsuedoObsNum + _PsuedoObsN + 4 * i + 1 + 0, 1), V_zhy(_PsuedoObsNum + _PsuedoObsN + 4 * i + 1 + 1, 1),
					V_zhy(_PsuedoObsNum + _PsuedoObsN + 4 * i + 1 + 2, 1), V_zhy(_PsuedoObsNum + _PsuedoObsN + 4 * i + 1 + 3, 1));
			}
			else if (_FreqNum == 1 && flag_epoch_count>1)
			{
				if (resfp)
					fprintf(resfp, "%c%02d %8.4f  %8.4f\n", sys, prn,
					V_zhy(_PsuedoObsNum + _PsuedoObsN + 2 * i + 1, 1), V_zhy(_PsuedoObsNum + _PsuedoObsN + 2 * i + 2, 1));
			}
		}
	}
	if (resfp) fclose(resfp);

	// Output Amb Residuals ===================
	/*if (_ResAmb.size())
	{
		FILE *ambresfp = fopen("out/amb_res.resl", "a+");
		if (ambresfp)
			fprintf(ambresfp, "#%04d %02d %02d %02d %02d %02d  satn: %2d (",
			popt.CurrentEpoch.m_Year, popt.CurrentEpoch.m_Month, popt.CurrentEpoch.m_Day,
			popt.CurrentEpoch.m_Hour, popt.CurrentEpoch.m_Min, (int)popt.CurrentEpoch.m_Sec,
			_ResAmb.size());
		for (int i = 0; i < _ResAmb.size(); ++i)
		{
			char sys = 0;
			int prn = 0;
			switch (sys = GetSystem(_Prn_common[i]))
			{
				case 'G':{prn = _Prn_common[i] - MIN_GPS_PRN + 1; break; }
				case 'R':{prn = _Prn_common[i] - MIN_GLO_PRN + 1; break; }
				case 'E':{prn = _Prn_common[i] - MIN_GAL_PRN + 1; break; }
				case 'C':{prn = _Prn_common[i] - MIN_BDS_PRN + 1; break; }
				default:break;
			}
			if (ambresfp) fprintf(ambresfp, "%c%02d", sys, prn);
		}
		if (ambresfp) fprintf(ambresfp, ")\n");
		
		for (int i = 0; i < _ResAmb.size(); ++i)
		{
			char sys = 0;
			unsigned int prn = 0;
			sys = GetSysPrn(_Prn_common[i], prn);

			if (resfp)
				fprintf(resfp, "%c%02d %8.4f\n", sys, prn, _ResAmb[i]);
		}

		if (ambresfp) fclose(ambresfp);
	}*/
	// ========================================

	return 1;
}

int LeastSquarePPP::PPP_STEC_Sat_Output(ppp_option_t &popt)
{
	unsigned int i = 0;
	int AmbIndexOffset = 0;
	for (i = 0; i < m_SatValidN; ++i)
	{
		char sys = 0;
		unsigned int prn;
		char ion_out[200] = { 0 };
		int slot = 0;

		if (_RefSatValidFlag && m_Prn[i] == _RefSatPrn){
			AmbIndexOffset = 1;
			continue;
		}

		sys = GetSysPrn(m_Prn[i], prn);

		double wavelength1 = 0.0, wavelength2 = 0.0;
		wavelength1 = Get_WaveLength(m_Prn[i], 1, 0);             // Wavelength of f1
		wavelength2 = Get_WaveLength(m_Prn[i], 2, 0);             // Wavelength of f2

		sprintf(ion_out, "%s%c%02d_%s.ion",
			popt.IOfile.out_f.out_ion, sys, prn, popt.IOfile.out_f.time_str);

		double ion_delay = X_Trans(3 * est_xyz + _CbiasNum + est_rclk*_SystemN + _TropNum + _PsuedoObsN + i + 1, 1);
		double amb1 = 0, amb2 = 0;
		amb1 = X_Trans(3 * est_xyz + _CbiasNum + est_rclk*_SystemN + _TropNum + m_SatValidN + _AmbNum * i + 1 - AmbIndexOffset, 1);
		if (_FreqNum==2)
			amb2 = X_Trans(3 * est_xyz + _CbiasNum + est_rclk*_SystemN + _TropNum + m_SatValidN + _AmbNum * i + 2, 1);

		// Convert meter to TECU
		/*if (sys == 'G')
		{
			ion_delay = Meter2TECu(ion_delay, GPS_FREQU_L1);
		}
		else if (sys == 'C')
		{
			if (strstr(CodeType1[BDS], "C2I"))
				ion_delay = Meter2TECu(ion_delay, BDS_FREQU_B1);
			else if (strstr(CodeType1[BDS], "C7I"))
				ion_delay = Meter2TECu(ion_delay, BDS_FREQU_B2);
		}
		else if (sys == 'R')
		{
			slot = m_Slot[i];

			ion_delay = Meter2TECu(ion_delay, GLO_FREQU_L1(slot));
		}
		else if (sys == 'E')
		{
			ion_delay = Meter2TECu(ion_delay, GAL_FREQU_E1);
		}*/
		// Convert meter to TECU end

		FILE *fp = fopen(ion_out, "a+");
		if (fp)
		{
			/*fprintf(fp, "%4d:%02d:%02d:%02d:%02d:%02d %12.4f\n",
				popt.CurrentEpoch.m_Year, popt.CurrentEpoch.m_Month, popt.CurrentEpoch.m_Day,
				popt.CurrentEpoch.m_Hour, popt.CurrentEpoch.m_Min, (int)popt.CurrentEpoch.m_Sec,
				ion_delay);*/
			fprintf(fp, "%4d:%02d:%02d:%02d:%02d:%02d %12.4f ",
				popt.CurrentEpoch.m_Year, popt.CurrentEpoch.m_Month, popt.CurrentEpoch.m_Day,
				popt.CurrentEpoch.m_Hour, popt.CurrentEpoch.m_Min, (int)popt.CurrentEpoch.m_Sec,
				amb1/wavelength1);
			if (_FreqNum > 1)
				fprintf(fp, "%12.4f ", amb2 / wavelength2);
			fprintf(fp, "\n");
		}
		if (fp) fclose(fp);
	}

	return 1;
}

/// Least square filting ppp result output
int LeastSquarePPP::LsqPPP_Output(ppp_option_t &popt, PPPObsData &pppdata)
{
	gnsstime t = popt.CurrentEpoch;

	// Output CRD and Trop results
	PPP_CRD_Trop_Output(popt);

	// Write obs residuals to file
	if (_traceMode)
	{
		PPP_Res_Output(popt);
	}

	// Output Troposphere Series
	if (0)
	{
		Coordinate neu;
		neu.XYZ.setXYZ(popt.StaPos_CRD[0], popt.StaPos_CRD[1], popt.StaPos_CRD[2]);
		neu._xyz2blh();//这个坐标是怎么给出的

		char trop_path[PATH_LENGTH] = { 0 };
		sprintf(trop_path, "%s%s/%s%s%s", popt.IOfile.out_f.out_fp, popt.ProjName, popt.StaName, popt.IOfile.out_f.time_str, ".trop");
		FILE *fp = fopen(trop_path, "a+");
		if (fp)
		{
			fprintf(fp, "%04d:%02d:%02d:%02d:%02d:%02d",
				t.m_Year, t.m_Month, t.m_Day, t.m_Hour, t.m_Min, (int)t.m_Sec);
			fprintf(fp, "%10.4f %10.4f %10.4f%15.3f\n",
				neu.BLH._B*180.0 / PI, neu.BLH._L*180.0 / PI, neu.BLH._H, m_TropZenith * 1000);
		}
		if (fp) fclose(fp);
	}

	// Output slant ionosphere slant delay
	if (_MathModel == MODE_UNCOMBINE||_MathModel==MODE_UNCOMBINE_ION_CONT)
	{
		if (_traceMode)
		{
			if (_IonModel == ION_MODEL_EST_WHITE)
			{
				if (popt.CurrentEpoch.m_Hour)
				PPP_STEC_Sat_Output(popt);                      // single satellite stec

				//PPP_STEC_Output(popt, pppdata);                  // total satellite stec
			}
		}
	}
	// Output slant ionosphere delay end

	return 1;
}

int LeastSquarePPP::UdState2KF(ppp_option_t* popt, unsigned int SatN, PPPObsData& pppdata, vector<double>ion_GIM, PreciseData& predata)
{

	int orbn = 0; int g_or_c = 0; double freq1 = 0; double freq2 = 0; double ion_ = 0; int g_5 = 0;
	double P1 = 0, P2 = 0, P3 = 0, P4 = 0, P5 = 0, L1 = 0, L2 = 0, L3 = 0, L4 = 0, L5 = 0; vector<int>prn_;
	vector<int>ion_index; vector<double>Wavelength1, Wavelength2, Wavelength3, Wavelength4, Wavelength5;
	//_GIM_P1P2 = 2;//测试
	double wavelength1 = 0,wavelength2 = 0,wavelength3 = 0,wavelength4 = 0,wavelength5 = 0;
	char systype[2]; double f11=0, f22=0;
	for (int i_sat = 0; i_sat < SatN; i_sat++)//一颗卫星 多个频率 那怎么记录呢 先判断系统  选择要用的 频率
	{
		unsigned int satno = pppdata.m_Prn[i_sat];
		if (GetSystem(satno) == 'G' ) { 
			orbn = 0;
			wavelength1 = Get_WaveLength(satno, 1, orbn);//L1
			wavelength2 = Get_WaveLength(satno, 2, orbn);//L2
			wavelength3 = Get_WaveLength(satno, 5, orbn);//L5
		}
		if (GetSystem(satno) ==  'E') 
		{
			orbn = 0; 
			wavelength1 = Get_WaveLength(satno, 1, orbn);//E1
			wavelength2 = Get_WaveLength(satno, 2, orbn);//E5A
			wavelength3 = Get_WaveLength(satno, 3, orbn);//E5b
			wavelength4 = Get_WaveLength(satno, 4, orbn);//E5(A+B)
			wavelength5 = Get_WaveLength(satno, 5, orbn);//E6
			f11=Get_Frequency(satno, 1, orbn);
			f22 = Get_Frequency(satno, 2, orbn);
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
		
		Wavelength1.push_back(wavelength1);
		Wavelength2.push_back(wavelength2);
		Wavelength3.push_back(wavelength3);
		Wavelength4.push_back(wavelength4);
		Wavelength5.push_back(wavelength5);
	}
	if (_GIM_P1P2 == P1P2_const && _MathModel != MODE_IF_COMBINE && _MathModel != MODE_UOFC_COMBINE && _MathModel != 4)//双频电离层
	{
		for (int i_sat = 0; i_sat < SatN; i_sat++)
		{
			P1 = pppdata.m_Pr1[i_sat]; P2 = pppdata.m_Pr2[i_sat];
			if (Wavelength1[i_sat] == 0 || Wavelength2[i_sat] == 0 || P1 == 0 || P2 == 0) {
				ION_.push_back(0);
				continue;
			}
			ion_ = (P1 - P2) / (1.0 - SQR(Wavelength2[i_sat] / Wavelength1[i_sat]));
			ION_.push_back(ion_);
		}
		double ion_over = 0;
		for (int i = 0; i < SatN; i++) ion_over += ION_[i];
		for (int sat_x = 0; sat_x < SatN; sat_x++)
		{
			if (ION_[sat_x] == 0.0)  ION_[sat_x] = ion_over / ION_.size();
		}
	}
	if (_MathModel == 4)
	{
		for (int i_sat = 0; i_sat < SatN; i_sat++)
		{
			P1 = pppdata.m_Pr1[i_sat]; P2 = pppdata.m_Pr2[i_sat];
			if (Wavelength1[i_sat] == 0 || Wavelength2[i_sat] == 0 || P1 == 0 || P2 == 0) {
				ION_.push_back(0);
				continue;
			}
		
			double ion_1 = (P1 - P2) / (1.0 - SQR(f11 / f22));  // 自己算的
			ion_ = (P1 - P2) / (1.0 - SQR(Wavelength2[i_sat] / Wavelength1[i_sat]));//原来没改之前
			//ion_ = ion_ * SQR(LIGHTSPEED / Wavelength1[i_sat]) / 40.3*1e-16;
			ION_.push_back(ion_);	
		}
		double ion_over = 0;
		for (int i = 0; i < SatN; i++) ion_over += ION_[i];
		for (int sat_x = 0; sat_x < SatN; sat_x++)
		{
			if (ION_[sat_x] == 0.0)  { 
				ION_[sat_x] = ion_over / ION_.size(); 
			} 
		}
	}
	if (_GIM_P1P2 == GIM_const)//格网电离层
	{
		for (int i = 0; i < ion_GIM.size(); i++)
		{
			ION_.push_back(ion_GIM[i]);
		}
	}
	if (_MathModel != MODE_IF_COMBINE && _MathModel != 4) {
		for (int i_sat = 0; i_sat < SatN; i_sat++)
		{
			L1 = pppdata.m_L1[i_sat];
			P1 = pppdata.m_Pr1[i_sat];
			L1 *= Wavelength1[i_sat];
			double bias_ = L1 - P1 + 2 * ION_[i_sat];
			BIAS_.push_back(bias_);//模糊度初值；
			prn_.push_back(pppdata.m_Prn[i_sat]);
		}
	}
	if (_MathModel == 4) {
		for (int i_sat = 0; i_sat < SatN; i_sat++)
		{
			L1 = pppdata.m_L1[i_sat];
			P1 = pppdata.m_Pr1[i_sat];
			L1 *= Wavelength1[i_sat];
			L2 = pppdata.m_L2[i_sat];
			P2 = pppdata.m_Pr2[i_sat];
			L2 *= Wavelength2[i_sat];
			L3 = pppdata.m_L3[i_sat];
			P3 = pppdata.m_Pr3[i_sat];
			L3 *= Wavelength3[i_sat];
			L4 = pppdata.m_L4[i_sat];
			P4 = pppdata.m_Pr4[i_sat];
			L4 *= Wavelength4[i_sat];
			L5 = pppdata.m_L5[i_sat];
			P5 = pppdata.m_Pr5[i_sat];
			L5 *= Wavelength5[i_sat];


			unsigned int satno1 = pppdata.m_Prn[i_sat];
			double bias_ = L1 - P1 + 2 * ION_[i_sat];
			double bias_1 = L2 - P2 + 2 * ION_[i_sat] * SQR(Wavelength2[i_sat] / Wavelength1[i_sat]);
			double bias_2 = (L3 == 0 || P3 == 0) ? 0 : (L3 - P3 + 2 * ION_[i_sat] * SQR(Wavelength3[i_sat] / Wavelength1[i_sat]));
			double bias_3 = (L4 == 0 || P4 == 0) ? 0 :( L4 - P4 + 2 * ION_[i_sat] * SQR(Wavelength4[i_sat] / Wavelength1[i_sat]));
			double bias_4 = (L5 == 0 || P5 == 0) ? 0 : (L5 - P5 + 2 * ION_[i_sat] * SQR(Wavelength5[i_sat] / Wavelength1[i_sat]));

			BIAS_.push_back(bias_);//模糊度初值；
			BIAS_1.push_back(bias_1);//模糊度初值；
			BIAS_2.push_back(bias_2);//模糊度初值；
			BIAS_3.push_back(bias_3);//模糊度初值；
			BIAS_4.push_back(bias_4);//模糊度初值；
			prn_.push_back(pppdata.m_Prn[i_sat]);
		}
	}
	if (_MathModel == MODE_IF_COMBINE)
	{
		for (int i = 0; i < SatN; i++)
		{
			if (pppdata.m_Pr1[i] * pppdata.m_Pr2[i] * pppdata.m_L2[i] * pppdata.m_L1[i] == 0)
			{
				--SatN; --m_SatValidN; --i;
				continue;
			}
			double C1 = SQR(Wavelength2[i]) / (SQR(Wavelength2[i]) - SQR(Wavelength1[i]));
			double C2 = -SQR(Wavelength1[i]) / (SQR(Wavelength2[i]) - SQR(Wavelength1[i]));
			P1 = C1 * pppdata.m_Pr1[i] + C2 * pppdata.m_Pr2[i];
			L1 = C1 * pppdata.m_L1[i] * Wavelength1[i] + C2 * pppdata.m_L2[i] * Wavelength2[i];
			double bias = L1 - P1;
			BIAS_.push_back(bias);
			prn_.push_back(pppdata.m_Prn[i]);
		}
	}
	popt->ifb_mode = 0;
	if (est_xyz)
	{
		Udpos(popt);
	}

	udtrop(popt);
	if (_MathModel == 1 || _MathModel == 3||_MathModel == 4)
	{
		udion(popt, SatN, prn_);
	}
	//更新模糊度
	udbias(popt,SatN,prn_);
	if (_SystemN > 1)
	{
		udisb(popt);
	}//加个ISB；
	
	if (_MathModel == 4) 
	{
		udifb(popt,prn_);
	}
	udclk(popt);
	prn_.clear();
	Wavelength1.clear();
	Wavelength2.clear();
	Wavelength3.clear();
	Wavelength4.clear();
	Wavelength5.clear();
	BIAS_.clear();
	BIAS_1.clear();
	BIAS_2.clear();
	BIAS_3.clear();
	BIAS_4.clear();
	X_Trans.Clear();
	P_Trans.Clear();
	//X_zhy.Output_to_File("E:\\x_zhy",1);
	//P_zhy.Output_to_File("E:\\P_zhy", 1);
	return 0;
}

int LeastSquarePPP::Udpos(ppp_option_t* popt)
{
	if (popt->mode_ppp == 1)
	{
		if (flag_epoch_count == 1)
		{
			for (int i = 1; i < 4; i++)
			{
				//X_zhy.SetValue(i, 1,popt->StaPos_CRD[i-1] );
				X_zhy.SetValue(i, 1, popt->ProXYZ[i - 1]);   //O文件头概略坐标初始化；
				P_zhy.SetValue(i, i, VAR_POS);
			}
		}	
		else if(flag_epoch_count > 1)
		{
			for (int i = 1; i < 4; i++)
			{
				double   a = X_Trans(i, 1);
				X_zhy.SetValue(i, 1, X_Trans(i,1));
				for (int j = 1; j < 4; j++)
				{
					P_zhy.SetValue(i, j, P_Trans(i, j));
				}
			}
		}
	}
	else if (popt->mode_ppp == 0)
	{
		for (int i = 1; i < 4; i++)
		{
			X_zhy.SetValue(i, 1, popt->StaPos_SPP[i - 1]);
			P_zhy.SetValue(i, i, VAR_POS);
		}
	}
	else return 0;
	return 0;
}

int LeastSquarePPP::udclk(ppp_option_t* popt)//ISB 没加  多系统时候在加把；
{

	if (flag_epoch_count==1)
	{
#if 0
		//double a = VAR_CLK;
		X_zhy.SetValue(3 * (est_xyz)+est_rclk, 1, 0);
		P_zhy.SetValue(3 * (est_xyz)+1, 3 * (est_xyz)+1, VAR_CLK);
#endif // 0

		
#if 1
		X_zhy.SetValue(3 * (est_xyz)+est_rclk, 1, (popt->Rclk));
		P_zhy.SetValue(3 * (est_xyz)+est_rclk, 3 * (est_xyz)+est_rclk, SQR(60));
#endif // 0

		

	}
	else
	{
		
		//X_zhy.SetValue(4, 1, (popt->Rclk));
		//P_zhy.SetValue(4, 4, VAR_CLK);
		//P_zhy.SetValue(4, 4, P_Trans(4, 4)+3600);
		//X_zhy.SetValue(4, 1, X_Trans(4, 1));
#if 0
		X_zhy.SetValue(3 * (est_xyz)+est_rclk, 1, 0);
		P_zhy.SetValue(3 * (est_xyz)+est_rclk, 3 * (est_xyz)+est_rclk, VAR_CLK);
#endif
#if 1
		X_zhy.SetValue(3 * (est_xyz)+est_rclk, 1, (popt->Rclk));
		P_zhy.SetValue(3 * (est_xyz)+est_rclk, 3 * (est_xyz)+est_rclk, SQR(60));
#endif 
		for (int i = 1; i < P_zhy.m_row+1; i++)
		{
			for (int j = 1; j < P_zhy.m_col+1; j++)
			{
				if ((i== 4||j==4)&&(i!=j))
				{
					P_zhy.SetValue(i, j, 0);
				}
			}
		}	
 	}	
	return 0;
}
int LeastSquarePPP::udisb(ppp_option_t* popt)
{
	int index_isb = 3 * (est_xyz) + est_rclk + _ISBNum;
	double varisb = SQR(60); double isb_rw = 0.05;
	if (flag_epoch_count == 1)
	{
		for (int i = 0; i < _ISBNum; i++)
		{
			double a = popt->ISB[i];
			/*X_zhy.SetValue(5 + i, 1, popt->ISB[i]);
			P_zhy.SetValue(5 + i, 5 + i, VAR_CLK);*/
			//X_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 1,0);
			X_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 1, popt->ISB[i]);
			P_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 3 * (est_xyz)+1 + est_rclk + i, varisb);
			if (_ISBMODE == 2)  //随机游走
			{
				X_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 1, popt->ISB[i]);
				P_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 3 * (est_xyz)+1 + est_rclk + i, SQR(20.0));
			}
		}
	}
	else
	{
		for (int i = 0; i < _ISBNum; i++)     //  有瑕疵    
		{
			//X_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 1, X_Trans(3 * (est_xyz)+1 + est_rclk + i,1));
			if (_ISBMODE == 0)//________________________________________________________________________________//白噪声
			{
				//X_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 1, 0);
				X_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 1, popt->ISB[i]);
				P_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 3 * (est_xyz)+1 + est_rclk + i, varisb);
				for (int i = 1; i < P_zhy.m_row + 1; i++)
				{
					for (int j = 1; j < P_zhy.m_col + 1; j++)
					{
						if ((i == 5 || j == 5)&&(j!=i))
						{
							P_zhy.SetValue(i, j, 0);
						}
					}
				}
			}
			else if (_ISBMODE == 1)//___________________________________________________________________________//时间常数   多长时间一个    
			{
				if (flag_epoch_count % 120 == 0)//代表被整除  
				{
					varisb = P_Trans(3 * (est_xyz)+1 + est_rclk + i, 3 * (est_xyz)+1 + est_rclk + i) * 100;
				}
				else
				{
					varisb = P_Trans(3 * (est_xyz)+1 + est_rclk + i, 3 * (est_xyz)+1 + est_rclk + i);
				}
				X_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 1, X_Trans(3 * (est_xyz)+1 + est_rclk + i, 1));
				P_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 3 * (est_xyz)+1 + est_rclk + i, varisb);
				for (int i = 1; i < index_isb; i++)  //有瑕疵   要是ISB >1的话那就是有问题了  少了协方差
				{
					P_zhy.SetValue(index_isb, i, P_Trans(index_isb, i));
					P_zhy.SetValue(i, index_isb, P_Trans(i, index_isb));
				}
			}
			else if (_ISBMODE == 2)//____________________________________________________________________________//随机游走
			{
				X_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 1, X_Trans(3 * (est_xyz)+1 + est_rclk + i, 1));
				P_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 3 * (est_xyz)+1 + est_rclk + i, (P_Trans(3 * (est_xyz)+1 + est_rclk + i, 3 * (est_xyz)+1 + est_rclk + i) + SQR(isb_rw) * popt->dt));
				for (int i = 1; i < index_isb; i++)  //有瑕疵   要是ISB >1的话那就是有问题了  少了协方差
				{
					P_zhy.SetValue(index_isb, i, P_Trans(index_isb, i));
					P_zhy.SetValue(i, index_isb, P_Trans(i, index_isb));
				}
			}

			//P_zhy.SetValue(3 * (est_xyz)+1 + est_rclk + i, 3 * (est_xyz)+1 + est_rclk + i, (P_Trans(3 * (est_xyz)+1 + est_rclk + i, 3 * (est_xyz)+1 + est_rclk + i)+0.3*0.3));

			//X_zhy.SetValue(5 + i, 1, X_Trans(5 + i, 1));
			//P_zhy.SetValue(5 + i, 5 + i, P_Trans(5 + i, 5 + i) + 100);
			/*X_zhy.SetValue(5 + i, 1, popt->ISB[i]);
			  P_zhy.SetValue(5 + i, 5 + i, VAR_CLK);*/
		}
		
	}
	return 0;
}

int LeastSquarePPP::udifb(ppp_option_t* popt,vector<int>& prn_)
{
	int num = popt->IFB_num;
	if (0)                  //白噪声
	{
		for (int i = 0; i < num; i++)
		{
			X_zhy.SetValue(_ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, 1, 0);
			P_zhy.SetValue(_ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, _ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, SQR(60));
		}
	}
	if (1)                  //常参数
	{
		if (flag_epoch_count==1)
		{
			for (int i = 0; i < num; i++)
			{
				X_zhy.SetValue(_ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, 1, 0);
				P_zhy.SetValue(_ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, _ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, SQR(60));
			}
		}
		else
		{
			for (int i = 0; i < num; i++)
			{
				for (int j = 0; j < num; j++)
				{
					X_zhy.SetValue(_ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, 1, X_Trans(_ParaN_const_pre + popt->PRE_num + _Prn_previous.size() + i + 1, 1));
					P_zhy.SetValue(_ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, _ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + j + 1, P_Trans(_ParaN_const_pre + popt->PRE_num + _Prn_previous.size() + i + 1, _ParaN_const_pre + popt->PRE_num + _Prn_previous.size() + j + 1));
				}
				
			}
			if (0)
			{
				for (int i = 0; i < _ParaN_const; i++)
				{
					for (int j = 0; j < popt->IFB_num; j++)
					{
						P_zhy.SetValue(_ParaN_const + popt->CUR_mum + m_SatValidN + j + 1, i + 1, P_Trans(_ParaN_const + popt->PRE_num + _Prn_previous.size() + j + 1, i + 1));
						P_zhy.SetValue(i + 1, _ParaN_const + popt->CUR_mum + m_SatValidN + j + 1, P_Trans(_ParaN_const + popt->PRE_num + _Prn_previous.size() + j + 1, i + 1));
					}
				}
			}
			

		}
	}
	if (0)                  //随机游走
	{
		if (flag_epoch_count == 1)
		{
			for (int i = 0; i < num; i++)
			{
				X_zhy.SetValue(_ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, 1, 0);
				P_zhy(_ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, _ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1)=SQR(60);
			}
		}
		else
		{
			for (int i = 0; i < num; i++)
			{
			
				X_zhy.SetValue(_ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, 1, X_Trans(_ParaN_const_pre + popt->PRE_num + _Prn_previous.size() + i + 1, 1));
				P_zhy.SetValue(_ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, _ParaN_const + m_SatValidN * _IonNum + popt->CUR_mum + i + 1, P_Trans(_ParaN_const_pre + popt->PRE_num + _Prn_previous.size() + i + 1, _ParaN_const_pre + popt->PRE_num + _Prn_previous.size() + i + 1)+30*SQR(0.01));
			}
			/*for (int i = 0; i < _ParaN_const; i++)
			{
				for (int j = 0; j < popt->IFB_num; j++)
				{
					P_zhy.SetValue(_ParaN_const + popt->CUR_mum + m_SatValidN + j + 1, i + 1, P_Trans(_ParaN_const + popt->PRE_num + _Prn_previous.size() + j + 1, i + 1));
					P_zhy.SetValue(i + 1, _ParaN_const + popt->CUR_mum + m_SatValidN + j + 1, P_Trans(_ParaN_const + popt->PRE_num + _Prn_previous.size() + j + 1, i + 1));
				}
			}*/
		}

	}


	return 0;
}

int LeastSquarePPP::udion(ppp_option_t* popt, unsigned int satN, vector<int>& prn_)//电离层初始值怎么给？  直接给零？？？
{	
	if (1)
	{
		if (flag_epoch_count ==1||_Prn_previous.size()==0)
		{
				for (int i = 0; i < satN; i++)//第一个历元直接赋初值 电离层 这个初值不准  格网 或者  再看看论文  
				{
					X_zhy.SetValue(_ParaN_const + est_ion + i, 1, ION_[i] );
					P_zhy.SetValue(_ParaN_const + est_ion + i, _ParaN_const + est_ion + i, VAR_IONO);
				}

				//for (int i = 0; i < satN; i++)//白噪声
				//{
				//	X_zhy.SetValue(_ParaN_const + est_ion + i, 1, 0);
				//	P_zhy.SetValue(_ParaN_const + est_ion + i, _ParaN_const + est_ion + i, VAR_IONO);
				//}
		}
		else
		{
			vector<unsigned int> Prn_Common, Prn_index_previous, Prn_index_current, bias_index_current;
			for (int i = 0; i < m_SatValidN; ++i)
			{
				for (int j = 0; j < _Prn_previous.size(); ++j)//设置前后历元观测值index  为了操作矩阵元素；
				{
					if (prn_[i] == _Prn_previous[j])
					{
						int k = 0;
						for (k = 0; k < _Cycle_Slip.size(); ++k)
						{
							if (prn_[i] == _Cycle_Slip[k])
								break;
						}
						if (k == _Cycle_Slip.size())
						{
							Prn_Common.push_back(prn_[i]);
							Prn_index_previous.push_back(j);
							Prn_index_current.push_back(i);
						}
					}
				}
			}
			_CommSatN = Prn_Common.size();
			for (int i = 0; i < m_SatValidN; ++i)
			{
				for (int j = 0; j < Prn_Common.size(); j++)
				{
					if (prn_[i] == Prn_Common[j])
					{
						break;
					}
					else if (j + 1 == Prn_Common.size())
					{
						bias_index_current.push_back(i);
					}
				}
			}
			if (bias_index_current.size() > 0)
			{
				int numb = bias_index_current.size();
				for (int k = 0; k < numb; k++)
				{
					if (X_zhy(_ParaN_const + bias_index_current[k] + 1, 1) == 0)
					{
						X_zhy(_ParaN_const + bias_index_current[k] + 1, 1) = ION_[bias_index_current[k]] ;
						P_zhy(_ParaN_const + bias_index_current[k] + 1, _ParaN_const + bias_index_current[k] + 1) = VAR_IONO;
					}
				}
			}
			for (int i = 0; i < _CommSatN; i++)
			{
				for (int j = 0; j < _AmbNum; j++)
				{
					X_zhy(_ParaN_const + _AmbNum * Prn_index_current[i] + j + 1, 1) = X_Trans(_ParaN_const_pre + _AmbNum * Prn_index_previous[i] + j + 1, 1);
					for (int k = 0; k < _CommSatN; k++)
					{

						P_zhy(_ParaN_const + _AmbNum * Prn_index_current[i] + j + 1, _ParaN_const + _AmbNum * Prn_index_current[k] + j + 1) =
							P_Trans(_ParaN_const_pre + _AmbNum * Prn_index_previous[i] + j + 1, _ParaN_const_pre + _AmbNum * Prn_index_previous[k] + j + 1);
						if (_ParaN_const + _AmbNum * Prn_index_current[i] + j + 1 == _ParaN_const + _AmbNum * Prn_index_current[k] + j + 1)
						{
							double psd =  0.04;
							//
							//{
#if 1
							if (1) {
								if (PSD[Prn_Common[i]] != 0)
								{
									if (PSD[Prn_Common[i]]< 0.0001)
									{
										//psd = PSD[Prn_Common[i]] * 100;
										psd = PSD[Prn_Common[i]] * 10;
									}
									else if(PSD[Prn_Common[i]] <0.004)
									{
										//psd = PSD[Prn_Common[i]] * 10;
										psd = PSD[Prn_Common[i]];
									}
									else
									{
										psd = PSD[Prn_Common[i]];
									}
								}
							}
							/*else
							{

								if (PSD[Prn_Common[i]] != 0)
								{
									if (PSD[Prn_Common[i]] > 0.04)
									{
										psd = 0.04;
									}
									else if (PSD[Prn_Common[i]] < 0.00001)
									{
										psd = PSD[Prn_Common[i]] * 100;
									}
									else
									{
										psd = PSD[Prn_Common[i]] * 10;
									}
								}
							}*/
#endif // 0

							
								
							//}
																																	//0.04    0.02  0.01 0.001
							if (_MathModel == 4)
							{
								P_zhy(_ParaN_const + _AmbNum * Prn_index_current[i] + j + 1, _ParaN_const + _AmbNum * Prn_index_current[k] + j + 1) += SQR(psd) * popt->interval;
							}//电离层随机游走
							else{ P_zhy(_ParaN_const + _AmbNum * Prn_index_current[i] + j + 1, _ParaN_const + _AmbNum * Prn_index_current[k] + j + 1) += SQR(_IonProcessNoise) * popt->interval; }
						}
					}
				}
			}

			//for (int i = 0; i < satN; i++)//白噪声
			//{
			//	X_zhy.SetValue(_ParaN_const + est_ion + i, 1, 0);
			//	P_zhy.SetValue(_ParaN_const + est_ion + i, _ParaN_const + est_ion + i, VAR_IONO);
			//}


			// Covariance of X,Y,Z,Trop with ION
			for (int i = 0; i < _CommSatN; ++i)
			{
				for (int j = 0; j < _ParaN_const; ++j)
				{
					for (int k = 0; k < _AmbNum; ++k)
					{
						P_zhy(_ParaN_const + _AmbNum * Prn_index_current[i] + 1, j + 1) =
							P_Trans(_ParaN_const_pre + _AmbNum * Prn_index_previous[i] + k + 1, j + 1);
					}
				}
			}
			for (int i = 0; i < _ParaN_const; ++i)
			{
				for (int j = 0; j < _CommSatN; ++j)
				{
					for (int k = 0; k < _AmbNum; ++k)
					{
						P_zhy(i + 1, _ParaN_const + _AmbNum * Prn_index_current[j] + 1) =
							P_Trans(i + 1, _ParaN_const_pre + _AmbNum * Prn_index_previous[j] + k + 1);
					}
				}
			}
			
			// Covariance of BIAS,Trop with ION
			if (_MathModel == 4)
			{
				for (int i = 0; i < _CommSatN; ++i)
				{

					for (int k = 0; k < _CommSatN; ++k)
					{
						for (int cur = 0; cur < num_zhy[Prn_index_current[k]]; cur++)
						{
							for (int pre = 0; pre < num_zhy_pre[Prn_index_previous[k]]; pre++)
							{
								if (type_[Prn_index_current[k]][cur + 1] == type_pre[Prn_index_previous[k]][pre + 1]) 
								{
									P_zhy(_ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + satN + cur, _ParaN_const + Prn_index_current[i] + 1) =
									P_Trans(_ParaN_const_pre + _AmbNum * num_zhy_accu_pre[Prn_index_previous[k]] + 1 + pre + _Prn_previous.size(), _ParaN_const_pre + Prn_index_previous[i] + 1);
									if (popt->ifb_mode==1)
									{
										for (int covifb = 0; covifb < popt->IFB_num; covifb++)
										{
											P_zhy(_ParaN_const + popt->CUR_mum + 1 + satN + covifb, _ParaN_const + Prn_index_current[i] + 1) =
												P_Trans(_ParaN_const_pre + popt->PRE_num + 1 + covifb + _Prn_previous.size(), _ParaN_const_pre + Prn_index_previous[i] + 1);
											P_zhy(_ParaN_const + Prn_index_current[i] + 1, _ParaN_const + popt->CUR_mum + 1 + satN + covifb ) =
												P_Trans(_ParaN_const_pre + popt->PRE_num + 1 + covifb + _Prn_previous.size(), _ParaN_const_pre + Prn_index_previous[i] + 1);

										}
										
									}
									//P_zhy.Output_to_File("E:\\p.txt", 0);
								}
							}
						}
						/*if (num_zhy[Prn_index_current[k]] == num_zhy_pre[Prn_index_previous[k]])
						{
							for (int t = 0; t < num_zhy[Prn_index_current[k]]; t++)
							{
								P_zhy(_ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + satN + j+t, _ParaN_const + j + Prn_index_current[i] + 1) =
									P_Trans(_ParaN_const + _AmbNum * num_zhy_accu_pre[Prn_index_previous[k]] + j + 1 + _Prn_previous.size() + t, _ParaN_const + j + Prn_index_previous[i] + 1);
							}
						}
						else if (num_zhy[Prn_index_current[k]] > num_zhy_pre[Prn_index_previous[k]])
						{
							for (int t = 0; t < num_zhy_pre[Prn_index_current[k]]; t++)
							{
								P_zhy(_ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + satN + j + t, _ParaN_const + j + Prn_index_current[i] + 1) =
									P_Trans(_ParaN_const + _AmbNum * num_zhy_accu_pre[Prn_index_previous[k]] + j + 1 + _Prn_previous.size() + t, _ParaN_const + j + Prn_index_previous[i] + 1);
							}
						}
						else if (num_zhy[Prn_index_current[k]] < num_zhy_pre[Prn_index_previous[k]])
						{
							for (int t = 0; t < num_zhy[Prn_index_current[k]]; t++)
							{
								P_zhy(_ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + satN + j + t, _ParaN_const + j + Prn_index_current[i] + 1) =
									P_Trans(_ParaN_const + _AmbNum * num_zhy_accu_pre[Prn_index_previous[k]] + j + 1 + _Prn_previous.size() + t, _ParaN_const + j + Prn_index_previous[i] + 1);
							}
						}*/
					}

					/*for (int i = 0; i < length; i++)
					{

					}*/




				}
				for (int i = 0; i < _CommSatN; ++i)
				{
					for (int j = 0; j < _AmbNum; ++j)
					{
						for (int k = 0; k < _CommSatN; ++k)
						{
							for (int cur = 0; cur < num_zhy[Prn_index_current[k]]; cur++)
							{
								for (int pre = 0; pre < num_zhy_pre[Prn_index_previous[k]]; pre++)
								{
									if (type_[Prn_index_current[k]][cur+1] == type_pre[Prn_index_previous[k]][pre+1])
									{

										P_zhy(_ParaN_const + j + Prn_index_current[i] + 1, _ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + satN + j + cur) =
											P_Trans(_ParaN_const_pre + j + Prn_index_previous[i] + 1, _ParaN_const_pre + _AmbNum * num_zhy_accu_pre[Prn_index_previous[k]] + j + 1 + _Prn_previous.size() + pre);

									}
								}
							}
							/*if (num_zhy[Prn_index_current[k]] == num_zhy_pre[Prn_index_previous[k]])
							{
								for (int t = 0; t < num_zhy[Prn_index_current[k]]; t++)
								{
									P_zhy(_ParaN_const + j + Prn_index_current[i] + 1, _ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + satN + j + t) =
										P_Trans(_ParaN_const + j + Prn_index_previous[i] + 1, _ParaN_const + _AmbNum * num_zhy_accu_pre[Prn_index_previous[k]] + j + 1 + _Prn_previous.size() + t);
								}
							}
							else if (num_zhy[Prn_index_current[k]] > num_zhy_pre[Prn_index_previous[k]])
							{
								for (int t = 0; t < num_zhy_pre[Prn_index_current[k]]; t++)
								{
									P_zhy(_ParaN_const + j + Prn_index_current[i] + 1, _ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + satN + j + t) =
										P_Trans(_ParaN_const + j + Prn_index_previous[i] + 1, _ParaN_const + _AmbNum * num_zhy_accu_pre[Prn_index_previous[k]] + j + 1 + _Prn_previous.size() + t);
								}
							}
							else if (num_zhy[Prn_index_current[k]] < num_zhy_pre[Prn_index_previous[k]])
							{
								for (int t = 0; t < num_zhy[Prn_index_current[k]]; t++)
								{
									P_zhy(_ParaN_const + j + Prn_index_current[i] + 1, _ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + satN + j + t) =
										P_Trans(_ParaN_const + j + Prn_index_previous[i] + 1, _ParaN_const + _AmbNum * num_zhy_accu_pre[Prn_index_previous[k]] + j + 1 + _Prn_previous.size() + t);
								}
							}*/
						}
					}
				}
			}
			else
			{
				for (int i = 0; i < _CommSatN; ++i)
				{
					for (int j = 0; j < _AmbNum; ++j)
					{
						for (int k = 0; k < _CommSatN; ++k)
						{
							P_zhy(_ParaN_const + _AmbNum * Prn_index_current[i] + 1 + satN + j, _ParaN_const + j + Prn_index_current[k] + 1) =
								P_Trans(_ParaN_const_pre + _AmbNum * Prn_index_previous[i] + j + 1 + _Prn_previous.size(), _ParaN_const_pre + j + Prn_index_previous[k] + 1);
						}
					}
				}
				for (int i = 0; i < _CommSatN; ++i)
				{
					for (int j = 0; j < _AmbNum; ++j)
					{
						for (int k = 0; k < _CommSatN; ++k)
						{
							P_zhy(_ParaN_const + j + Prn_index_current[k] + 1, _ParaN_const + _AmbNum * Prn_index_current[i] + 1 + satN + j) =
								P_Trans(_ParaN_const_pre + j + Prn_index_previous[k] + 1, _ParaN_const_pre + _AmbNum * Prn_index_previous[i] + j + 1 + _Prn_previous.size());
						}
					}
				}
			}
			
		}
	}
	return 0;
}
int LeastSquarePPP::udtrop(ppp_option_t* popt)//对流层梯度估计 没写上 写上是不是参数太多？？？ 方程病态？
{
	int dtrpindex = 3 * est_xyz + est_rclk * _SystemN + est_trop;
	int dtrpindex_pre= 3 * est_xyz + est_rclk * _SystemN_pre + est_trop;
	if (flag_epoch_count == 1 || _Prn_previous.size() == 0)
	{
		X_zhy.SetValue(dtrpindex, 1, 0.15);//0.15为经验初值；
		P_zhy.SetValue(dtrpindex, dtrpindex, VAR_ZTD);
	}
	if (flag_epoch_count>1)
	{
		X_zhy(dtrpindex, 1) = X_Trans(dtrpindex, 1);
		double dtrp=P_Trans(dtrpindex_pre, dtrpindex_pre) + SQR(_TropProcessNoise) * popt->interval;
		P_zhy.SetValue(dtrpindex, dtrpindex, dtrp);
		for (int i = 1; i < 3 * est_xyz+1; i++)
		{
			P_zhy.SetValue(i, dtrpindex, P_Trans(i, dtrpindex_pre));
			P_zhy.SetValue(dtrpindex, i, P_Trans(dtrpindex_pre, i));
		}
	}
	return 0;
}

int LeastSquarePPP::udbias(ppp_option_t* popt, unsigned int satN_,vector<int>& prn_)
{
	
	if (_MathModel==4)
	{
		if (flag_epoch_count == 1||_Prn_previous.size() == 0)
		{
			for (int i = 0; i < satN_; i++)
			{
				for (int j = 0; j < num_zhy[i]; j++)
				{
					if (type_int[i][j] == 1)X_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[i] + j + 1, 1) = BIAS_[i];
					else if (type_int[i][j] == 2)X_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[i] + j + 1, 1) = BIAS_1[i];
					else if (type_int[i][j] == 3)
					{
						X_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[i] + j + 1, 1) = BIAS_2[i];
					}
					else if (type_int[i][j] == 4)
					{
						X_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[i] + j + 1, 1) = BIAS_3[i];
					}
					else if (type_int[i][j] == 5)
					{
						X_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[i] + j + 1, 1) = BIAS_4[i];
					}

					P_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[i] + j + 1, _ParaN_const + satN_ * _IonNum + num_zhy_accu[i] + j + 1) = VAR_BIAS;

				}
			}
		}
		else
		{
			vector<unsigned int> Prn_Common, Prn_index_previous, Prn_index_current, bias_index_current;
			for (int i = 0; i < m_SatValidN; ++i)
			{
				for (int j = 0; j < _Prn_previous.size(); ++j)
				{
					if (prn_[i] == _Prn_previous[j])
					{
						int k = 0;
						for (k = 0; k < _Cycle_Slip.size(); ++k)
						{
							if (prn_[i] == _Cycle_Slip[k])
								break;
						}
						if (k == _Cycle_Slip.size())
						{
							Prn_Common.push_back(prn_[i]);
							Prn_index_previous.push_back(j);
							Prn_index_current.push_back(i);
						}
					}
				}
			}
			for (int i = 0; i < m_SatValidN; ++i)
			{
				for (int j = 0; j < Prn_Common.size(); j++)
				{
					if (prn_[i] == Prn_Common[j])
					{
						break;
					}
					else if (j + 1 == Prn_Common.size())
					{
						bias_index_current.push_back(i);
					}
				}
			}
			
			_CommSatN = Prn_Common.size();

			for (int i = 0; i < _CommSatN; i++)
			{
				for (int j = 0; j < _AmbNum; j++)
				{
					for (int cur = 0; cur < num_zhy[Prn_index_current[i]]; cur++)
					{
						int pre = 0;
						for (pre = 0; pre < num_zhy_pre[Prn_index_previous[i]]; pre++)
						{
							if (type_[Prn_index_current[i]][cur + 1] == type_pre[Prn_index_previous[i]][pre + 1])
							{
								X_zhy(_ParaN_const + m_SatValidN * _IonNum + _AmbNum * num_zhy_accu[Prn_index_current[i]] + j + 1 + cur, 1) = X_Trans(_ParaN_const_pre + _Prn_previous.size() * _IonNum + _AmbNum * num_zhy_accu_pre[Prn_index_previous[i]] + j + 1 + pre, 1);
								break;
							} 
						}
						if (pre == num_zhy_pre[Prn_index_previous[i]])
						{
							if (type_int[Prn_index_current[i]][cur] == 1) 
							{
								X_zhy(_ParaN_const + m_SatValidN * _IonNum + _AmbNum * num_zhy_accu[Prn_index_current[i]] + j + 1 + cur, 1) = BIAS_[Prn_index_current[i]];
							}
							else if (type_int[Prn_index_current[i]][cur] == 2)
							{
								X_zhy(_ParaN_const + m_SatValidN * _IonNum + _AmbNum * num_zhy_accu[Prn_index_current[i]] + j + 1 + cur, 1) = BIAS_1[Prn_index_current[i]];
							}
							else if (type_int[Prn_index_current[i]][cur] == 3)
							{
								X_zhy(_ParaN_const + m_SatValidN * _IonNum + _AmbNum * num_zhy_accu[Prn_index_current[i]] + j + 1 + cur, 1) = BIAS_2[Prn_index_current[i]];
							}
							else if (type_int[Prn_index_current[i]][cur] == 4)
							{
								X_zhy(_ParaN_const + m_SatValidN * _IonNum + _AmbNum * num_zhy_accu[Prn_index_current[i]] + j + 1 + cur, 1) = BIAS_3[Prn_index_current[i]];
							}
							else if (type_int[Prn_index_current[i]][cur] == 5)
							{
								X_zhy(_ParaN_const + m_SatValidN * _IonNum + _AmbNum * num_zhy_accu[Prn_index_current[i]] + j + 1 + cur, 1) = BIAS_4[Prn_index_current[i]];
							}
						}
					}
				}
			}
			for (int i = 0; i < _CommSatN; i++){
				for (int cur = 0; cur < num_zhy[Prn_index_current[i]]; cur++)
				{
					for (int pre = 0; pre < num_zhy_pre[Prn_index_previous[i]]; pre++)
					{
						if (type_[Prn_index_current[i]][cur + 1] == type_pre[Prn_index_previous[i]][pre + 1])
						{
							for (int k = 0; k < _CommSatN; k++) {
								for (int cur1 = 0; cur1 < num_zhy[Prn_index_current[k]]; cur1++)
								{
									int pre1 = 0;
									for ( pre1 = 0; pre1 < num_zhy_pre[Prn_index_previous[k]]; pre1++)
									{
										if (type_[Prn_index_current[k]][cur1 + 1] == type_pre[Prn_index_previous[k]][pre1 + 1])
										{
											P_zhy(_ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[i]] + cur + 1 + m_SatValidN * _IonNum, _ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + cur1 + m_SatValidN * _IonNum) =
												P_Trans(_ParaN_const_pre + _AmbNum * num_zhy_accu_pre[Prn_index_previous[i]] + 1 + pre + _Prn_previous.size() * _IonNum, _ParaN_const_pre + _AmbNum * num_zhy_accu_pre[Prn_index_previous[k]] + pre1 + 1 + _Prn_previous.size() * _IonNum);
											if (popt->ifb_mode==1)
											{
												for (int covifb = 0; covifb < popt->IFB_num; covifb++)
												{
													P_zhy(_ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[i]] + cur + 1 + m_SatValidN * _IonNum, _ParaN_const + popt->CUR_mum + 1  + covifb + m_SatValidN * _IonNum) =
														P_Trans(_ParaN_const_pre + _AmbNum * num_zhy_accu_pre[Prn_index_previous[i]] + 1 + pre + _Prn_previous.size() * _IonNum, _ParaN_const_pre + popt->PRE_num  + covifb + 1 + _Prn_previous.size() * _IonNum);
												}
												for (int covifb = 0; covifb < popt->IFB_num; covifb++)
												{
													P_zhy(_ParaN_const + popt->CUR_mum + covifb  + 1 + m_SatValidN * _IonNum, _ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + cur1 + m_SatValidN * _IonNum) =
														P_Trans(_ParaN_const_pre + popt->PRE_num + covifb + 1 + _Prn_previous.size() * _IonNum, _ParaN_const_pre + _AmbNum * num_zhy_accu_pre[Prn_index_previous[k]] + pre1 + 1 + _Prn_previous.size() * _IonNum);
												}
											}
											break;
										}
									}                                                                                           
									if (pre1 == num_zhy_pre[Prn_index_previous[k]])
									{
										P_zhy(_ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + cur1 + m_SatValidN * _IonNum, _ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[k]] + 1 + cur1 + m_SatValidN * _IonNum) = VAR_BIAS;
									}
								}
							}
						}
					}
				}
			}
			
			// Covariance of X,Y,Z,Trop with Ambituities
			for (int i = 0; i < _CommSatN; ++i)
			{
				for (int k = 0; k < _ParaN_const; ++k)
				{
					for (int j = 0; j < _AmbNum; ++j)
					{
						for (int cur = 0; cur < num_zhy[Prn_index_current[i]]; cur++)
						{
							for (int pre = 0; pre < num_zhy_pre[Prn_index_previous[i]]; pre++)
							{
								if (type_[Prn_index_current[i]][cur+1] == type_pre[Prn_index_previous[i]][pre+1])
								{
									P_zhy(_ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[i]] + 1 + m_SatValidN + j + cur, k + 1) =
										P_Trans(_ParaN_const_pre + _AmbNum * num_zhy_accu_pre[Prn_index_previous[i]] + j + 1 + _Prn_previous.size() + pre, k + 1);
								}
							}
						}
					
					}
				}
			}
			for (int i = 0; i < _CommSatN; ++i)
			{
				for (int k = 0; k < _ParaN_const; ++k)
				{
					for (int j = 0; j < _AmbNum; ++j)
					{
						for (int cur = 0; cur < num_zhy[Prn_index_current[i]]; cur++)
						{
							for (int pre = 0; pre < num_zhy_pre[Prn_index_previous[i]]; pre++)
							{
								if (type_[Prn_index_current[i]][cur+1] == type_pre[Prn_index_previous[i]][pre+1]) 
								{
									P_zhy(k + 1, _ParaN_const + _AmbNum * num_zhy_accu[Prn_index_current[i]] + 1 + m_SatValidN + j + cur) =
										P_Trans(k + 1, _ParaN_const_pre + _AmbNum * num_zhy_accu_pre[Prn_index_previous[i]] + j + 1 + _Prn_previous.size() + pre);
								}
							}
						}
						
					}
				}
			}
			if (bias_index_current.size() > 0)//窜行
			{
				int numb = bias_index_current.size();
				for (int k = 0; k < numb; k++)
				{
					for (int t = 0; t < num_zhy[bias_index_current[k]]; t++)
					{

						if (type_int[bias_index_current[k]][t] == 1)
						{

							X_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + 1, 1) = BIAS_[bias_index_current[k]];
							P_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + 1, _ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + 1) = VAR_BIAS;
						}
						else if (type_int[bias_index_current[k]][t] == 2)
						{

							X_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1, 1) = BIAS_1[bias_index_current[k]];
							P_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1, _ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1) = VAR_BIAS;

						}
						else if (type_int[bias_index_current[k]][t] == 3)//在这
						{
							X_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1, 1) = BIAS_2[bias_index_current[k]];
							P_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1, _ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1) = VAR_BIAS;
						}
						else if (type_int[bias_index_current[k]][t] == 4)
						{
							X_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1, 1) = BIAS_3[bias_index_current[k]];
							P_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1, _ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1) = VAR_BIAS;
						}
						else if (type_int[bias_index_current[k]][t] == 5)
						{
							X_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1, 1) = BIAS_4[bias_index_current[k]];
							P_zhy(_ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1, _ParaN_const + satN_ * _IonNum + num_zhy_accu[bias_index_current[k]] + t + 1) = VAR_BIAS;
						}
					}

				}
			}

			
			for (int zzz =0 ; zzz < popt->epochnum; zzz++)
			{
				if (P_zhy(_ParaN_const + satN_ * _IonNum + 1 + zzz, _ParaN_const + satN_ * _IonNum + 1 + zzz) == 0) 
					P_zhy(_ParaN_const + satN_ * _IonNum + 1 + zzz, _ParaN_const + satN_ * _IonNum + 1 + zzz)== VAR_BIAS;
			}

		}
	}
	else
	{
		if (flag_epoch_count == 1||_Prn_previous.size() == 0)
		{
			for (int i = 0; i < satN_; i++)
			{
				if (_MathModel == 2)
				{
					X_zhy(_ParaN_const + satN_ * _IonNum + i + 1, 1) = BIAS_[i] / 2;
				}
				else
				{
					X_zhy(_ParaN_const + satN_ * _IonNum + i + 1, 1) = BIAS_[i];
				}
				P_zhy(_ParaN_const + satN_ * _IonNum + i + 1, _ParaN_const + satN_ * _IonNum + i + 1) = VAR_BIAS;
			}
		}
		else
		{
			vector<unsigned int> Prn_Common, Prn_index_previous, Prn_index_current, bias_index_current;
			for (int i = 0; i < m_SatValidN; ++i)
			{
				for (int j = 0; j < _Prn_previous.size(); ++j)
				{
					if (prn_[i] == _Prn_previous[j])
					{
						int k = 0;
						for (k = 0; k < _Cycle_Slip.size(); ++k)
						{
							if (prn_[i] == _Cycle_Slip[k])
								break;
						}
						if (k == _Cycle_Slip.size())
						{
							Prn_Common.push_back(prn_[i]);
							Prn_index_previous.push_back(j);
							Prn_index_current.push_back(i);
						}
					}
				}
			}
			for (int i = 0; i < m_SatValidN; ++i)
			{
				for (int j = 0; j < Prn_Common.size(); j++)
				{
					if (prn_[i] == Prn_Common[j])
					{
						break;
					}
					else if (j + 1 == Prn_Common.size())
					{
						bias_index_current.push_back(i);
					}
				}
			}
			if (bias_index_current.size() > 0)
			{
				int numb = bias_index_current.size();
				for (int k = 0; k < numb; k++)
				{
					if (X_zhy(_ParaN_const + satN_ * _IonNum + bias_index_current[k] + 1, 1) == 0)
					{
						if (_MathModel == MODE_UOFC_COMBINE)
						{
							X_zhy(_ParaN_const + satN_ * _IonNum + bias_index_current[k] + 1, 1) = BIAS_[bias_index_current[k]] / 2;
							P_zhy(_ParaN_const + satN_ * _IonNum + bias_index_current[k] + 1, _ParaN_const + satN_ * _IonNum + bias_index_current[k] + 1) = VAR_BIAS;
						}
						else
						{
							X_zhy(_ParaN_const + satN_ * _IonNum + bias_index_current[k] + 1, 1) = BIAS_[bias_index_current[k]];
							P_zhy(_ParaN_const + satN_ * _IonNum + bias_index_current[k] + 1, _ParaN_const + satN_ * _IonNum + bias_index_current[k] + 1) = VAR_BIAS;
						}

					}
				}
			}
			_CommSatN = Prn_Common.size();
			for (int i = 0; i < _CommSatN; i++)
			{

				for (int j = 0; j < _AmbNum; j++)
				{
					X_zhy(_ParaN_const + m_SatValidN * _IonNum + _AmbNum * Prn_index_current[i] + j + 1, 1) = X_Trans(_ParaN_const_pre + _Prn_previous.size() * _IonNum + _AmbNum * Prn_index_previous[i] + j + 1, 1);
					for (int k = 0; k < _CommSatN; k++)
					{
						P_zhy(_ParaN_const + _AmbNum * Prn_index_current[i] + j + 1 + m_SatValidN * _IonNum, _ParaN_const + _AmbNum * Prn_index_current[k] + j + 1 + m_SatValidN * _IonNum) =
							P_Trans(_ParaN_const_pre + _AmbNum * Prn_index_previous[i] + j + 1 + _Prn_previous.size() * _IonNum, _ParaN_const_pre + _AmbNum * Prn_index_previous[k] + j + 1 + _Prn_previous.size() * _IonNum);
					}

				}
			}
			// Covariance of X,Y,Z,Trop with Ambituities
			for (int i = 0; i < _CommSatN; ++i)
			{
				for (int j = 0; j < _ParaN_const; ++j)
				{
					for (int k = 0; k < _AmbNum; ++k)
					{
						P_zhy(_ParaN_const + _AmbNum * Prn_index_current[i] + 1 + m_SatValidN * _IonNum, j + 1) =
							P_Trans(_ParaN_const_pre + _AmbNum * Prn_index_previous[i] + k + 1 + _Prn_previous.size() * _IonNum, j + 1);
					}
				}
			}

			for (int i = 0; i < _ParaN_const; ++i)
			{
				for (int j = 0; j < _CommSatN; ++j)
				{
					for (int k = 0; k < _AmbNum; ++k)
					{
						P_zhy(i + 1, _ParaN_const + _AmbNum * Prn_index_current[j] + 1 + m_SatValidN * _IonNum) =
							P_Trans(i + 1, _ParaN_const_pre + _AmbNum * Prn_index_previous[j] + 1 + _Prn_previous.size() * _IonNum);
					}
				}
			}
		}
	}
	

	return 0;
}


