#include <string.h>
#include "SUPREME_DownLoad.h"
#include "SUPREME_CommonFunction.h"

/// GNSS data download Web Site
const static char cddis_nasa_ftp[50] = "ftp://cddis.nasa.gov/pub/";
const static char product_cddis_nasa_ftp[100] = "ftp://cddis.nasa.gov/pub/gps/products/";
const static char product_cddis_nasa_gbm_ftp[100] = "ftp://cddis.nasa.gov/pub/gps/products/mgex/";
const static char obsnav_cddis_nasa_ftp[100] = "ftp://cddis.nasa.gov/pub/gps/data/daily/";
const static char obsnav_cddis_nasa_mgex_ftp[100] = "ftp://cddis.nasa.gov/pub/gps/data/campaign/mgex/daily/rinex3/";
const static char unibe_code_ftp[50] = "ftp://ftp.aiub.unibe.ch/";
const static char crd_unibe_ftp[100] = "ftp://ftp.aiub.unibe.ch/BSWUSER52/STA/";

/// Download Navigation file
int Download_Nav(char* sta, char *des_path, gnsstime &t)
{
	char nav_fname[50] = { 0 }, nav_fname_rnx2[100], nav_fname_rnx3[100];
	char download_command_nav[300] = { 0 };
	char decompression_command_nav[100] = { 0 };

	if (strstr(sta, "brdc"))
	{
		if (t.m_Year >= 2017)
		{
			sprintf(nav_fname, "BRDC00IGS_R_%04d%03d0000_01D_MN.rnx.gz", t.m_Year, t.m_Doy);
			sprintf(nav_fname_rnx3, "%sBRDC00IGS_R_%04d%03d0000_01D_MN.rnx", des_path, t.m_Year, t.m_Doy);
			sprintf(nav_fname_rnx2, "%sbrdm%03d0.%02dp", des_path, t.m_Doy, t.m_Year % 100);
		}
		
	}
	else
	{
		sprintf(nav_fname, "%s%03d0.%02d%c.Z", sta, t.m_Doy, t.m_Year % 100, 'p');
	}

	sprintf(download_command_nav, ".\\bin\\wget.exe -nv -P%s %s%4d/%03d/%02d%c/%s -t 2",
		des_path, obsnav_cddis_nasa_ftp, t.m_Year, t.m_Doy, t.m_Year % 100, 'p', nav_fname);

	//sprintf(download_command_nav, ".\\bin\\wget.exe -nv -P%s %s%4d/%03d/%02dp/%s -t 2",
	//	des_path, obsnav_cddis_nasa_ftp, t.m_Year, t.m_Doy, t.m_Year % 100, nav_fname);

	sprintf(decompression_command_nav, ".\\bin\\gzip.exe -N -d %s%s", des_path, nav_fname);

	system(download_command_nav);

	system(decompression_command_nav);

	rename(nav_fname_rnx3, nav_fname_rnx2);

	return 1;
}

/// Download Observation file
int Download_Obs(char* sta, char *des_path, gnsstime &t, unsigned int fnamemode)
{
	char obs_fname[50] = { 0 }, obs_fname_rnx2[MAX_DOWNLOAD_COMMAND_LENGTH], obs_fname_rnx3[MAX_DOWNLOAD_COMMAND_LENGTH];
	char obs_fname_rnx2o[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char download_command_obs[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char decompression_command_obs[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char d2o_command[100] = { 0 };

	char sta_name[5] = { 0 };
	char country[4] = { 0 };

	
	sprintf(obs_fname_rnx2, "%s%s%03d0.%02dd", des_path,sta, t.m_Doy, t.m_Year % 100);
	sprintf(obs_fname_rnx2o, "%s%s%03d0.%02do", des_path, sta, t.m_Doy, t.m_Year % 100);

	if (FileExist(obs_fname_rnx2o))
		return 1;
	
	if (fnamemode == _LONGFNAME)
	{
		char flag = 0;
		if (!GetCountry("conf/Station_MGEX_List.dat",sta, country, flag))
		{
			printf("Warrning: Can not find country information in the list.\n");
			return 0;
		}

		sprintf(obs_fname, "%s00%s_%c_%04d%03d0000_01D_30S_MO.crx.gz", sta, country,flag, t.m_Year, t.m_Doy);	
		sprintf(obs_fname_rnx3, "%s%s00%s_R_%04d%03d0000_01D_30S_MO.crx", des_path, sta, country, t.m_Year, t.m_Doy, t.m_Year % 100);
	}

	if (fnamemode == _SHORTFNAME)
	{
		Cap2Lower(sta, sta_name);
		sprintf(obs_fname, "%s%03d0.%02dd.Z", sta_name, t.m_Doy, t.m_Year % 100);
	}

	sprintf(download_command_obs, ".\\bin\\wget.exe -nv -P%s %s%4d/%03d/%02dd/%s -t 2",
		des_path, obsnav_cddis_nasa_ftp, t.m_Year, t.m_Doy, t.m_Year % 100, obs_fname);

	sprintf(decompression_command_obs, ".\\bin\\gzip.exe -N -d %s%s", des_path, obs_fname);

	system(download_command_obs);

	system(decompression_command_obs);

	if (fnamemode == _LONGFNAME)
		rename(obs_fname_rnx3, obs_fname_rnx2);

	Crx2Rnx(obs_fname_rnx2, 1);

	return 1;
}

/// Download precise ephemeris
int DownloadEph(char *des_path, char *Center, gnsstime &t)
{
	char ephf_name[20] = { 0 }, ephf2_name[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char download_command_eph[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char decompression_command_eph[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };

	if (!FileExist(des_path))
		CreatePath(des_path);

	if (strcmp(Center, "cod") == 0 || strcmp(Center, "COD") == 0
		|| strcmp(Center, "cof") == 0 || strcmp(Center, "COF") == 0)
	{
		sprintf(ephf_name, "%s%4d%1d.eph.Z", Center, t.m_Week, t.m_Dow);
		sprintf(ephf2_name, "%s%s%4d%1d.eph",des_path, Center, t.m_Week, t.m_Dow);
	}
	else
	{
		sprintf(ephf_name, "%s%4d%1d.sp3.Z", Center, t.m_Week, t.m_Dow);
		sprintf(ephf2_name, "%s%s%4d%1d.sp3",des_path, Center, t.m_Week, t.m_Dow);
	}

	if (FileExist(ephf2_name))
		return 0;

	// Multi-GNSS product
	if (strcmp(Center, "gbm") == 0 || strcmp(Center, "GBM") == 0
		|| strcmp(Center, "wum") == 0 || strcmp(Center, "WUM") == 0
		|| strcmp(Center, "com") == 0 || strcmp(Center, "COM") == 0
		|| strcmp(Center, "grm") == 0 || strcmp(Center, "GRM") == 0
		|| strcmp(Center, "tum") == 0 || strcmp(Center, "TUM") == 0)
		sprintf(download_command_eph, ".\\bin\\wget.exe -nv -P%s %s%4d/%s -t 2",
		des_path, product_cddis_nasa_gbm_ftp, t.m_Week, ephf_name);
	else // GPS and GLONASS product
		sprintf(download_command_eph, ".\\bin\\wget.exe -nv -P%s %s%4d/%s -t 2",
		des_path, product_cddis_nasa_ftp, t.m_Week, ephf_name);

	sprintf(decompression_command_eph, ".\\bin\\gzip.exe -N -d %s%s", des_path, ephf_name);

	system(download_command_eph);
	system(decompression_command_eph);

	return 1;
}

/// Download precise clocks
int DownloadClk(char *des_path, char *Center, gnsstime &t)
{
	char clkf_name[20] = { 0 }, clkf2_name[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char download_command_clk[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char decompression_command_clk[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };

	if (strcmp(Center, "cod") == 0 || strcmp(Center, "COD") == 0
		|| strcmp(Center, "cof") == 0 || strcmp(Center, "COF") == 0)
	{
		sprintf(clkf_name, "%s%4d%1d.clk.Z", Center, t.m_Week, t.m_Dow);
		sprintf(clkf2_name, "%s%s%4d%1d.clk",des_path, Center, t.m_Week, t.m_Dow);
	}
	else
	{
		sprintf(clkf_name, "%s%4d%1d.clk.Z", Center, t.m_Week, t.m_Dow);
		sprintf(clkf2_name, "%s%s%4d%1d.clk",des_path, Center, t.m_Week, t.m_Dow);
	}

	if (FileExist(clkf2_name))
		return 0;

	// Multi-GNSS product
	if (strcmp(Center, "gbm") == 0 || strcmp(Center, "GBM") == 0 
		|| strcmp(Center, "wum") == 0 || strcmp(Center, "WUM") == 0
		|| strcmp(Center, "com") == 0 || strcmp(Center, "COM") == 0
		|| strcmp(Center, "grm") == 0 || strcmp(Center, "GRM") == 0
		|| strcmp(Center, "tum") == 0 || strcmp(Center, "TUM") == 0)
		sprintf(download_command_clk, ".\\bin\\wget.exe -nv -P%s %s%4d/%s -t 2",
		des_path, product_cddis_nasa_gbm_ftp, t.m_Week, clkf_name);
	else // GPS and GLONASS product
		sprintf(download_command_clk, ".\\bin\\wget.exe -nv -P%s %s%4d/%s -t 2",
		des_path, product_cddis_nasa_ftp, t.m_Week, clkf_name);

	sprintf(decompression_command_clk, ".\\bin\\gzip.exe -N -d %s%s", des_path, clkf_name);

	system(download_command_clk);
	system(decompression_command_clk);

	return 1;
}

int Download_EphClk(char *des_path, char *Center, gnsstime &t, unsigned int typeflag)
{
	if (typeflag == EPH || typeflag == EPHCLK)
		DownloadEph(des_path, Center, t);

	if (typeflag == CLK || typeflag == EPHCLK)
		DownloadClk(des_path, Center, t);

	return 1;
}

/// DCB data download
int Download_DCB(char *des_path, gnsstime &t, unsigned int typeflag)
{
	char dcbf_name[20] = { 0 };
	char download_command[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char decompression_command[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };

	if (!FileExist(des_path))
		CreatePath(des_path);

	if ((typeflag == P1C1DCB) || (typeflag == P1C1_P1P2) || (typeflag == DCBALL))
	{
		sprintf(dcbf_name, "P1C1%02d%02d.DCB.Z", t.m_Year % 100, t.m_Month);

		sprintf(download_command, ".\\bin\\wget.exe -nv -P%s %s/CODE/%4d/%s -t 2",
			des_path, unibe_code_ftp, t.m_Year, dcbf_name);

		sprintf(decompression_command, ".\\bin\\gzip.exe -N -d %s%s", des_path, dcbf_name);

		system(download_command);
		system(decompression_command);
	}
	if ((typeflag == P1P2DCB) || (typeflag == P1C1_P1P2) || (typeflag == DCBALL))
	{
		memset(dcbf_name, 0, sizeof(char)* 20);
		memset(download_command, 0, sizeof(char)* 100);
		memset(decompression_command, 0, sizeof(char)* 100);

		sprintf(dcbf_name, "P1P2%02d%02d.DCB.Z", t.m_Year % 100, t.m_Month);

		sprintf(download_command, ".\\bin\\wget.exe -nv -P%s %s/CODE/%4d/%s -t 2",
			des_path, unibe_code_ftp, t.m_Year, dcbf_name);

		sprintf(decompression_command, ".\\bin\\gzip.exe -N -d %s%s", des_path, dcbf_name);

		system(download_command);
		system(decompression_command);
	}
	if ((typeflag == P2C2DCB) || (typeflag == DCBALL))
	{
		memset(dcbf_name, 0, sizeof(char)* 20);
		memset(download_command, 0, sizeof(char)* 100);
		memset(decompression_command, 0, sizeof(char)* 100);

		sprintf(dcbf_name, "P2C2%02d%02d_RINEX.DCB.Z", t.m_Year % 100, t.m_Month);

		sprintf(download_command, ".\\bin\\wget.exe -nv -P%s %s/CODE/%4d/%s -t 2",
			des_path, unibe_code_ftp, t.m_Year, dcbf_name);

		sprintf(decompression_command, ".\\bin\\gzip.exe -N -d %s%s", des_path, dcbf_name);

		system(download_command);
		system(decompression_command);
	}

	return 1;
}

/// Erp data download
int Download_Erp(char *des_path, gnsstime &t)
{
	char erpf_name[20] = { 0 };
	char download_command[100] = { 0 };
	char decompression_command[100] = { 0 };

	if (!FileExist(des_path))
		CreatePath(des_path);

	sprintf(erpf_name, "igs%04d7.erp.Z", t.m_Week);

	sprintf(download_command, ".\\bin\\wget.exe -nv -P%s %s%4d/%s -t 2",
		des_path, product_cddis_nasa_ftp, t.m_Week, erpf_name);

	sprintf(decompression_command, ".\\bin\\gzip.exe -N -d %s%s", des_path, erpf_name);

	system(download_command);
	system(decompression_command);

	return 1;
}

/// Download SNX file
int Download_SNX(char *des_path, char *Center, gnsstime &t)
{
	char snxf_name[20] = { 0 };
	char download_command_snx[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char decompression_command_snx[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };

	if (!FileExist(des_path)) CreatePath(des_path);

	sprintf(snxf_name, "%s%02dP%4d.snx", Center, t.m_Year % 100, t.m_Week);

	char snxpath[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	sprintf(snxpath, "%s%s", des_path, snxf_name);
	if (FileExist(snxpath)) return 2;

	sprintf(download_command_snx, ".\\bin\\wget.exe -nv -P%s %s%4d/%s.Z -t 2",
		des_path, product_cddis_nasa_ftp, t.m_Week, snxf_name);

	sprintf(decompression_command_snx, ".\\bin\\gzip.exe -N -d %s%s.Z", des_path, snxf_name);

	system(download_command_snx);
	system(decompression_command_snx);

	return 1;
}

/// Download STA file
int Download_STA(char *des_path)
{
	char staf_name[20] = "CODE.STA";
	char download_command_sta[100] = { 0 };

	if (!FileExist(des_path))
		CreatePath(des_path);

	sprintf(download_command_sta, ".\\bin\\wget.exe -nv -P%s %sCODE.STA -t 2", des_path, crd_unibe_ftp);

	system(download_command_sta);

	return 1;
}

int Download_GIM(char *des_path, char *Center, gnsstime &t)
{
	char ionf_name[50] = { 0 };
	char download_command[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char decompression_command[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };

	if (!FileExist(des_path))
		CreatePath(des_path);

	sprintf(ionf_name, "%sg%03d0.%02di.Z", Center, t.m_Doy, t.m_Year % 100);

	sprintf(download_command, ".\\bin\\wget.exe -nv -P%s %s/gps/products/ionex/%4d/%03d/%s -t 2",
		des_path, cddis_nasa_ftp, t.m_Year, t.m_Doy, ionf_name);

	sprintf(decompression_command, ".\\bin\\gzip.exe -N -d %s%s", des_path, ionf_name);

	system(download_command);
	system(decompression_command);

	return 1;
}

/// Download trop zenith path delay product
int Download_ZPD(char *des_path, char *StaName, gnsstime &t)
{
	char trpf_name[50] = { 0 };
	char trpf_name2[50] = { 0 };
	char download_command[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };
	char decompression_command[MAX_DOWNLOAD_COMMAND_LENGTH] = { 0 };

	if (!FileExist(des_path))
		CreatePath(des_path);

	sprintf(trpf_name, "%s%03d0.%02dzpd.gz", StaName, t.m_Doy, t.m_Year % 100);

	sprintf(trpf_name2, "%s%s%03d0.%02dzpd", des_path, StaName, t.m_Doy, t.m_Year % 100);
	if (FileExist(trpf_name2))
		return 1;

	sprintf(download_command, ".\\bin\\wget.exe -nv -P%s %stroposphere/zpd/%4d/%03d/%s -t 2",
		des_path, product_cddis_nasa_ftp, t.m_Year, t.m_Doy, trpf_name);

	sprintf(decompression_command, ".\\bin\\gzip.exe -N -d %s%s", des_path, trpf_name);

	system(download_command);
	system(decompression_command);

	return 1;
}
