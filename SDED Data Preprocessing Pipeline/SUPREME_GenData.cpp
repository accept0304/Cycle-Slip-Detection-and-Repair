#include <string.h>
#include <math.h>

#include "SUPREME_GenData.h"
#include "SUPREME_CommonFunction.h"
#include "SUPREME_Coordinate.h"
#include "SUPREME_DownLoad.h"

static pcvs_t pcvss = { 0 };        /* receiver antenna parameters */
static pcvs_t pcvsr = { 0 };        /* satellite antenna parameters */

/// Set Antenna parameters
void GenData::setpcv(ppp_option_t &popt, const pcvs_t *pcvs, const pcvs_t *pcvr)
{
	char sys = 0;
	unsigned int prn = 0;
	unsigned int i = 0, j = 0, k = 0;
	double pos[3] = { 0.0 }, del[3] = { 0.0 }, dt = 0.0;
	char id[64] = { 0 };

	pcv_t *pcv = NULL;

	if (!m_pcvs) return;

	// Satellite antenna parameter
	for (i = 0; i<GNSS_SATNO_NUM; i++)
	{
		if (!(pcv = searchpcv(i + 1, "", popt.CurrentEpoch.m_gtime, pcvs)))
		{
			continue;
		}
		//m_pcvs[i] = *pcv;
		m_pcvs[pcv->sat] = *pcv;  // modified by zcb

		if (pcv->dzen == 0.0) j = 10;
		else j = round((pcv->zen2 - pcv->zen1) / pcv->dzen);
	}

	// Receiver antenna parameter
	for (i = 0; i<1; i++)
	{
		//if (!strcmp(popt.opt3.antenna_type, "*")) 
		//{ 
		//	/* set by station parameters */
		//	strcpy(popt.opt3.antenna_type, sta[i].antdes);
		//	if (sta[i].deltype == 1) 
		//	{ 
		//		/* xyz */
		//		if (Norm(popt.ProXYZ, 3)>0.0) 
		//		{
		//			ecef2pos(popt.ProXYZ, pos);
		//			ecef2enu(pos, popt.opt3.pro_ant, del);
		//			for (j = 0; j<3; j++) popt.opt3.pro_ant[j] = del[j];
		//		}
		//	}
		//	else 
		//	{ 
		//		/* enu */
		//		for (j = 0; j<3; j++) popt->antdel[j] = stas[i].del[j];
		//	}
		//}

		if (!(pcv = searchpcv(0, popt.antenna_type, popt.CurrentEpoch.m_gtime, pcvr)))
		{
			continue;
		}

		m_pcvr = *pcv;
	}
}

/// Check General data file
//int GenData::GenDataCheck(ppp_option_t &popt)
//{
//	unsigned int SysNum = 0;
//	char igsCenter[5] = { 0 };
//	gnsstime t = popt.CurrentEpoch;
//
//	if (popt.System&PRO_SYS_GPS)
//		SysNum++;
//	if (popt.System&PRO_SYS_GLO)
//		SysNum++;
//	if (popt.System&PRO_SYS_BDS)
//		SysNum++;
//	if (popt.System&PRO_SYS_GAL)
//		SysNum++;
//
//	if (SysNum > 1) // multi-gnss
//		sprintf(igsCenter, "gbm");
//	else
//		sprintf(igsCenter, "igs");
//
//	/// ANTEX
//	if (!FileExist(popt.IOfile.gen_f.atx_fl))
//	{
//		if (FileExist("gen/igs08.atx"))
//		{
//			memset(popt.IOfile.gen_f.atx_fl, 0, sizeof(char)*PATH_LENGTH);
//			strcpy(popt.IOfile.gen_f.atx_fl, "gen/igs08.atx");
//		}
//		else if (FileExist("gen/igs14.atx"))
//		{
//			memset(popt.IOfile.gen_f.atx_fl, 0, sizeof(char)*PATH_LENGTH);
//			strcpy(popt.IOfile.gen_f.atx_fl, "gen/igs14.atx");
//		}
//		else
//		{
//			printf("Error:No PCO and PCV file input.\n");
//			return 0;
//		}
//	}

	//if (popt.process_mode == POST_PPP_MODE)
	//{
	//	/// Precise product
	//	if (!FileExist(popt.IOfile.rinex_f.prephf))
	//	{
	//		char eph_path[PATH_LENGTH] = { 0 };
	//		sprintf(eph_path, "orb/%s%4d%d.sp3", igsCenter, t.m_Week, t.m_Dow);
	//		memset(popt.IOfile.rinex_f.prephf, 0, sizeof(char)*PATH_LENGTH);
	//		sprintf(popt.IOfile.rinex_f.prephf, eph_path);

	//		EphClk_download("orb/", t, igsCenter, EPH);
	//	}
	//	if (!FileExist(popt.IOfile.rinex_f.preclkf))
	//	{
	//		char clk_path[PATH_LENGTH] = { 0 };
	//		sprintf(clk_path, "orb/%s%4d%d.clk", igsCenter, t.m_Week, t.m_Dow);
	//		memset(popt.IOfile.rinex_f.preclkf, 0, sizeof(char)*PATH_LENGTH);
	//		sprintf(popt.IOfile.rinex_f.preclkf, clk_path);

	//		EphClk_download("orb/", t, igsCenter, CLK);
	//	}

		///// Erp data
		//if (!FileExist(popt.IOfile.gen_f.erp_fl))
		//{
		//	char erp_path[PATH_LENGTH] = { 0 };
		//	sprintf(erp_path, "gen/erp/%s%4d%d.erp", "igs", t.m_Week, t.m_Dow);
		//	memset(popt.IOfile.gen_f.erp_fl, 0, sizeof(char)*PATH_LENGTH);
		//	sprintf(popt.IOfile.gen_f.erp_fl, erp_path);

		//	Erp_download("gen/erp/", t);
		//}

		///// DCB data
		//if (!FileExist(popt.IOfile.gen_f.dcb_p1c1_fl))
		//{
		//	char p1c1dcb_path[PATH_LENGTH] = { 0 };
		//	sprintf(p1c1dcb_path, "gen/dcb/%s%02d%0d.DCB", "P1C1", t.m_Year % 100, t.m_Month);
		//	memset(popt.IOfile.gen_f.dcb_p1c1_fl, 0, sizeof(char)*PATH_LENGTH);
		//	sprintf(popt.IOfile.gen_f.dcb_p1c1_fl, p1c1dcb_path);

		//	Dcb_download("gen/dcb/", t, P1C1DCB);
		//}
	//}

	//return 1;
//}

/// Read GNSS General Data
int GenData::LoadGenData(ppp_option_t &popt)
{
	// Gen data file check
	//GenDataCheck(popt);

	// Get Satellite Antenna PCO and PCV
	if (ReadAntex(popt.IOfile.gen_f.atx_fl, &pcvss) == 0) return 0;
	setpcv(popt, &pcvss, &pcvss);

	// Read DCB data: P1-C1  P1-P2
	m_Dcb.Read_DCB(popt.IOfile.gen_f.dcb_p1c1_fl);
	m_Dcb.Read_DCB(popt.IOfile.gen_f.dcb_p1p2_fl);
	m_Dcb.Read_DCB(popt.IOfile.gen_f.dcb_p2c2_fl);
	//m_Dcb.Read_mgex_DCB(popt.IOfile.gen_f.dcb_mgex);
	m_Dcb.Read_mgex_DCB(popt.IOfile.gen_f.dcb_mgex);
	IFCL_bias.read_ifcb(popt.IOfile.gen_f.ifcb_f);
	// Read ERP data
	ReadErpFile(popt.IOfile.gen_f.erp_fl, &m_Erp);

	// Read BLQ data
	ReadBlqFile(popt.IOfile.gen_f.blq_fl, popt.StaName, OceanDisp);

	// Read GIM data if use ION_MODEL_GIM model
	//popt.ion_model = 3;
	//if (popt.ion_model >= ION_MODEL_GIM_CORRECT)//¶ÁµçÀë²ãGIMÎÄ¼þ£»
	//{
	ReadTecGridFile(popt.IOfile.gen_f.ion_fl, &GIMdata, 0);

		//IonData.SetGimFile(popt.IOfile.gen_f.ion_fl);
		//IonData.ReadIonexGimData();
	//}

	return 1;
}

GenData::GenData()
{
	m_pcvs = (pcv_t*)malloc(sizeof(pcv_t)*GNSS_SATNO_NUM);
	memset(m_pcvs, 0, sizeof(pcv_t)*GNSS_SATNO_NUM);
	
	//GIMdata.tec = (ionex_tec_t*)malloc(sizeof(ionex_tec_t));
		
}

GenData::~GenData()
{
	if (m_pcvs) free(m_pcvs);

	if (pcvss.pcv) { free(pcvss.pcv); pcvss.pcv = NULL; pcvss.n = pcvss.nmax = 0; }
	if (pcvsr.pcv) { free(pcvsr.pcv); pcvsr.pcv = NULL; pcvsr.n = pcvsr.nmax = 0; }

}
