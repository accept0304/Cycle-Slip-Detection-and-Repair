/* -------------------------------------------------------------------------
* SUPREME_GeneralData.cpp : GNSS General Data
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2019.07.30
* ------------------------------------------------------------------------- */
#ifndef SUPREME_GENERAL_DATA_H_HH
#define SUPREME_GENERAL_DATA_H_HH

#include "SUPREME_Options.h"
#include "SUPREME_Anntenna.h"
#include "SUPREME_DCB.h"
#include "SUPREME_Tide.h"
#include "SUPREME_Ionosphere.h"
#include "SUPREME_Ifcb.h"
/// General data class
class GenData
{
public:
	GenData();
	~GenData();

	double OceanDisp[6 * 11];               // Ocean tide loading parameters

	erp_t m_Erp;                            // Erp data

	pcv_t *m_pcvs;                          // Satellite PCV
	pcv_t m_pcvr;                           // Receiver PCV
	SatelliteAntenna AntSat;                // Satellite Antenna Infomation
	ReceiverAntenna AntRec;                 // Receiver Antenna Infomation

	CodeDCB m_Dcb;                          // Gnss DCB Data

	GIM_t GIMdata;                          // GIM ionosphere data
	//IonexModel IonData;
	IFCB  IFCL_bias;

	//int GenDataCheck(ppp_option_t &popt);        // Check general file
	int LoadGenData(ppp_option_t &popt);         // Load general data

private:
	void setpcv(ppp_option_t &popt, const pcvs_t *pcvs, const pcvs_t *pcvr);
};

#endif


