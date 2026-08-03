#ifndef SUPREME_QUALITY_CHECK_H_HH
#define SUPREME_QUALITY_CHECK_H_HH

#include "SUPREME_Options.h"
#include "SUPREME_DataProcess.h"
#include "SUPREME_Matrix.h"
#include "SUPREME_OrbClk.h"

void CycleSlipInit();

int RecClkSlipRepair(ppp_option_t &popt, unsigned int Satn, unsigned int *Prn, double *Pr1, double *Pr2, double *L1, double *L2);

int RC_Repair(ppp_option_t& popt, PPPObsData& inputdata, NavData* navdata, double* Ele, unsigned int flag_epoch_count);

int QC_CycleSlip_ClockSlip(ppp_option_t &popt, PPPObsData &inputdata, NavData *navdata);


int PseudoRange_C1P1_Detect(const double &c1, const double &p1);

int QC_L(MatrixT L);



#endif