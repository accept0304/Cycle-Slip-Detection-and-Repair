/* -------------------------------------------------------------------------
* SUPREME_Corrections.cpp : GNSS Corrections
*
*           Copyright (C) by C. Zhao  All rights reserved
*           Contact caszcb@163.com
*
* Create : 2018.03.16
* ------------------------------------------------------------------------- */
#ifndef SUPREME_GNSS_CORRECTIONS_H_HH
#define SUPREME_GNSS_CORRECTIONS_H_HH

double PseudoRange_PCV_Correct(double pcv_rec, double pcv_sat, double Pr);
double CarrierPhase_PCV_Correct(double wavelength, double pcv_rec, double pcv_sat, double L);

double RelativeEffect(double SatPos[3], double SatVel[3]);

double SagnacEffect(double satpos[3], double sta[3]);

double GravitationEffect(double ss, double rr, double rs);

double BDSMultipathCorr(unsigned int satno, double elev, double dmp[3]);

#endif
