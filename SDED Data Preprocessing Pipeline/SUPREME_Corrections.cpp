#include "SUPREME_Corrections.h"
#include "SUPREME_Constant.h"
#include "SUPREME_CommonFunction.h"
#include "SUPREME_Coordinate.h"
#include "SUPREME_SunMoon.h"

/// Pseudo range corrected by receiver and satellite pcv
double PseudoRange_PCV_Correct(double pcv_rec, double pcv_sat, double pr)
{
	if (pr == 0)return 0;
	else return pr - pcv_rec - pcv_sat;
}

/// Phase corrected by receiver and satellite pcv (Unit:cycle)
double CarrierPhase_PCV_Correct(double wavelength, double pcv_rec, double pcv_sat, double L)
{
	if (L == 0)return 0;
	else return (wavelength > 0) ? (L - (pcv_rec + pcv_sat) / wavelength) : L;
}

/// Relative effection
double RelativeEffect(double SatPos[3], double SatVel[3])
{
	return (-2.0 * (SatPos[0] * SatVel[0] + SatPos[1] * SatVel[1] + SatPos[2] * SatVel[2]) / LIGHTSPEED);
}

/// Sagnac effection
double SagnacEffect(double satpos[3], double sta[3])
{
	return OMEGA_GPS*(satpos[0] * sta[1] - satpos[1] * sta[0]) / LIGHTSPEED;
}

/// Gravitation effect
double GravitationEffect(double ss, double rr, double rs)
{
	return (2 * E_GM / LIGHTSPEED / LIGHTSPEED*log((ss + rr + rs) / (ss + rr - rs)));
}

/// BDS satellite-induced code pseudorange variations correct
double BDSMultipathCorr(unsigned int satno, double elev, double dmp[3])
{
	const static double IGSOCOEF[3][10] = {		/* m */
		{ -0.55, -0.40, -0.34, -0.23, -0.15, -0.04, 0.09, 0.19, 0.27, 0.35 },	//B1
		{ -0.71, -0.36, -0.33, -0.19, -0.14, -0.03, 0.08, 0.17, 0.24, 0.33 },	//B2
		{ -0.27, -0.23, -0.21, -0.15, -0.11, -0.04, 0.05, 0.14, 0.19, 0.32 },	//B3
	};
	const static double MEOCOEF[3][10] = {		/* m */
		{ -0.47, -0.38, -0.32, -0.23, -0.11, 0.06, 0.34, 0.69, 0.97, 1.05 },	//B1
		{ -0.40, -0.31, -0.26, -0.18, -0.06, 0.09, 0.28, 0.48, 0.64, 0.69 },	//B2
		{ -0.22, -0.15, -0.13, -0.10, -0.04, 0.05, 0.14, 0.27, 0.36, 0.47 },	//B3
	};
	const static double IGSOCOEF_LIXIN[3][17]
	{
		{-0.238, -0.54, -0.43 ,-0.274 ,-0.255, -0.265, -0.168, -0.137, -0.064, -0.019,	0.025,	0.13,	0.175,	0.238,	0.234,	0.272,	0.302 },
		{-0.246, -0.421, -0.296 ,-0.272,-0.26, -0.227, -0.154, -0.115, -0.068, -0.041,	0.037,	0.079,	0.116,	0.167,	0.202	,0.25,	0.259},
		{-0.406, -0.211 ,-0.169 ,-0.157, -0.138 ,-0.253, -0.089, -0.118, -0.073, -0.01 ,-0.008,	0.049,	0.084,	0.127,	0.152,	0.207,0.155}
	};
	const static double MEOCOEF_LIXIN[3][17]
	{
		{-0.375, -0.405 ,-0.314, -0.22, -0.18, -0.144 ,-0.179, -0.091,	0.005,	0.081,	0.222,	0.322,	0.463,	0.631,	0.716,	0.918,	0.955},
		{-0.4, -0.238, -0.213, -0.228, -0.131, -0.113, -0.114, -0.065,	0.006,	0.068,	0.187,	0.224,	0.326,	0.44,	0.485,	0.583,	0.628 },
		{-0.134, -0.174, -0.112 ,-0.131 ,-0.047, -0.052, -0.081, -0.03,	0.044,	0.048,	0.083,	0.133,	0.201,	0.287,	0.288,	0.367,	0.393}

	};


	char sys = 0;
	unsigned int prn = 0;
	int i = 0, j = 0, b = 0;
	double a = 0.0;

	sys = GetSysPrn(satno, prn);

	if (sys != 'C') return 0;

	if (elev <= 0.0) return 0;

	a = elev * 180 / PI;
	b = (int)(a / 5);
	a /= 5;
	if (prn >= 6 && prn < 11) { // IGSO(C06, C07, C08, C09, C10)
		if (b <= 0) {
			for (j = 0; j < 3; j++) dmp[j] = IGSOCOEF_LIXIN[j][0];
		}
		else if (b >= 16) {
			for (j = 0; j < 3; j++) dmp[j] = IGSOCOEF_LIXIN[j][16];
		}
		else {
			for (j = 0; j < 3; j++) dmp[j] = IGSOCOEF_LIXIN[j][b] * (1.0 - a + b) + IGSOCOEF_LIXIN[j][b + 1] * (a - b);
		}
	}
	else if (prn >= 11 && prn < 15) {   // MEO(C11, C12, C13, C14)
		if (b <= 0) {
			for (j = 0; j < 3; j++) dmp[j] = MEOCOEF_LIXIN[j][0];
		}
		else if (b >= 16) {
			for (j = 0; j < 3; j++) dmp[j] = MEOCOEF_LIXIN[j][16];
		}
		else {
			for (j = 0; j < 3; j++) dmp[j] = MEOCOEF_LIXIN[j][b] * (
				-a + b) + MEOCOEF[j][b + 1] * (a - b);
		}
	}

	return 1;
}