#include "SUPREME_Meteorological.h"

#ifndef MAX_PATH_LEN
#define MAX_PATH_LEN     400        // Max line length
#endif

/// Read GPT2w Grid Data
int Trop_GPT2w::ReadGPT2wGrdFile(char *grdFileName)
{
	char strline[MAX_PATH_LEN] = { 0 };
	int strlen = MAX_PATH_LEN;
	FILE *fp = fopen(grdFileName, "r");
	if (!fp)
	{
		printf("Warrning: Trop grid file can not be open...\n");
		return 0;
	}

	if (fp)
	{
		GetFileLine(fp, strline, strlen);
		while (!feof(fp))
		{
			if (GetFileLine(fp, strline, strlen) == 0) continue;

			GPT2w_GrdVar GrdvarT;
			GrdvarT.lat = str2num(strline, 0, 6);
			GrdvarT.lon = str2num(strline, 7, 6);
			GrdvarT.pgrid[0] = str2num(strline, 14, 6);
			GrdvarT.pgrid[1] = str2num(strline, 21, 5);
			GrdvarT.pgrid[2] = str2num(strline, 27, 4);
			GrdvarT.pgrid[3] = str2num(strline, 32, 4);
			GrdvarT.pgrid[4] = str2num(strline, 37, 4);

			GrdvarT.Tgrid[0] = str2num(strline, 42, 5);
			GrdvarT.Tgrid[1] = str2num(strline, 48, 5);
			GrdvarT.Tgrid[2] = str2num(strline, 54, 4);
			GrdvarT.Tgrid[3] = str2num(strline, 59, 4);
			GrdvarT.Tgrid[4] = str2num(strline, 64, 4);

			GrdvarT.Qgrid[0] = str2num(strline, 69, 5) / 1000;
			GrdvarT.Qgrid[1] = str2num(strline, 75, 5) / 1000;
			GrdvarT.Qgrid[2] = str2num(strline, 81, 5) / 1000;
			GrdvarT.Qgrid[3] = str2num(strline, 87, 5) / 1000;
			GrdvarT.Qgrid[4] = str2num(strline, 93, 5) / 1000;

			GrdvarT.dTgrid[0] = str2num(strline, 99, 5) / 1000;
			GrdvarT.dTgrid[1] = str2num(strline, 105, 5) / 1000;
			GrdvarT.dTgrid[2] = str2num(strline, 111, 4) / 1000;
			GrdvarT.dTgrid[3] = str2num(strline, 116, 4) / 1000;
			GrdvarT.dTgrid[4] = str2num(strline, 121, 4) / 1000;

			GrdvarT.u = str2num(strline, 126, 7);
			GrdvarT.Hs = str2num(strline, 134, 8);

			GrdvarT.ahgrid[0] = str2num(strline, 143, 6) / 1000;
			GrdvarT.ahgrid[1] = str2num(strline, 150, 7) / 1000;
			GrdvarT.ahgrid[2] = str2num(strline, 158, 7) / 1000;
			GrdvarT.ahgrid[3] = str2num(strline, 166, 7) / 1000;
			GrdvarT.ahgrid[4] = str2num(strline, 174, 7) / 1000;

			GrdvarT.awgrid[0] = str2num(strline, 182, 7) / 1000;
			GrdvarT.awgrid[1] = str2num(strline, 190, 7) / 1000;
			GrdvarT.awgrid[2] = str2num(strline, 198, 7) / 1000;
			GrdvarT.awgrid[3] = str2num(strline, 206, 7) / 1000;
			GrdvarT.awgrid[4] = str2num(strline, 214, 7) / 1000;

			m_GPT2w_Grd.push_back(GrdvarT);
		}

		_isReadGridData = true;
	}

	if (fp) fclose(fp);

	return 1;
}

/// Get Meteorological Parameter by GPT2w Model
/* Parameter:
* dmjd:  Modified Julian date (scalar, only one epoch per call is possible)
* dlat:  Ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
* dlon:  Longitude in radians [-pi:pi] or [0:2pi] (vector)
* hell:  Ellipsoidal height in m (vector)
* it:    case 1: no time variation but static quantities
case 0: with time variation (annual and semiannual terms)
* ------------------------------------------------------------------------------ */
GPT2_OutputPara Trop_GPT2w::GetGPT2wMeteorologicalPara(double dmjd, double dlat, double dlon, double hell)
{
	GPT2_OutputPara m_result;

	//% change the reference epoch to January 1 2000
	double dmjd1 = dmjd - 51544.5;
	//% mean gravity in m/s**2
	double gm = 9.80665;
	//% molar mass of dry air in kg/mol
	double dMtr = 28.965e-3;
	//% universal gas constant in J/K/mol
	double Rg = 8.3143;
	double cosfy = 0, coshy = 0, sinfy = 0, sinhy = 0;

	//if (it == 1) //% then  constant parameters
	//{
	//	cosfy = 0;
	//	coshy = 0;
	//	sinfy = 0;
	//	sinhy = 0;
	//}
	//else
	{
		cosfy = cos(dmjd1 / 365.25 * 2 * PI);
		coshy = cos(dmjd1 / 365.25 * 4 * PI);
		sinfy = sin(dmjd1 / 365.25 * 2 * PI);
		sinhy = sin(dmjd1 / 365.25 * 4 * PI);
	}

	double plon = 0, ppod = 0;
	//% only positive longitude in degrees
	if (dlon < 0)
		plon = (dlon + 2 * PI) * 180 / PI;
	else
		plon = dlon * 180 / PI;
	// % transform to polar distance in degrees
	ppod = (-dlat + PI / 2) * 180 / PI;
	//% find the index (line in the grid file) of the nearest point
	double ipod = floor((ppod + 5) / 5);
	double ilon = floor((plon + 5) / 5);
	//% normalized (to one) differences, can be positive or negative
	double diffpod = (ppod - (ipod * 5 - 2.5)) / 5;
	double difflon = (plon - (ilon * 5 - 2.5)) / 5;
	// % added by HCY
	if (ipod == 37)
		ipod = 36;
	//% get the number of the corresponding line
	int indx1 = (ipod - 1) * 72 + ilon - 1;
	//% near the poles: nearest neighbour interpolation, otherwise: bilinear
	double bilinear = 0;
	if (ppod > 2.5 && ppod < 177.5)
		bilinear = 1;

	//% case of nearest neighbourhood
	if (bilinear == 0)
	{
		int ix = indx1;
		GPT2w_GrdVar tempGrdVar = m_GPT2w_Grd[ix];
		// % transforming ellipsoidial height to orthometric height
		m_result.undu = tempGrdVar.u;
		double hgt = hell - tempGrdVar.u;
		// % pressure, temperature at the heigtht of the grid
		double T0 = tempGrdVar.Tgrid[0] +
			tempGrdVar.Tgrid[1] * cosfy + tempGrdVar.Tgrid[2] * sinfy +
			tempGrdVar.Tgrid[3] * coshy + tempGrdVar.Tgrid[4] * sinhy;
		double p0 = tempGrdVar.pgrid[0] +
			tempGrdVar.pgrid[1] * cosfy + tempGrdVar.pgrid[2] * sinfy +
			tempGrdVar.pgrid[3] * coshy + tempGrdVar.pgrid[4] * sinhy;
		//% specific humidity
		double Q = tempGrdVar.Qgrid[0] +
			tempGrdVar.Qgrid[1] * cosfy + tempGrdVar.Qgrid[2] * sinfy +
			tempGrdVar.Qgrid[3] * coshy + tempGrdVar.Qgrid[4] * sinhy;

		//% lapse rate of the temperature
		m_result.dT = tempGrdVar.dTgrid[0] +
			tempGrdVar.dTgrid[1] * cosfy + tempGrdVar.dTgrid[2] * sinfy +
			tempGrdVar.dTgrid[3] * coshy + tempGrdVar.dTgrid[4] * sinhy;
		//% station height - grid height
		double redh = hgt - tempGrdVar.Hs;
		//  % temperature at station height in Celsius
		m_result.T = T0 + m_result.dT*redh - 273.15;
		// % temperature lapse rate in degrees / km
		m_result.dT = m_result.dT * 1000;
		// % virtual temperature in Kelvin
		double Tv = T0*(1 + 0.6077*Q);

		double c = gm*dMtr / (Rg*Tv);

		//% pressure in hPa
		m_result.p = (p0*exp(-c*redh)) / 100;

		//% water vapour pressure in hPa
		m_result.e = (Q*m_result.p) / (0.622 + 0.378*Q);

		//% hydrostatic coefficient ah 
		m_result.ah = tempGrdVar.ahgrid[0] +
			tempGrdVar.ahgrid[1] * cosfy + tempGrdVar.ahgrid[2] * sinfy +
			tempGrdVar.ahgrid[3] * coshy + tempGrdVar.ahgrid[4] * sinhy;

		//% wet coefficient aw
		m_result.aw = tempGrdVar.awgrid[0] +
			tempGrdVar.awgrid[1] * cosfy + tempGrdVar.awgrid[2] * sinfy +
			tempGrdVar.awgrid[3] * coshy + tempGrdVar.awgrid[4] * sinhy;

	}
	else if (bilinear == 1)
	{// % bilinear interpolation
		double signpod = 0;
		double signlon = 0;
		if (diffpod > 0)
			signpod = 1;
		else if (diffpod < 0)
			signpod = -1;
		if (difflon > 0)
			signlon = 1;
		else if (difflon < 0)
			signlon = -1;

		double ipod1 = ipod + signpod;// sign(diffpod);
		double ilon1 = ilon + signlon;// sign(difflon);
		if (ilon1 == 73)
			ilon1 = 1;
		else if (ilon1 == 0)
			ilon1 = 72;
		// % get the number of the line
		double indx2 = (ipod1 - 1) * 72 + ilon - 1;//  % along same longitude
		double indx3 = (ipod - 1) * 72 + ilon1 - 1;//% along same polar distance
		double indx4 = (ipod1 - 1) * 72 + ilon1 - 1;// % diagonal
		//定义变量
		double undul[4] = { 0 }, Ql[4] = { 0 }, dTl[4] = { 0 }, Tl[4] = { 0 },
			pl[4] = { 0 }, ahl[4] = { 0 }, awl[4] = { 0 };
		double Indexflag[4] = { indx1, indx2, indx3, indx4 };
		for (int l = 0; l < 4; l++)
		{
			GPT2w_GrdVar tempGrdVar = m_GPT2w_Grd[Indexflag[l]];
			// % transforming ellipsoidial height to orthometric height:
			//% Hortho = -N + Hell
			undul[l] = tempGrdVar.u;
			double hgt = hell - undul[l];

			//% pressure, temperature at the heigtht of the grid
			double T0 = tempGrdVar.Tgrid[0] +
				tempGrdVar.Tgrid[1] * cosfy + tempGrdVar.Tgrid[2] * sinfy +
				tempGrdVar.Tgrid[3] * coshy + tempGrdVar.Tgrid[4] * sinhy;
			double p0 = tempGrdVar.pgrid[0] +
				tempGrdVar.pgrid[1] * cosfy + tempGrdVar.pgrid[2] * sinfy +
				tempGrdVar.pgrid[3] * coshy + tempGrdVar.pgrid[4] * sinhy;

			//% humidity 
			Ql[l] = tempGrdVar.Qgrid[0] +
				tempGrdVar.Qgrid[1] * cosfy + tempGrdVar.Qgrid[2] * sinfy +
				tempGrdVar.Qgrid[3] * coshy + tempGrdVar.Qgrid[4] * sinhy;

			// % reduction = stationheight - gridheight
			double Hs1 = tempGrdVar.Hs;
			double redh = hgt - Hs1;

			//% lapse rate of the temperature in degree / m
			dTl[l] = tempGrdVar.dTgrid[0] +
				tempGrdVar.dTgrid[1] * cosfy + tempGrdVar.dTgrid[2] * sinfy +
				tempGrdVar.dTgrid[3] * coshy + tempGrdVar.dTgrid[4] * sinhy;

			// % temperature reduction to station height
			Tl[l] = T0 + dTl[l] * redh - 273.15;

			// % virtual temperature
			double Tv = T0*(1 + 0.6077*Ql[l]);
			double c = gm*dMtr / (Rg*Tv);

			//% pressure in hPa
			pl[l] = (p0*exp(-c*redh)) / 100;

			// % hydrostatic coefficient ah
			ahl[l] = tempGrdVar.ahgrid[0] +
				tempGrdVar.ahgrid[1] * cosfy + tempGrdVar.ahgrid[2] * sinfy +
				tempGrdVar.ahgrid[3] * coshy + tempGrdVar.ahgrid[4] * sinhy;

			//% wet coefficient aw
			awl[l] = tempGrdVar.awgrid[0] +
				tempGrdVar.awgrid[1] * cosfy + tempGrdVar.awgrid[2] * sinfy +
				tempGrdVar.awgrid[3] * coshy + tempGrdVar.awgrid[4] * sinhy;
		}

		double dnpod1 = abs(diffpod); //% distance nearer point
		double dnpod2 = 1 - dnpod1; //  % distance to distant point
		double dnlon1 = abs(difflon);
		double dnlon2 = 1 - dnlon1;

		// % pressure
		double R1 = dnpod2*pl[0] + dnpod1*pl[1];
		double R2 = dnpod2*pl[2] + dnpod1*pl[3];
		m_result.p = dnlon2*R1 + dnlon1*R2;

		//  % temperature
		R1 = dnpod2*Tl[0] + dnpod1*Tl[1];
		R2 = dnpod2*Tl[2] + dnpod1*Tl[3];
		m_result.T = dnlon2*R1 + dnlon1*R2;

		// % temperature in degree per km
		R1 = dnpod2*dTl[0] + dnpod1*dTl[1];
		R2 = dnpod2*dTl[2] + dnpod1*dTl[3];
		m_result.dT = (dnlon2*R1 + dnlon1*R2) * 1000;

		// % humidity
		R1 = dnpod2*Ql[0] + dnpod1*Ql[1];
		R2 = dnpod2*Ql[2] + dnpod1*Ql[3];
		double Q = dnlon2*R1 + dnlon1*R2;
		m_result.e = (Q*m_result.p) / (0.622 + 0.378*Q);

		//% hydrostatic
		R1 = dnpod2*ahl[0] + dnpod1*ahl[1];
		R2 = dnpod2*ahl[2] + dnpod1*ahl[3];
		m_result.ah = dnlon2*R1 + dnlon1*R2;

		//% wet
		R1 = dnpod2*awl[0] + dnpod1*awl[1];
		R2 = dnpod2*awl[2] + dnpod1*awl[3];
		m_result.aw = dnlon2*R1 + dnlon1*R2;

		//% undulation
		R1 = dnpod2*undul[0] + dnpod1*undul[1];
		R2 = dnpod2*undul[2] + dnpod1*undul[3];
		m_result.undu = dnlon2*R1 + dnlon1*R2;
	}

	return m_result;
}

/// Get GNSS-based PWV
double GNSS_PWV::GetGNSSPWV()
{
	const double k1 = 77.604;
	const double k2 = 64.79;
	const double k3 = 3.776e5;       //  K ^ 2 / hPa
	const double k2_p = 17;          //  K / hPa
	const double rou = 1000;         //  kg / m ^ 3
	const double Rv = 461.51;        //  J / K / kg

	m_PI = (1e6 / (rou*Rv*(k3 / m_Tm + k2_p))) * 100;
	m_ZWD = m_ZTD - m_ZHD;
	m_PWV = m_PI*m_ZWD;

	return m_PWV;
}

