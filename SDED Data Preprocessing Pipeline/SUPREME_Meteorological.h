/* -------------------------------------------------------------------------
* SUPREME_Meteorological.h : SUPREME for Meteorological Applications
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2019.01.08
* ------------------------------------------------------------------------- */
#include "SUPREME_Troposphere.h"
#include "SUPREME_CommonFunction.h"

typedef struct GPT2_OutputPara
{
	double p;               // Pressure in hPa
	double T;               // Temperature in degrees
	double dT;              // Temperature lapse rate in degrees per km
	double e;               // Water vapour pressure in hPa
	double ah;              // Hydrostatic mapping function coefficient at zero height (VMF1)
	double aw;              // Wet mapping function coefficient (VMF1)
	double undu;            // Geoid undulation in meter
};

typedef struct GPT2w_GrdVar
{
	double lat;             // Latitude
	double lon;             // Lontitude
	double pgrid[5];        // Pressure in Pascal
	double Tgrid[5];        // Temperature in Kelvin
	double Qgrid[5];        // Specific humidity in kg/kg
	double dTgrid[5];       // Temperature lapse rate in Kelvin/m
	double u;               // Geoid undulation in m
	double Hs;              // Orthometric grid height in m
	double ahgrid[5];       // Hydrostatic mapping function coefficient, dimensionless
	double awgrid[5];       // Wet mapping function coefficient, dimensionless
};

class Trop_GPT2w
{
public:
	Trop_GPT2w();
	~Trop_GPT2w();

	vector<GPT2w_GrdVar> m_GPT2w_Grd;
	int ReadGPT2wGrdFile(char *grdFileName);
	GPT2_OutputPara GetGPT2wMeteorologicalPara(double dmjd, double dlat, double dlon, double hell);

private:
	bool _isReadGridData;
};

class GNSS_PWV
{
	GNSS_PWV();
	~GNSS_PWV();

	double m_Ts, m_Tm, m_Ps, m_PI;
	double m_lat, m_lon, m_h;               // GNSS station position (geodetic coordinate)
	double m_ZTD, m_ZHD, m_ZWD, m_PWV;
	gnsstime m_t;

	double GetGNSSPWV();

	unsigned int _TmModel;
private:
	double rou, Rv, k3, k2;
};
