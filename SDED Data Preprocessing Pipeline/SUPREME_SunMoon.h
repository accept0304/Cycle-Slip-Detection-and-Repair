/* -------------------------------------------------------------------------
* SUPREME_SunMoon.cpp : Sun and Moon Position
*
*           Copyright (C) by C. Zhao
*           Contact caszcb@163.com
*
* Create : 2018.03.16
* ------------------------------------------------------------------------- */

#ifndef SUPREME_SUN_MOON_H_HH
#define SUPREME_SUN_MOON_H_HH

#include "SUPREME_GnssTime.h"

#define SUN_POSITION   1              // sun position flag
#define MOON_POSITION  2              // moom position flag

#define PLANET_MAX_EPOCH       1000
#define SUN_MOON_MAX_EPOCH     1000

/// Planet Position
class SunMoonPos
{
public:
	SunMoonPos()
	{
		m_Mjd = 0;
		m_SunPos[0] = 0; m_SunPos[1] = 0; m_SunPos[2] = 0;
		m_MoonPos[0] = 0; m_MoonPos[1] = 0; m_MoonPos[2] = 0;
	}

	double m_Mjd;
	double m_SunPos[3];
	double m_MoonPos[3];
};

class SunMoonEph
{
public:
	SunMoonEph(){};

	int count;
	double m_Interval;
	SunMoonPos Eph[SUN_MOON_MAX_EPOCH];
};

void GetSunMoonPos(gnsstime &obst, double *rsun, double *rmoon, double *gmst);

#endif
