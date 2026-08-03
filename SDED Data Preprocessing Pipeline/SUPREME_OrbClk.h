/* -------------------------------------------------------------------------
* SUPREME_OribtClock.cpp : GNSS Orbit and Clock
*
*           Copyright (C) by C. Zhao
*           Contact caszcb@163.com
*
* Create : 2019.07.29
* ------------------------------------------------------------------------- */

#ifndef SUPREME_ORBITCLOCK_H_HH
#define SUPREME_ORBITCLOCK_H_HH

#include "SUPREME_GnssTime.h"
#include "SUPREME_NavRNX.h"
#include "SUPREME_ClkRNX.h"
#include "SUPREME_Anntenna.h"

/// Satellite information class -----------------------------------------------------
// Function: Temporary storage of satellite Postioin,Velocity,Clock,Element,Azimuth
// Parameter:
//            flag   :   0:only navigation result   1:precise result
//            m_Size :   usefule satellite count
//            m_Prn  :   Vehicle prn
//            m_Xs   :   Satellite postion
//            m_Vs   :   Satellite velocity
//            m_Clk  :   Satellite clock error
//            m_Ele  :   Satellite Elevation
//            m_Azi  :   Satellite Azimuth
//
// Add in 2017-9-2
// ----------------------------------------------------------------------------------
class SatInfo
{
public:
	SatInfo();
	SatInfo(unsigned int satn);
	~SatInfo();

	unsigned int m_Size;
	unsigned int *m_Prn, *flag;

	double *m_Xs, *m_Vs, *m_Clk, *m_Ele, *m_Azi;
	int DeleSat(int index);
	int Initialize(unsigned int satn);
	int Clear();
private:
	unsigned int release;
};

int GetBdsIODE_by_CRC24(BdsEph *ptr);
int GetBdsIODE_by_TOE(BdsEph *ptr);
int GetGloIODE(GloEph *ptr);

/// Get Navigation Ephemeris by Prn and IODE
GpsEph* Find_GpsNav_by_Iode(unsigned int prn, int iode, NavData* navdata);
GloEph* Find_GloNav_by_Iode(unsigned int prn, int iode, NavData* navdata);
GalEph* Find_GalNav_by_Iode(unsigned int prn, int iode, NavData* navdata);
BdsEph* Find_BdsNav_by_Iode(unsigned int prn, int iode, NavData* navdata, int flag);
void* Find_NavEph_by_Iode(unsigned int satno, int iode, NavData* navdata, int flag);

GpsEph* Find_GPS_NavEph(unsigned int prn, gnsstime& t, NavData* navdata);
GloEph* Find_GLO_NavEph(unsigned int prn, gnsstime& t, NavData* navdata);
GalEph* Find_GAL_NavEph(unsigned int prn, gnsstime& t, NavData* navdata);
BdsEph* Find_BDS_NavEph(unsigned int prn, gnsstime& t, NavData* navdata);
void* Find_NavEph(unsigned int prn, gnsstime& t, NavData* navdata);
void* Find_NavEph(PreciseData &predata, unsigned int Prn, gnsstime &obs_t);

double motion_equationX(double x, double y, double z, double Vy, double Ax);
double motion_equationY(double x, double y, double z, double Vx, double Ay);
double motion_equationZ(double x, double y, double z, double Az);

/// Calculate Satellite Position and Velocity by Navigation
int Cal_GPS_PosVel(gnsstime& t, GpsEph* gpseph, double pos[3], double vel[3]);
int Cal_GLO_PosVel(gnsstime& t, GloEph* gloeph, double pos[3], double vel[3]);
int Cal_GAL_PosVel(gnsstime& t, GalEph* galeph, double pos[3], double vel[3]);
int Cal_BDS_PosVel(gnsstime& t, BdsEph* bdseph, double pos[3], double vel[3]);

/// Calculate Satellite Clock Offset by Navigation
double Cal_GPS_Clk(gnsstime& t, GpsEph* gpseph);
double Cal_GLO_Clk(gnsstime& t, GloEph* gloeph);
double Cal_GAL_Clk(gnsstime& t, GalEph* galeph);
double Cal_BDS_Clk(gnsstime& t, BdsEph* bdseph);

double Precise_Interp_Clock(gnsstime& time, Clk_EpochData *ptr, int len, int n);

/// Broadcast satellite pos and clk to precise satellite pos and clk
int Nav2Pre(double pos[3], double vel[3], double dr, double da, double dc, double dx[3]);
int Eph2PreEph(gnsstime& obst, double pos[3], double vel[3], void* ptr_eph, int sys);
int Clk2Preclk(gnsstime& obst, double& dclk, void* ptreph, int sys);

/// Get satellite pos, vel, clk realtime
double Cal_Satellite_PosVel_byTime(void* eph, gnsstime& obst, double pos[3], double vel[3],
	double& SatClk, char Sys, unsigned int flag, double& TGD1, double& TGD2);
int Cal_Sat_PosVel_Real(unsigned int prn, gnsstime& t, double pr, double Sta[3], void* eph,
	double pos[3], double vel[3], double& SatClk, int flag, double& TGD1, double& TGD2);


/// Get satelltie pos, vel, clk post
int Cal_Sat_PosVel_Post(unsigned int prn, gnsstime& t, double pr, double Sta[3],
	PreEpochData* ptr_pre, int len_pre, Clk_EpochData* ptr_clk, int len_clk,
	double satpos[3], double satvel[3], double& SatClk);

/// Get satellite elevatioin and azimuth
int Cal_Sat_EleAzimuth(double sta[3], double sat[3], double& el, double& azimuth);

int EarthShadow(double sun[3], double sat[3], char *sat_type);


typedef struct zhyinfo 
{
	unsigned int m_Prn[100];
	double m_Xs[300];
	unsigned int m_Size;
	double Ele[300];
}zhyinfo;
#endif