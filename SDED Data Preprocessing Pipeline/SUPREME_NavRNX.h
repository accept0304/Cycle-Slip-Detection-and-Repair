/* -------------------------------------------------------------------------
* SUPREME_NavRNX.h : SUPREME for Navigation information
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2019.07.29
* ------------------------------------------------------------------------- */
#ifndef SUPREME_NAVIGATION_RINEX_H_HH
#define SUPREME_NAVIGATION_RINEX_H_HH

#include "SUPREME_GnssTime.h"

#define  SAT_NUMBER                 160

/// Gps Ephemeris of one epoch
class GpsEph
{
public:
	int    m_prn;
	double m_clock_bias;
	double m_clock_drift;
	double m_clock_driftrate;
	double m_IODE;            // 1-orbit
	double m_Crs;
	double m_Delta_n;
	double m_M0;
	double m_Cuc;             // 2-orbit
	double m_e;
	double m_Cus;
	double m_sqrt_A;
	double m_TOE;             // 3-orbit
	double m_Cic;
	double m_OMEGA0;
	double m_Cis;
	double m_i0;              // 4-orbit
	double m_Crc;
	double m_omega;
	double m_OMEGADOT;
	double m_IDOT;            // 5-orbit
	double m_L2code;
	double m_GPSweek;
	double m_L2P;
	double m_precision;       // 6-orbit
	double m_SVhealth;
	double m_TGD;
	double m_IODC;
	double m_TOW;             // 7-orbit
	double m_h;

	GpsEph* m_Next;
	gnsstime m_Toc_Gps;

	gnsstime m_RtOrbitTime;       // Real Time Orbit Time
	gnsstime m_RtClockTime;       // Real Time Clock Time

	int Nav_ana_GpsEph(FILE* fp, char* strline, unsigned int len, int flag);

	void Display_GpsEph();
	void Write_GpsEph(FILE *fp);
};

class GloEph
{
public:
	int   m_prn;
	double m_tau;
	double m_gamma;
	double m_tk;
	double m_X;                 // 1-orbit
	double m_Vx;
	double m_Ax;
	double m_health;
	double m_Y;                 // 2-orbit
	double m_Vy;
	double m_Ay;
	double m_frequency_number;
	double m_Z;                 // 3-orbit
	double m_Vz;
	double m_Az;
	double m_E;

	GloEph* m_Next;
	gnsstime m_Toc_Glo;
	gnsstime m_Toc_Gps;

	gnsstime m_RtOrbitTime;       // Real Time Orbit Time
	gnsstime m_RtClockTime;       // Real Time Clock Time

	int Nav_ana_GloEph(FILE* fp, char* strline, unsigned int len, int flag);

	void Display_GloEph();
};

class BdsEph
{
public:
	int    m_prn;
	double m_clock_bias;
	double m_clock_drift;
	double m_clock_driftrate;
	double m_IODE;         // 1-orbit
	double m_Crs;
	double m_Delta_n;
	double m_M0;
	double m_Cuc;          // 2-orbit
	double m_e;
	double m_Cus;
	double m_sqrt_A;
	double m_TOE;          // 3-orbit
	double m_Cic;
	double m_OMEGA0;
	double m_Cis;
	double m_i0;           // 4-orbit
	double m_Crc;
	double m_omega;
	double m_OMEGADOT;
	double m_IDOT;         // 5-orbit
	double m_reserved1;
	double m_bdweek;
	double m_reserved2;
	double m_precision;    // 6-orbit
	double m_svhealth;
	double m_TGD1;
	double m_TGD2;
	double m_IODC;         // 7-orbit

	BdsEph* m_Next;

	gnsstime m_Toc_Bds;
	gnsstime m_Toc_Gps;

	gnsstime m_RtOrbitTime;       // Real Time Orbit Time
	gnsstime m_RtClockTime;       // Real Time Clock Time

	int Nav_ana_BdsEph(FILE* fp, char* strline, unsigned int len, int flag);

	void Display_BdsEph();
};

class GalEph
{
public:
	int   m_prn;
	double m_clock_bias;
	double m_clock_drift;
	double m_clock_driftrate;
	double m_IODnav;          // 1-orbit
	double m_Crs;
	double m_Delta_n;
	double m_M0;
	double m_Cuc;             // 2-orbit
	double m_e;
	double m_Cus;
	double m_sqrt_A;
	double m_TOE;             // 3-orbit
	double m_Cic;
	double m_OMEGA0;
	double m_Cis;
	double m_i0;              // 4-orbit
	double m_Crc;
	double m_omega;
	double m_OMEGADOT;
	double m_IDOT;            // 5-orbit
	double m_data_source;
	double m_GALweek;         // TOE
	double m_reserved1;
	double m_SISA;            // 6-orbit
	double m_SVhealth;
	double m_BGD_E1_E5a;
	double m_BGD_E1_E5b;
	double m_TOW;             // 7-orbit
	double m_reserved2;
	double m_reserved3;
	double m_reserved4;

	GalEph* m_Next;

	gnsstime m_Toc_Gal;
	gnsstime m_Toc_Gps;

	gnsstime m_RtOrbitTime;       // Real Time Orbit Time
	gnsstime m_RtClockTime;       // Real Time Clock Time

	int Nav_ana_GalEph(FILE* fp, char* strline, unsigned int len, int flag);

	void Display_GalEph();
};

/// Navigation Header information
class NavHeader
{
public:
	NavHeader();

	char   m_Sys;               // System falg
	char   m_Type;
	int   m_Gps_Ion_flag;
	int   m_Gal_Ion_flag;
	double m_Version;           // Navigation version
	double m_Gps_ionA[4];
	double m_Gps_ionB[4];
	double m_Gal_ionA[4];

	int Read_NavigationHeader(char* filename);

	void Display_NavHeadr();
};

/// GPS navigation ephemeris all epoch
class Nav_GpsEph
{
public:
	unsigned int m_Prn;
	unsigned int m_EpochCount;

	GpsEph* m_GpsEph;
	GpsEph* m_Tail;
	Nav_GpsEph* m_Next;

};

/// GLO navigation ephemeris all epoch
class Nav_GloEph
{
public:
	unsigned int m_Prn;
	unsigned int m_EpochCount;

	GloEph* m_GloEph;
	GloEph* m_Tail;

	Nav_GloEph* m_Next;
};

/// BDS navigation ephemeris all epoch
class Nav_BdsEph
{
public:
	unsigned int m_Prn;
	unsigned int m_EpochCount;

	BdsEph* m_BdsEph;
	BdsEph* m_Tail;

	Nav_BdsEph* m_Next;
};

/// GAL navigation ephemeris all epoch
class Nav_GalEph
{
public:
	unsigned int m_Prn;
	unsigned int m_EpochCount;

	GalEph* m_GalEph;
	GalEph* m_Tail;

	Nav_GalEph* m_Next;
};

/// Navigation data
class NavData
{
public:
	NavData();
	~NavData();

	unsigned int m_Nav_Buff_Size;     // Navigation buffer size (num of epoch to save) 0:save all

	unsigned int m_GPS_PRN_Count;     // GPS Satellite number
	unsigned int m_GLO_PRN_Count;     // GLO Satellite number
	unsigned int m_BDS_PRN_Count;     // GAL Satellite number
	unsigned int m_GAL_PRN_Count;     // BDS Satellite number

	NavHeader m_Header;               // Navigation header

	Nav_GpsEph* m_Nav_GpsEph;         // GPS Navigation ephemeris
	Nav_GloEph* m_Nav_GloEph;         // GLO Navigation ephemeris
	Nav_BdsEph* m_Nav_BdsEph;         // GAL Navigation ephemeris
	Nav_GalEph* m_Nav_GalEph;         // BDS Navigation ephemeris

	Nav_GpsEph* gpstail;
	Nav_GloEph* glotail;
	Nav_GalEph* galtail;
	Nav_BdsEph* bdstail;

	int Nav_Add_GpsEph(GpsEph* eph);
	int Nav_Add_GloEph(GloEph* eph);
	int Nav_Add_GalEph(GalEph* eph);
	int Nav_Add_BdsEph(BdsEph* eph);

	int Read_NavigationFile_Body(char* filename);
	int Read_NavigationFile(char* filename);

	int Clear();
};

struct PreHeaderSat
{
	int m_Precise;
	unsigned int m_Prn;
};

/// Precise file header information
class PreFileHeader
{
public:
	char m_Version;                           // Precise file version
	char m_FileType;                          // File type
	char m_DataType[6];                       // Data type
	char m_CoordSystem[6];                    // Coordinate System
	char m_OrbitType[4];                      // Orbit Type
	char m_Organization[5];                   // Organization
	unsigned int m_TotalEpochNum;             // Total Epoch Number
	unsigned int m_FirstEpoch_Week;           // First Epoch Week
	unsigned int m_StaCount;                  // Station Count
	double m_FirstEpoch_Sec;                  // First Epoch Second
	double m_Interval;                        // Interval
	double m_Mjd_u;
	double m_Mjd_m;

	PreHeaderSat m_Satllite[SAT_NUMBER];
	gnsstime m_FirstEpoch;

	int Read_PreciseFileHeader(char *filename);
};

class PreEpochData
{
public:
	double m_X, m_Y, m_Z;
	double n_Vx, m_Vy, m_Vz;
	double m_ClockCorr;

	PreEpochData *m_Previous;
	PreEpochData *m_Next;

	gnsstime m_GnssTime;

	unsigned int Pre_ana_Stringline_p(char *strline, unsigned int len);
};

class PreDataSat
{
public:
	unsigned int m_Prn;
	unsigned int m_EpochCount;         // satellite epoch count
	PreEpochData *m_FirstEpoch;
	PreEpochData *m_MidEpoch;
	PreEpochData *LastEpoch;
	PreEpochData **m_IndexTable;
	PreDataSat *m_Next;
};

typedef PreEpochData * PreEpochdata_ptr;

/// Precise Data
class PreciseData
{
public:
	PreciseData();
	~PreciseData();

	int m_BDS_IODE_flag;

	unsigned int m_PrnCount;
	PreFileHeader m_Header;
	NavData *m_NavData;
	PreDataSat *m_Predata;

	PreDataSat* Precise_Search_Prn(unsigned int prn);
	PreEpochData* Precise_Search_Nearest_Epoch(unsigned int prn, gnsstime &time);
	PreEpochData* Precise_Search_Nearest_Epoch_w(unsigned int prn, gnsstime &time);
	PreEpochData* Get_Precise_InterpData(unsigned int prn, gnsstime& time, int& len, int n);
	PreEpochData* Precise_Search_by_Epoch(unsigned int prn, gnsstime &time);

	int Precise_Interp_Pos_w(unsigned int prn, gnsstime& time, double pos[3], int n);
	int Precise_Interp_Vel_w(unsigned int prn, gnsstime& time, double vel[3], int n);
	int Precise_Interp_PosVel_w(unsigned int prn, gnsstime& time, double pos[3], double vel[3], int n);

	int Precise_Add_Data_p(unsigned int prn, PreEpochData *epochdata, PreDataSat **tail_sat);

	int Read_PreciseFileBody(char* filename);
	int Read_PreciseFile(char* filename, int flag);

	int Extract_PreciseClock(char* filename);
	int Clear();
};

/// Find Ephemeris by prn
Nav_GpsEph* Nav_GPS_FindPrn(unsigned int prn, NavData* navdata);
Nav_GloEph* Nav_GLO_FindPrn(unsigned int prn, NavData* navdata);
Nav_GalEph* Nav_GAL_FindPrn(unsigned int prn, NavData* navdata);
Nav_BdsEph* Nav_BDS_FindPrn(unsigned int prn, NavData* navdata);

void Pre_EpochData_ptr_Advance(PreEpochData **ptr, int count);

int Precise_Interp_Pos(gnsstime& time, PreEpochData* ptr, int len, double pos[3], int n);
int Precise_Interp_Vel(gnsstime& time, PreEpochData* ptr, int len, double vel[3], int n);

int Get_Nav_Orbitn(unsigned int prn, NavData* navdata);

#endif
