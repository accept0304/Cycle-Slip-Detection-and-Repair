/* -------------------------------------------------------------------------
* SUPREME_Antenna.cpp : Satellite and Receiver Antenna Information
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2018.03.16
* ------------------------------------------------------------------------- */
#ifndef SUPREME_ANTENNA_H_HH
#define SUPREME_ANTENNA_H_HH

#include "SUPREME_GnssTime.h"

#define MAXANT        64                         // Max length of station name/antenna type
#define NSYS_USED     5                          // Number of used satellite systems, GPS/GLO/BDS/GAL/QZS
#define NFREQ         3                          // Number of frequency

/// BDS PCO Correction by Wuhan University
static const double BDS_PCO_WUM[35][3] = {
	{ 0.600, 0, 1.100 }, { 0.600, 0, 1.100 }, { 0.600, 0, 1.100 }, { 0.600, 0, 1.100 }, { 0.600, 0, 1.100 },
	{ 0.5864, 0, 2.5137 }, { 0.5864, 0, 2.7219 }, { 0.5864, 0, 3.4400 },
	{ 0.5864, 0, 3.5519 }, { 0.5864, 0, 4.0870 }, { 0.5750, 0, 1.9907 },
	{ 0.5750, 0, 2.2491 }, { 0.5750, 0, 2.0259 }, { 0.5750, 0, 2.1443 },
	{ 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 },
};

/// BDS PCO Correction by ESA
static const double BDS_PCO_ESA[35][3] = {
	{ 0.600, 0, 1.100 }, { 0.600, 0, 1.100 }, { 0.600, 0, 1.100 }, { 0.600, 0, 1.100 }, { 0.600, 0, 1.100 },
	{ 0.5490, 0, 3.0490 }, { 0.5490, 0, 3.2367 }, { 0.5490, 0, 3.8426 },
	{ 0.5490, 0, 3.9736 }, { 0.5490, 0, 3.8821 }, { 0.5490, 0, 2.0695 },
	{ 0.5490, 0, 2.3135 }, { 0.5490, 0, 2.2018 }, { 0.5490, 0, 2.3117 },
	{ 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 },
	{ 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 },
};

/// reference: rtklib
typedef struct {                                 // antenna parameter type
	int sat;                                     // satellite number (0:receiver)
	char type[MAXANT];                           // antenna type
	char code[MAXANT];                           // serial number or satellite code
	gtime_t ts, te;                              // valid time start and end
	double off[NSYS_USED*NFREQ][3];              // phase center offset e/n/u or x/y/z (m)
	double var[NSYS_USED*NFREQ][80 * 20];        // phase center variation (m)
	// el=90,85,...,0 or nadir=0,1,2,3,... (deg)
	double dazi;                                 // Increment of the azimuth£º0 to 360 with increment 'DAZI'(in degrees).
	double zen1, zen2, dzen;                     // Receiver antenna:Definition of the grid in zenith angle.
	// Satellite antenna:Definition of the grid in nadir angle.
} pcv_t;

typedef struct {                                 // antenna parameters type
	int n, nmax;                                 // number of data/allocated
	pcv_t *pcv;                                  // antenna parameters data
} pcvs_t;

// Phase Center Variance Data
class PCV_Data
{
public:
	PCV_Data();
	~PCV_Data();

	unsigned int m_ZenithNum;
	unsigned int m_AzimuthNum;
	unsigned int m_ZenithDiff;
	unsigned int m_AzimuthDiff;

	double **m_PCV1;
	double **m_PCV2;

	double *m_IndexCol;
	double *m_IndexRow;
};

/// Satellite Antenna Infomation
class AntennaInfo_Sat
{
public:
	AntennaInfo_Sat();
	~AntennaInfo_Sat();

	char m_Ant_Type[21];                 // Satellite Antenna Type
	char m_Sat_Type[4];                  // Satellite Type
	unsigned int m_Prn;                  // Satellite PRN

	AntennaInfo_Sat *m_Next;

	double m_PCO[3];                     // dx, dy, dz

	PCV_Data m_PCV;

	int Clear();
};

/// Satellite Antenna Data
class SatelliteAntenna
{
public:
	SatelliteAntenna();
	~SatelliteAntenna();

	unsigned int m_PrnCount;             // Satellite prn counts

	AntennaInfo_Sat* m_AntInfo;
	AntennaInfo_Sat** m_IndexTable;

	int Clear();
};

/// Receiver Antenna Data
class ReceiverAntenna
{
public:
	ReceiverAntenna();
	~ReceiverAntenna();

	char m_AntennaType[21];             // Receiver Antenna type

	double m_PCO1[3];                   // Phase Center Offset on L1
	double m_PCO2[3];                   // Phase Center Offset on L2

	PCV_Data m_PCV;                     // Phase Center Variable
};

/// Read antx file: igs08.atx  igs14.atx
int ReadAntex(const char *file, pcvs_t *pcvs);

pcv_t *searchpcv(int sat, const char *type, gtime_t time, const pcvs_t *pcvs);

/// Read Antenna Infomation File
int Read_Satellite_Antenna_PCO(char *SatAntFile, double jd, SatelliteAntenna& sat_ant_info);
int Read_Satellite_Antenna_PCV(char *PhaseFile, SatelliteAntenna& sat_ant_info);
int Read_Receiver_Antenna(char* PhaseFile, char* ant_type, ReceiverAntenna& rec_ant_info);

AntennaInfo_Sat* Find_SatAntennaInfo_by_Prn(SatelliteAntenna& SatAntInfo, unsigned int prn);

/// Satellite antenna PCO and PCV
void SatAntPCV(char sys, int sat, const double *rs, const double *rr, const pcv_t *pcv, double *dant);
void SatAntPCO(pcv_t *pcvs, const double *rsun, const double *rs, char sys, int sat,int orb_n, double *dant);
int SatAntPCO(char sys, int sat, double satpos[3], double sunpos[3], pcv_t *pcvs, double SatPCO[3]);
int SatAntPCO(double satpos[3], double sunpos[3], AntennaInfo_Sat& sat_info, double SatPCO[3]);
int SatAntPCO_BDS(unsigned int prn, double satpos[3], double sunpos[3], double SatPCO[3]);
int SatAntPCV(double sta[3], double sat[3], double sun[3], AntennaInfo_Sat& sat_info, unsigned int freq, double& pcv);

// Receiver anntenna PCO and PCV
void RecAntModel(char sys, const pcv_t *pcv, const double *del, const double *azel, int opt, double *dant);
int RecAntPCO(double a, double e, double sta[3], double ant[3], ReceiverAntenna& rec_info, double dv1[3], double dv2[3]);
int RecAntPCV(double a, double e, double sta[3], double sat[3], ReceiverAntenna& rec_info, unsigned int freq, double& pcv);

int RecAntOffset(double a, double e, double sta[3], double ant[3], double dv[3]);

int WindupCorr(double rsun[3], double rsat[3], double rsta[3], double& phw);

#endif
