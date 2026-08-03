/* -------------------------------------------------------------------------
* SUPREME_Ionosphere.cpp : SUPREME Ionosphere Function
*
*           Copyright (C) by C. Zhao  All rights reserved
*           Contact caszcb@163.com
*
* Create : 2018.03.14
* ------------------------------------------------------------------------- */
#ifndef SUPREME_IONOSPHERE_H_HH
#define SUPREME_IONOSPHERE_H_HH

#include "SUPREME_GnssTime.h"
#include "SUPREME_Constant.h"
#include "SUPREME_Options.h"


/// Ion Setting
#define VAR_NOTEC               30.0*30.0      // Variance of no tec */
#define MIN_EL                  0.0            // Min elevation angle (rad) */
#define MIN_HGT                -1000.0         // Min user height (m) */

typedef struct {        /* TEC grid type */
	gtime_t time;       /* epoch time (GPST) */
	int ndata[3];       /* TEC grid data size {nlat,nlon,nhgt} */
	double rb;          /* earth radius (km) */
	double lats[3];     /* latitude start/interval (deg) */
	double lons[3];     /* longitude start/interval (deg) */
	double hgts[3];     /* heights start/interval (km) */
	double *data;       /* TEC grid data (tecu) */
	float *rms;         /* RMS values (tecu) */
	double *grid_lat;
	double *grid_lon;
} tec_t;


struct ionex_tec_t {     /* TEC grid type */
	gtime_t time;       /* epoch time (GPST) */
	int ndata[3];       /* TEC grid data size {nlat,nlon,nhgt} */
	double rb;          /* earth radius (km) */
	double lats[3];     /* latitude start/interval (deg) */
	double lons[3];     /* longitude start/interval (deg) */
	double hgts[3];     /* heights start/interval (km) */
	double* data;       /* TEC grid data (tecu) */
	double* rms;        /* RMS values (tecu) */
	double* grid_lat;
	double* grid_lon;
};

struct GIM_t
{
	int nt, ntmax;                       /* number of tec grid data */
	double cbias[GNSS_SATNO_NUM][3];     /* satellite dcb (0:p1-p2,1:p1-c1,2:p2-c2) (m) */
	ionex_tec_t *tec;
};



class IonexModel
{
public:
	IonexModel();
	~IonexModel();

	void SetGimFile(char *file);
	bool ReadIonexGimData();
	int iontec(gtime_t time, double ipp_b, double ipp_l, double *tec, double *var);
	double GetGridValue(int imap, int m, double &lat, double &lon, double &tec, double &rms);
	int GetGridCount(int imap);
	gtime_t GetMapTime(int imap);
	int GetMapCount();
	ionex_tec_t* ionex_tec;
	int nt, nt_max;
private:
	char reason[PATH_LENGTH];
	char ionex_file[PATH_LENGTH];

	

	FILE *fp_ionex;

	int getindex(double value, const double *range);
	int nitem(const double *range);
	int dataindex(int i, int j, int k, const int *ndata);
	ionex_tec_t* addtec(const double *lats, const double *lons, const double *hgts, double rb);
	double readionexh(double *lats, double *lons, double *hgts, double *rb, double *nexp);
	int readionexb(const double *lats, const double *lons, const double *hgts, double rb, double nexp);
	void combtec();
	int interptec(int k, const double *posp, double *value, double *rms);
	int interptec(int k, const double *posp, int isession, double *value, double *rms);
	int interptec(int k, const double *posp, const ionex_tec_t *tec, double *value, double *rms);
	int ion_tec_one_map(gtime_t time, const ionex_tec_t *tec, double ipp_b, double ipp_l, double *tec_sum, double *var);
};

double Meter2TECu(double iondelay, double gnssfreq);
double TECu2Meter(double Mtecu, double gnssfreq);

double IonMapf(const double *pos, const double *azel);

double IonPPP(const double *pos, const double *azel, double re, double hion, double *posp);

double IonModel_K8(gnsstime& t, const double *ion_A, const double *ion_B, const double *pos, const double& ele, const double& azi);

int IonModel_Grid(gtime_t time, double freq, const GIM_t *gim, const double *pos, const double *azel, int opt, double *delay, double *var, double *s_tec);

void ReadTecGridFile(char *file, GIM_t *gim, int opt);

#endif
