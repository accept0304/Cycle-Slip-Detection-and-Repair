/* -------------------------------------------------------------------------
* SUPREME_Tide.cpp : SUPREME Multi-GNSS Tide Corrections
*
*           Copyright (C) by C. Zhao
*           Contact caszcb@163.com
*
* Note: Solid Tide, Ocean Tide, Pole Tide
*
* Create : 2018.03.16
* ------------------------------------------------------------------------- */

#ifndef SUPREME_TIDE_H_HH
#define SUPREME_TIDE_H_HH

#include "SUPREME_GnssTime.h"

typedef struct {        /* earth rotation parameter data type */
	double mjd;         /* mjd (days) */
	double xp, yp;      /* pole offset (rad) */
	double xpr, ypr;    /* pole offset rate (rad/day) */
	double ut1_utc;     /* ut1-utc (s) */
	double lod;         /* length of day (s/day) */
} erpd_t;

typedef struct {        /* earth rotation parameter type */
	int n, nmax;        /* number and max number of data */
	erpd_t *data;       /* earth rotation parameter data */
} erp_t;

int ReadBlqFile(const char *file, const char *sta, double *odisp);

int ReadErpFile(const char *file, erp_t *erp);

void SolidTideEffect(const double *rsun, const double *rmoon, const double *pos, const double *E, double gmst, int opt, double *dr);

void OceanTideEffect(gtime_t tut, const double *odisp, double *denu);

void PoleTideEffect(gtime_t tut, const double *pos, const double *erpv, double *denu);

//int TideDisp(gnsstime &tutc, const double *rr, int opt, const double *odisp, double *dr);
int TideDisp(gnsstime &tutc, const double *rr, int opt, const erp_t *erp, const double *odisp, double *dr);

#endif
