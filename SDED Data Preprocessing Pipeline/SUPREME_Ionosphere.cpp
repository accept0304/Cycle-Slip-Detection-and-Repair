/* -------------------------------------------------------------------------
* SUPREME_Ionosphere.cpp : SUPREME Ionosphere Function
*
*           Copyright (C) by C. Zhao 2016-2018 All rights reserved
*           Contact zhaoC.@whigg.ac.cn
*
* Create : 2018.03.14
* ------------------------------------------------------------------------- */
#include <math.h>

#include "SUPREME_Ionosphere.h"
#include "SUPREME_CommonFunction.h"

/// Get index
int getindex(double value, const double *range)
{
	if (range[2] == 0.0) return 0;
	if (range[1]>0.0 && (value<range[0] || range[1]<value)) return -1;
	if (range[1]<0.0 && (value<range[1] || range[0]<value)) return -1;
	return (int)floor((value - range[0]) / range[2] + 0.5);
}

/// Get number of items
int nitem(const double *range)
{
	return getindex(range[1], range) + 1;
}

/// Data index (i:lat,j:lon,k:hgt)
int dataindex(int i, int j, int k, const int *ndata)
{
	if (i<0 || ndata[0] <= i || j<0 || ndata[1] <= j || k<0 || ndata[2] <= k) return -1;
	return i + ndata[0] * (j + ndata[1] * k);
}

/* -----------------------------------------------------------------------------------
* Name: Meter2TECu
* Description : convert meter to TECu
* Parameters  :
*		double iondelay		I		ionospheric delay [unit: meter]
*		double gnssfreq		I		common for GPS L1, BDS B1 or Galileo E1 frequency
* ----------------------------------------------------------------------------------- */
double Meter2TECu(double iondelay, double gnssfreq)
{
	double Mtecu = 0.0;
	Mtecu = iondelay*pow(gnssfreq, 2) / ION_FACTEC;
	return Mtecu;
}

double TECu2Meter(double Mtecu, double gnssfreq)
{
	double meter = 0.0;
	meter = Mtecu / pow(gnssfreq, 2) * ION_FACTEC;

	return meter;
}

/* ionosphere mapping function -------------------------------------------
* compute ionospheric delay mapping function by single layer model
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric mapping function
* ----------------------------------------------------------------------- */
double IonMapf(const double *pos, const double *azel)
{
	if (pos[2] >= HION) return 1.0;
	return 1.0 / cos(asin((WGS84_A + pos[2]) / (WGS84_A + HION)*sin(PI / 2.0 - azel[1])));
}

/* ionospheric pierce point position -------------------------------------------
* compute ionospheric pierce point (ipp) position and slant factor
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double re        I   earth radius (km)
*          double hion      I   altitude of ionosphere (km)
*          double *posp     O   pierce point position {lat,lon,h} (rad,m)
* return : slant factor
* notes  : see ref [2], only valid on the earth surface
*          fixing bug on ref [2] A.4.4.10.1 A-22,23
* ----------------------------------------------------------------------------- */
double IonPPP(const double *pos, const double *azel, double re, double hion, double *posp)
{
	double cosaz, rp, ap, sinap, tanap;

	rp = re / (re + hion)*cos(azel[1]);
	ap = PI / 2.0 - azel[1] - asin(rp);
	sinap = sin(ap);
	tanap = tan(ap);
	cosaz = cos(azel[0]);
	posp[0] = asin(sin(pos[0])*cos(ap) + cos(pos[0])*sinap*cosaz);

	if ((pos[0]> 70.0*D2R&& tanap*cosaz>tan(PI / 2.0 - pos[0])) ||
		(pos[0]<-70.0*D2R&&-tanap*cosaz>tan(PI / 2.0 + pos[0]))) {
		posp[1] = pos[1] + PI - asin(sinap*sin(azel[0]) / cos(posp[0]));
	}
	else {
		posp[1] = pos[1] + asin(sinap*sin(azel[0]) / cos(posp[0]));
	}
	return 1.0 / sqrt(1.0 - rp*rp);
}

/// Get ionospheric delay by broadcast ionosphere model (Klobuchar Model)
double IonModel_K8(gnsstime& t, const double *ion_A, const double *ion_B, const double *pos, const double& ele, const double& azi)
{
	double tt, f, psi, phi, lam, amp, per, x;
	int week;

	//if (pos[2]<-1E3 || ele <= 0) return 0.0;

	/* earth centered angle (semi-circle) */
	psi = 0.0137 / (ele / PI + 0.11) - 0.022;

	/* subionospheric latitude/longitude (semi-circle) */
	phi = pos[0] / PI + psi*cos(azi);
	if (phi> 0.416) phi = 0.416;
	else if (phi<-0.416) phi = -0.416;
	lam = pos[1] / PI + psi*sin(azi) / cos(phi*PI);

	/* geomagnetic latitude (semi-circle) */
	phi += 0.064*cos((lam - 1.617)*PI);

	/* local time (s) */
	tt = 43200.0*lam + t.m_SecondofWeek;
	tt -= floor(tt / 86400.0)*86400.0; /* 0<=tt<86400 */

	/* slant factor */
	f = 1.0 + 16.0*pow(0.53 - ele / PI, 3.0);

	/* ionospheric delay */
	amp = ion_A[0] + phi*(ion_A[1] + phi*(ion_A[2] + phi*ion_A[3]));
	per = ion_B[0] + phi*(ion_B[1] + phi*(ion_B[2] + phi*ion_B[3]));
	amp = amp<    0.0 ? 0.0 : amp;
	per = per<72000.0 ? 72000.0 : per;
	x = 2.0*PI*(tt - 50400.0) / per;

	return LIGHTSPEED*f*(fabs(x)<1.57 ? 5E-9 + amp*(1.0 + x*x*(-0.5 + x*x / 24.0)) : 5E-9);
}


//// Get ionosphere by GIM ============================================

/// Read Ionex dcb aux data
void readionexdcb(FILE *fp, double *dcb, double *rms)
{
	int i, sat;
	char buff[1024], id[32], *label;

	for (i = 0; i<GNSS_SATNO_NUM; i++) dcb[i] = rms[i] = 0.0;

	while (fgets(buff, sizeof(buff), fp)) {
		if (strlen(buff)<60) continue;
		label = buff + 60;

		if (strstr(label, "PRN / BIAS / RMS") == label) {

			strncpy(id, buff + 3, 3); id[3] = '\0';

			if (!(sat = GetSatNo(id, 3, 0)))
				continue;

			dcb[sat - 1] = str2num(buff, 6, 10);
			rms[sat - 1] = str2num(buff, 16, 10);
		}
		else if (strstr(label, "END OF AUX DATA") == label) break;
	}
}

/// Read Ionex header
double readionexh(FILE *fp, double *lats, double *lons, double *hgts,
	double *rb, double *nexp, double *dcb, double *rms)
{
	double ver = 0.0;
	char buff[1024], *label;

	while (fgets(buff, sizeof(buff), fp)) {

		if (strlen(buff)<60) continue;
		label = buff + 60;

		if (strstr(label, "IONEX VERSION / TYPE") == label) {
			if (buff[20] == 'I') ver = str2num(buff, 0, 8);
		}
		else if (strstr(label, "BASE RADIUS") == label) {
			*rb = str2num(buff, 0, 8);
		}
		else if (strstr(label, "HGT1 / HGT2 / DHGT") == label) {
			hgts[0] = str2num(buff, 2, 6);
			hgts[1] = str2num(buff, 8, 6);
			hgts[2] = str2num(buff, 14, 6);
		}
		else if (strstr(label, "LAT1 / LAT2 / DLAT") == label) {
			lats[0] = str2num(buff, 2, 6);
			lats[1] = str2num(buff, 8, 6);
			lats[2] = str2num(buff, 14, 6);
		}
		else if (strstr(label, "LON1 / LON2 / DLON") == label) {
			lons[0] = str2num(buff, 2, 6);
			lons[1] = str2num(buff, 8, 6);
			lons[2] = str2num(buff, 14, 6);
		}
		else if (strstr(label, "EXPONENT") == label) {
			*nexp = str2num(buff, 0, 6);
		}
		else if (strstr(label, "START OF AUX DATA") == label&&
			strstr(buff, "DIFFERENTIAL CODE BIASES")) {
			readionexdcb(fp, dcb, rms);
		}
		else if (strstr(label, "END OF HEADER") == label) {
			return ver;
		}
	}
	return 0.0;
}

/// Add tec data to navigation data
ionex_tec_t *addtec(const double *lats, const double *lons, const double *hgts,
	double rb, GIM_t *gim)
{
	ionex_tec_t *p, *gim_tec;
	gtime_t time0 = { 0 };
	int i, n, ndata[3];

	ndata[0] = nitem(lats);
	ndata[1] = nitem(lons);
	ndata[2] = nitem(hgts);
	if (ndata[0] <= 1 || ndata[1] <= 1 || ndata[2] <= 0) return NULL;

	if (gim->nt >= gim->ntmax)
	{
		gim->ntmax += 256;
		if (!(gim_tec = (ionex_tec_t*)realloc(gim->tec, sizeof(ionex_tec_t)*gim->ntmax)))
		{
			free(gim->tec); gim->tec = NULL; gim->nt = gim->ntmax = 0;
			return NULL;
		}
		gim->tec = gim_tec;
	}
	p = gim->tec + gim->nt;
	p->time = time0;
	p->rb = rb;
	for (i = 0; i<3; i++)
	{
		p->ndata[i] = ndata[i];
		p->lats[i] = lats[i];
		p->lons[i] = lons[i];
		p->hgts[i] = hgts[i];
	}
	n = ndata[0] * ndata[1] * ndata[2];

	if (!(p->data = (double *)malloc(sizeof(double)*n)) ||
		!(p->rms = (double*)malloc(sizeof(double)*n)))
	{
		return NULL;
	}
	for (i = 0; i<n; i++)
	{
		p->data[i] = 0.0;
		p->rms[i] = 0.0f;
	}
	gim->nt++;
	return p;
}

/// Read Ionex body
int readionexb(FILE *fp, const double *lats, const double *lons,
	const double *hgts, double rb, double nexp, GIM_t *gim)
{
	ionex_tec_t *p = NULL;
	gtime_t time = { 0 };
	double lat, lon[3], hgt, x;
	int i, j, k, n, m, index, type = 0;
	char buff[1024], *label = buff + 60;

	while (fgets(buff, sizeof(buff), fp))
	{
		if (strlen(buff)<60) continue;

		if (strstr(label, "START OF TEC MAP") == label)
		{
			if ((p = addtec(lats, lons, hgts, rb, gim))) type = 1;//出问题！！
		}
		else if (strstr(label, "END OF TEC MAP") == label)
		{
			type = 0;
			p = NULL;
		}
		else if (strstr(label, "START OF RMS MAP") == label)
		{
			type = 2;
			p = NULL;
		}
		else if (strstr(label, "END OF RMS MAP") == label)
		{
			type = 0;
			p = NULL;
		}
		else if (strstr(label, "EPOCH OF CURRENT MAP") == label)
		{
			if (str2time(buff, 0, 36, &time))
				continue;

			if (type == 2)
			{
				for (i = gim->nt - 1; i >= 0; i--)
				{
					if (fabs(timediff(time, gim->tec[i].time)) >= 1.0) continue;
					p = gim->tec + i;
					break;
				}
			}
			else if (p) p->time = time;
		}
		else if (strstr(label, "LAT/LON1/LON2/DLON/H") == label&&p)
		{
			lat = str2num(buff, 2, 6);
			lon[0] = str2num(buff, 8, 6);
			lon[1] = str2num(buff, 14, 6);
			lon[2] = str2num(buff, 20, 6);
			hgt = str2num(buff, 26, 6);

			i = getindex(lat, p->lats);
			k = getindex(hgt, p->hgts);
			n = nitem(lon);

			for (m = 0; m<n; m++) {
				if (m % 16 == 0 && !fgets(buff, sizeof(buff), fp)) break;

				j = getindex(lon[0] + lon[2] * m, p->lons);
				if ((index = dataindex(i, j, k, p->ndata))<0) continue;

				if ((x = str2num(buff, m % 16 * 5, 5)) == 9999.0) continue;

 				if (type == 1) p->data[index] = x*pow(10.0, nexp);
				else p->rms[index] = (float)(x*pow(10.0, nexp));
			}

		}
	}
	return 1;
}

/// Combine tec grid data
void combtec(GIM_t *gim)
{
	ionex_tec_t tmp;
	int i, j, n = 0;

	for (i = 0; i<gim->nt - 1; i++)
	{
		for (j = i + 1; j<gim->nt; j++)
		{
			if (timediff(gim->tec[j].time, gim->tec[i].time)<0.0)
			{
				tmp = gim->tec[i];
				gim->tec[i] = gim->tec[j];
				gim->tec[j] = tmp;
			}
		}
	}
	for (i = 0; i<gim->nt; i++) {
		if (i>0 && timediff(gim->tec[i].time, gim->tec[n - 1].time) == 0.0)
		{
			free(gim->tec[n - 1].data);
			free(gim->tec[n - 1].rms);
			gim->tec[n - 1] = gim->tec[i];
			continue;
		}
		gim->tec[n++] = gim->tec[i];
	}
	gim->nt = n;
}

/* Read Ionex Tec Grid File ------------------------------------------------------
* Read ionex ionospheric tec grid file
* args   : char   *file       I   ionex tec grid file
*                                 (wind-card * is expanded)
*          GIM_t  *gim        IO  GIM data type
*                                 nav->nt, nav->ntmax and nav->tec are modified
*          int    opt         I   read option (1: no clear of tec data,0:clear)
* return : none
* notes  : see ref [1]
* ------------------------------------------------------------------------------- */
void ReadTecGridFile(char *file, GIM_t *gim, int opt)
{
	FILE *fp;
	double lats[3] = { 0 }, lons[3] = { 0 }, hgts[3] = { 0 }, rb = 0.0, nexp = -1.0;
	double dcb[GNSS_SATNO_NUM] = { 0 }, rms[GNSS_SATNO_NUM] = { 0 };
	int i = 0, n = 1;;

	//% clear of tec grid data option
	if ((!opt) && (gim->tec))
	{
		if (gim->tec) free(gim->tec); gim->tec = NULL; gim->nt = gim->ntmax = 0;
	}

	for (i = 0; i<n; i++)
	{
		if (!(fp = fopen(file, "r")))
			continue;

		//% read ionex header
		if (readionexh(fp, lats, lons, hgts, &rb, &nexp, dcb, rms) <= 0.0)
			continue;

		//% read ionex body
		readionexb(fp, lats, lons, hgts, rb, nexp, gim);

		fclose(fp);
	}

	//% combine tec grid data
	if (gim->nt>0) combtec(gim);

	//% P1-P2 dcb
	for (i = 0; i<GNSS_SATNO_NUM; i++)
	{
		gim->cbias[i][0] = LIGHTSPEED*dcb[i] * 1E-9; /* ns->m */
	}
}

/// Interpolate tec grid data
int interptec(const ionex_tec_t *tec, int k, const double *posp, double *value, double *rms)
{
	double dlat, dlon, a, b, d[4] = { 0 }, r[4] = { 0 };
	int i, j, n, index;

	*value = *rms = 0.0;

	if (tec->lats[2] == 0.0 || tec->lons[2] == 0.0) return 0;

	dlat = posp[0] * R2D - tec->lats[0];
	dlon = posp[1] * R2D - tec->lons[0];
	if (tec->lons[2]>0.0) dlon -= floor(dlon / 360)*360.0; /*  0<=dlon<360 */
	else                  dlon += floor(-dlon / 360)*360.0; /* -360<dlon<=0 */

	a = dlat / tec->lats[2];
	b = dlon / tec->lons[2];
	i = (int)floor(a); a -= i;
	j = (int)floor(b); b -= j;

	//% get gridded tec data
	for (n = 0; n<4; n++)
	{
		if ((index = dataindex(i + (n % 2), j + (n<2 ? 0 : 1), k, tec->ndata))<0) continue;
		d[n] = tec->data[index];
		r[n] = tec->rms[index];
	}
	if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0)
	{

		//% bilinear interpolation (inside of grid)
		*value = (1.0 - a)*(1.0 - b)*d[0] + a*(1.0 - b)*d[1] + (1.0 - a)*b*d[2] + a*b*d[3];
		*rms = (1.0 - a)*(1.0 - b)*r[0] + a*(1.0 - b)*r[1] + (1.0 - a)*b*r[2] + a*b*r[3];
	}
	//% nearest-neighbour extrapolation (outside of grid)
	else if (a <= 0.5&&b <= 0.5&&d[0]>0.0) { *value = d[0]; *rms = r[0]; }
	else if (a> 0.5&&b <= 0.5&&d[1]>0.0) { *value = d[1]; *rms = r[1]; }
	else if (a <= 0.5&&b> 0.5&&d[2]>0.0) { *value = d[2]; *rms = r[2]; }
	else if (a> 0.5&&b> 0.5&&d[3]>0.0) { *value = d[3]; *rms = r[3]; }
	else
	{
		i = 0;
		for (n = 0; n<4; n++) if (d[n]>0.0) { i++; *value += d[n]; *rms += r[n]; }
		if (i == 0) return 0;
		*value /= i; *rms /= i;
	}
	return 1;
}

/// Ionosphere delay by tec grid data
int IonDelay_GIM(gtime_t time, double freq, const ionex_tec_t *tec, const double *pos,
	const double *azel, int opt, double *delay, double *var, double* s_tec1)
{
	const double fact = 40.30E16 / freq / freq; /* tecu->L1 iono (m) */
	double fs, posp[3] = { 0 }, vtec, rms, hion, rp;
	int i;

	*delay = *var = *s_tec1 = 0.0;

	for (i = 0; i<tec->ndata[2]; i++)  //% for a layer
	{

		hion = tec->hgts[0] + tec->hgts[2] * i;

		//% ionospheric pierce point position
		fs = IonPPP(pos, azel, tec->rb, hion, posp);

		if (opt & 2)
		{
			//% modified single layer mapping function (M-SLM) ref [2]
			rp = tec->rb / (tec->rb + hion)*sin(0.9782*(PI / 2.0 - azel[1]));
			fs = 1.0 / sqrt(1.0 - rp*rp);
		}
		if (opt & 1)
		{
			//% earth rotation correction (sun-fixed coordinate)
			posp[1] += 2.0*PI*timediff(time, tec->time) / 86400.0;
		}
		//% interpolate tec grid data
		if (!interptec(tec, i, posp, &vtec, &rms)) return 0;

		*s_tec1 += fs * vtec * 10E16;
		*delay += fact * fs * vtec;
		*var += fact * fact * fs * fs * rms * rms;

	}

	return 1;
}

/* Ionosphere model by tec grid data --------------------------------------------
* compute ionospheric delay by tec grid data
* args   : gtime_t time     I   time (gpst)
*          GIM_t  *gim      I   GIM data
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    opt       I   model option
*                                bit0: 0:earth-fixed,1:sun-fixed
*                                bit1: 0:single-layer,1:modified single-layer
*          double *delay    O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
* notes  : before calling the function, read tec grid data by calling readtec()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
* ------------------------------------------------------------------------------ */
int IonModel_Grid(gtime_t time, double freq, const GIM_t *gim, const double *pos,
	const double *azel, int opt, double *delay, double *var, double *s_tec)
{
	double dels[2], vars[2], s_tec1[2], a, tt;
	int i, stat[2];

	if (azel[1]<MIN_EL || pos[2]<MIN_HGT)
	{
		*delay = 0.0;
		*var = VAR_NOTEC;
		*s_tec = 0.0;
		return 1;
	}
	for (i = 0; i<gim->nt; i++)
	{
		if (timediff(gim->tec[i].time, time)>0.0) break;
	}
	if (i == 0 || i >= gim->nt)
		return 0;

	if ((tt = timediff(gim->tec[i].time, gim->tec[i - 1].time)) == 0.0)
		return 0;

	//% ionospheric delay by tec grid data
	stat[0] = IonDelay_GIM(time, freq, gim->tec + i - 1, pos, azel, opt, dels, vars, s_tec1);
	stat[1] = IonDelay_GIM(time, freq, gim->tec + i, pos, azel, opt, dels + 1, vars + 1, s_tec1 + 1);

	if (!stat[0] && !stat[1])
		return 0;

	if (stat[0] && stat[1])   //% linear interpolation by time
	{
		a = timediff(time, gim->tec[i - 1].time) / tt;
		*delay = dels[0] * (1.0 - a) + dels[1] * a;
		*var = vars[0] * (1.0 - a) + vars[1] * a;
		*s_tec = s_tec1[0] * (1.0 - a) + s_tec1[1] * a;
	}
	else if (stat[0])         //% nearest-neighbour extrapolation by time
	{
		*delay = dels[0];
		*var = vars[0];
		*s_tec = s_tec1[0];
	}
	else
	{
		*delay = dels[1];
		*var = vars[1];
		*s_tec = s_tec1[1];
	}

	return 1;
}

IonexModel::IonexModel()
{
	memset(this->reason, 0, sizeof(this->reason));
	this->ionex_tec = NULL;
	this->fp_ionex = NULL;
	memset(this->ionex_file, 0, sizeof(this->ionex_file));

	this->nt = 0;
	this->nt_max = 0;
}

IonexModel::~IonexModel()
{
	int imap;

	for (imap = 0; imap<this->nt_max; imap++)
	{
		if (this->ionex_tec[imap].data != NULL)
		{
			free(this->ionex_tec[imap].data);
			this->ionex_tec[imap].data = NULL;
		}
		if (this->ionex_tec[imap].rms != NULL)
		{
			free(this->ionex_tec[imap].data);
			this->ionex_tec[imap].data = NULL;
		}
		if (this->ionex_tec[imap].grid_lat != NULL)
		{
			free(this->ionex_tec[imap].grid_lat);
			this->ionex_tec[imap].grid_lat = NULL;
		}
		if (this->ionex_tec[imap].grid_lon != NULL)
		{
			free(this->ionex_tec[imap].grid_lon);
			this->ionex_tec[imap].grid_lon = NULL;
		}

	}

	if (this->ionex_tec != NULL)
		free(this->ionex_tec);

}

void IonexModel::SetGimFile(char *file)
{
	memset(this->ionex_file, 0, sizeof(this->ionex_file));
	strcpy(this->ionex_file, file);

	if (this->fp_ionex != NULL)
		free(this->fp_ionex);

	return;
}

/// Read Ionex GIM data (by LM)
bool IonexModel::ReadIonexGimData()
{
	char subroutine[20];
	double lats[3] = { 0 }, lons[3] = { 0 }, hgts[3] = { 0 }, rb = 0.0, nexp = -1.0;

	memset(subroutine, 0, sizeof(subroutine));
	strcpy(subroutine, "ReadIonexGimData");

	this->fp_ionex = fopen(this->ionex_file, "rt");

	if (this->fp_ionex == NULL)
	{
		memset(this->reason, 0, sizeof(this->reason));
		sprintf(this->reason, "IONEX file open failed.FileName:%-s", this->ionex_file);
		return false;
	}

	/* read ionex header */
	if (readionexh(lats, lons, hgts, &rb, &nexp) <= 0.0)
	{
		memset(this->reason, 0, sizeof(this->reason));
		sprintf(this->reason, "IONEX file read failed.FileName:%-s", this->ionex_file);
		if (this->fp_ionex != NULL) fclose(this->fp_ionex);
		return false;
	}
	/* read ionex body */
	if (this->readionexb(lats, lons, hgts, rb, nexp) != 1)
	{
		if (this->fp_ionex != NULL)
			fclose(this->fp_ionex);
		return false;

	}

	if (this->fp_ionex != NULL)
		fclose(this->fp_ionex);

	return true;
}

/* get index -----------------------------------------------------------------*/
int IonexModel::getindex(double value, const double *range)
{
	if (range[2] == 0.0) return 0;
	if (range[1]>0.0 && (value<range[0] || range[1]<value)) return -1;
	if (range[1]<0.0 && (value<range[1] || range[0]<value)) return -1;
	return (int)floor((value - range[0]) / range[2] + 0.5);
}

/* get number of items -------------------------------------------------------*/
int IonexModel::nitem(const double *range)
{
	return getindex(range[1], range) + 1;
}

/* data index (i:lat,j:lon,k:hgt) --------------------------------------------*/
int IonexModel::dataindex(int i, int j, int k, const int *ndata)
{
	if (i<0 || ndata[0] <= i || j<0 || ndata[1] <= j || k<0 || ndata[2] <= k) return -1;
	return i + ndata[0] * (j + ndata[1] * k);
}

/* add tec data to navigation data -------------------------------------------*/
ionex_tec_t* IonexModel::addtec(const double *lats, const double *lons, const double *hgts,
	double rb)
{
	gtime_t time0 = { 0 };
	int i, n, ndata[3];

	ndata[0] = nitem(lats);
	ndata[1] = nitem(lons);
	ndata[2] = nitem(hgts);
	if (ndata[0] <= 1 || ndata[1] <= 1 || ndata[2] <= 0) return NULL;

	if (this->nt >= this->nt_max)
	{
		this->nt_max += 1;
		this->ionex_tec = (ionex_tec_t *)realloc(this->ionex_tec, this->nt_max * sizeof(ionex_tec_t));
		if (this->ionex_tec == NULL)
		{
			memset(this->reason, 0, sizeof(this->reason));
			sprintf(this->reason, "ionex_tec malloc error ntmax=%d", this->nt_max);
			free(this->ionex_tec);
			this->ionex_tec = NULL;
			this->nt = this->nt_max = 0;
			return NULL;
		}
		else
		{
			for (i = 1; i <= 1; i++)
			{
				memset(this->ionex_tec[this->nt_max - i].hgts, 0, 3 * sizeof(double));
				memset(this->ionex_tec[this->nt_max - i].lats, 0, 3 * sizeof(double));
				memset(this->ionex_tec[this->nt_max - i].lons, 0, 3 * sizeof(double));
				memset(this->ionex_tec[this->nt_max - i].ndata, 0, 3 * sizeof(int));
				this->ionex_tec[this->nt_max - i].rb = 0.0;
				this->ionex_tec[this->nt_max - i].rms = NULL;
				this->ionex_tec[this->nt_max - i].data = NULL;
				this->ionex_tec[this->nt_max - i].time = time0;

			}

		}
	}
	this->ionex_tec[this->nt].time = time0;
	this->ionex_tec[this->nt].rb = rb;
	for (i = 0; i<3; i++)
	{
		this->ionex_tec[this->nt].ndata[i] = ndata[i];
		this->ionex_tec[this->nt].lats[i] = lats[i];
		this->ionex_tec[this->nt].lons[i] = lons[i];
		this->ionex_tec[this->nt].hgts[i] = hgts[i];
	}
	n = ndata[0] * ndata[1] * ndata[2];

	if (!(this->ionex_tec[this->nt].data = (double*)malloc(sizeof(double)*n)) ||
		!(this->ionex_tec[this->nt].rms = (double*)malloc(sizeof(double)*n)) ||
		!(this->ionex_tec[this->nt].grid_lat = (double *)malloc(sizeof(double)*n)) ||
		!(this->ionex_tec[this->nt].grid_lon = (double *)malloc(sizeof(double)*n)))
	{
		return NULL;
	}

	memset(this->ionex_tec[this->nt].data, 0, sizeof(double)*n);
	memset(this->ionex_tec[this->nt].rms, 0, sizeof(double)*n);
	memset(this->ionex_tec[this->nt].grid_lat, 0, sizeof(double)*n);
	memset(this->ionex_tec[this->nt].grid_lon, 0, sizeof(double)*n);

	return this->ionex_tec + (this->nt++);
}

/* read ionex header ---------------------------------------------------------*/
double IonexModel::readionexh(double *lats, double *lons, double *hgts, double *rb, double *nexp)
{
	double ver = 0.0;
	char buff[1024], *label;


	while (fgets(buff, sizeof(buff), this->fp_ionex))
	{

		if (strlen(buff)<60) continue;
		label = buff + 60;

		if (strstr(label, "IONEX VERSION / TYPE") == label)
		{
			if (buff[20] == 'I') ver = str2num(buff, 0, 8);
		}
		else if (strstr(label, "BASE RADIUS") == label)
		{
			*rb = str2num(buff, 0, 8);
		}
		else if (strstr(label, "HGT1 / HGT2 / DHGT") == label)
		{
			hgts[0] = str2num(buff, 2, 6);
			hgts[1] = str2num(buff, 8, 6);
			hgts[2] = str2num(buff, 14, 6);
		}
		else if (strstr(label, "LAT1 / LAT2 / DLAT") == label)
		{
			lats[0] = str2num(buff, 2, 6);
			lats[1] = str2num(buff, 8, 6);
			lats[2] = str2num(buff, 14, 6);
		}
		else if (strstr(label, "LON1 / LON2 / DLON") == label)
		{
			lons[0] = str2num(buff, 2, 6);
			lons[1] = str2num(buff, 8, 6);
			lons[2] = str2num(buff, 14, 6);
		}
		else if (strstr(label, "EXPONENT") == label)
		{
			*nexp = str2num(buff, 0, 6);
		}
		else if (strstr(label, "END OF HEADER") == label)
		{
			return ver;
		}
	}
	return 0.0;
}

/* read ionex body -----------------------------------------------------------*/
int IonexModel::readionexb(const double *lats, const double *lons, const double *hgts, double rb, double nexp)
{
	ionex_tec_t *p = NULL;
	gtime_t time = { 0 };
	double lat, lon[3], hgt, x;
	int i, j, k, n, m, index, type = 0;
	char buff[1024], *label = buff + 60;
	char subroutine[20];

	memset(subroutine, 0, sizeof(subroutine));
	strcpy(subroutine, "readionexb");


	while (fgets(buff, sizeof(buff), this->fp_ionex))
	{

		if (strlen(buff)<60) continue;

		if (strstr(label, "START OF TEC MAP") == label)
		{
			if ((p = addtec(lats, lons, hgts, rb)))
				type = 1;
		}
		else if (strstr(label, "END OF TEC MAP") == label) {
			type = 0;
			p = NULL;
		}
		else if (strstr(label, "START OF RMS MAP") == label) {
			type = 2;
			p = NULL;
		}
		else if (strstr(label, "END OF RMS MAP") == label)
		{
			type = 0;
			p = NULL;
		}
		else if (strstr(label, "EPOCH OF CURRENT MAP") == label)
		{
			if (str2time(buff, 0, 36, &time))
			{
				memset(this->reason, 0, sizeof(this->reason));
				sprintf(this->reason, "IONEX epoch invalid: %-36.36s", buff);
				continue;
			}
			if (type == 2)
			{
				for (i = this->nt - 1; i >= 0; i--)
				{
					if (fabs(timediff(time, this->ionex_tec[i].time)) >= 1.0)
						continue;
					p = this->ionex_tec + i;
					break;
				}
			}
			else if (p)
				p->time = time;
		}
		else if (strstr(label, "LAT/LON1/LON2/DLON/H") == label&&p)
		{
			lat = str2num(buff, 2, 6);
			lon[0] = str2num(buff, 8, 6);
			lon[1] = str2num(buff, 14, 6);
			lon[2] = str2num(buff, 20, 6);
			hgt = str2num(buff, 26, 6);

			i = getindex(lat, p->lats);
			k = getindex(hgt, p->hgts);
			n = nitem(lon);

			for (m = 0; m<n; m++)
			{
				if (m % 16 == 0 && !fgets(buff, sizeof(buff), this->fp_ionex))
					break;

				j = getindex(lon[0] + lon[2] * m, p->lons);
				if ((index = dataindex(i, j, k, p->ndata))<0)
					continue;

				if ((x = str2num(buff, m % 16 * 5, 5)) == 9999.0)
					continue;

				if (x<5.1 || x > 1800.0)
					continue;

				if (type == 1)
				{
					p->data[index] = x*pow(10.0, nexp);
					p->grid_lat[index] = lat;
					p->grid_lon[index] = lon[0] + m*lon[2];
				}
				else
					p->rms[index] = (float)(x*pow(10.0, nexp));
			}
		}
	}
	return 1;
}

/* combine tec grid data -----------------------------------------------------*/
void IonexModel::combtec()
{
	ionex_tec_t tmp;
	int i, j, n = 0;

	for (i = 0; i<this->nt - 1; i++) {
		for (j = i + 1; j<this->nt; j++)
		{
			if (timediff(this->ionex_tec[j].time, this->ionex_tec[i].time)<0.0)
			{
				tmp = this->ionex_tec[i];
				this->ionex_tec[i] = this->ionex_tec[j];
				this->ionex_tec[j] = tmp;
			}
		}
	}
	for (i = 0; i<this->nt; i++)
	{
		if (i>0 && timediff(this->ionex_tec[i].time, this->ionex_tec[n - 1].time) == 0.0)
		{
			free(this->ionex_tec[n - 1].data);
			free(this->ionex_tec[n - 1].rms);
			this->ionex_tec[n - 1] = this->ionex_tec[i];
			continue;
		}
		this->ionex_tec[n++] = this->ionex_tec[i];
	}
	this->nt = n;

}

/* read ionex tec grid file ----------------------------------------------------
* read ionex ionospheric tec grid file
* args   : char   *file       I   ionex tec grid file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*                                 nav->nt, nav->ntmax and nav->tec are modified
*          int    opt         I   read option (1: no clear of tec data,0:clear)
* return : none
* notes  : see ref [1]
*-----------------------------------------------------------------------------*/
/* interpolate tec grid data -------------------------------------------------*/

//这里有点问题，这里获取的tec值都是第一个时段的TEC值，所以会有错误
int IonexModel::interptec(int k, const double *posp, double *value, double *rms)
{
	double dlat, dlon, a, b, d[4] = { 0 }, r[4] = { 0 };
	int i, j, n, index;

	*value = *rms = 0.0;

	if (this->ionex_tec->lats[2] == 0.0 || this->ionex_tec->lons[2] == 0.0)
		return 0;

	dlat = posp[0] * R2D - this->ionex_tec->lats[0];
	dlon = posp[1] * R2D - this->ionex_tec->lons[0];
	if (this->ionex_tec->lons[2]>0.0)
		dlon -= floor(dlon / 360)*360.0; /*  0<=dlon<360 */
	else
		dlon += floor(-dlon / 360)*360.0; /* -360<dlon<=0 */

	a = dlat / this->ionex_tec->lats[2];
	b = dlon / this->ionex_tec->lons[2];
	i = (int)floor(a); a -= i;
	j = (int)floor(b); b -= j;

	/* get gridded tec data */
	for (n = 0; n<4; n++)
	{
		if ((index = dataindex(i + (n % 2), j + (n<2 ? 0 : 1), k, this->ionex_tec->ndata))<0)
			continue;
		d[n] = this->ionex_tec->data[index];
		r[n] = this->ionex_tec->rms[index];
	}
	if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0)
	{

		/* bilinear interpolation (inside of grid) */
		*value = (1.0 - a)*(1.0 - b)*d[0] + a*(1.0 - b)*d[1] + (1.0 - a)*b*d[2] + a*b*d[3];
		*rms = (1.0 - a)*(1.0 - b)*r[0] + a*(1.0 - b)*r[1] + (1.0 - a)*b*r[2] + a*b*r[3];
	}
	/* nearest-neighbour extrapolation (outside of grid) */
	else if (a <= 0.5&&b <= 0.5&&d[0]>0.0) { *value = d[0]; *rms = r[0]; }
	else if (a> 0.5&&b <= 0.5&&d[1]>0.0) { *value = d[1]; *rms = r[1]; }
	else if (a <= 0.5&&b> 0.5&&d[2]>0.0) { *value = d[2]; *rms = r[2]; }
	else if (a> 0.5&&b> 0.5&&d[3]>0.0) { *value = d[3]; *rms = r[3]; }
	else {
		i = 0;
		for (n = 0; n<4; n++) if (d[n]>0.0) { i++; *value += d[n]; *rms += r[n]; }
		if (i == 0) return 0;
		*value /= i; *rms /= i;
	}
	return 1;
}

// add by limin
int IonexModel::interptec(int k, const double *posp, const ionex_tec_t *tec, double *value, double *rms)
{
	double dlat, dlon, a, b, d[4] = { 0 }, r[4] = { 0 };
	int i, j, n, index;

	*value = *rms = 0.0;

	if (tec->lats[2] == 0.0 || tec->lons[2] == 0.0)
		return 0;

	dlat = posp[0] * R2D - tec->lats[0];
	dlon = posp[1] * R2D - tec->lons[0];
	if (tec->lons[2]>0.0)
		dlon -= floor(dlon / 360)*360.0; /*  0<=dlon<360 */
	else
		dlon += floor(-dlon / 360)*360.0; /* -360<dlon<=0 */

	a = dlat / tec->lats[2];
	b = dlon / tec->lons[2];
	i = (int)floor(a); a -= i;
	j = (int)floor(b); b -= j;

	/* get gridded tec data */
	for (n = 0; n<4; n++)
	{
		if ((index = dataindex(i + (n % 2), j + (n<2 ? 0 : 1), k, tec->ndata))<0)
			continue;
		d[n] = tec->data[index];
		r[n] = tec->rms[index];
	}
	if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0)
	{

		/* bilinear interpolation (inside of grid) */
		*value = (1.0 - a)*(1.0 - b)*d[0] + a*(1.0 - b)*d[1] + (1.0 - a)*b*d[2] + a*b*d[3];
		*rms = (1.0 - a)*(1.0 - b)*r[0] + a*(1.0 - b)*r[1] + (1.0 - a)*b*r[2] + a*b*r[3];
	}
	/* nearest-neighbour extrapolation (outside of grid) */
	else if (a <= 0.5&&b <= 0.5&&d[0]>0.0) { *value = d[0]; *rms = r[0]; }
	else if (a> 0.5&&b <= 0.5&&d[1]>0.0) { *value = d[1]; *rms = r[1]; }
	else if (a <= 0.5&&b> 0.5&&d[2]>0.0) { *value = d[2]; *rms = r[2]; }
	else if (a> 0.5&&b> 0.5&&d[3]>0.0) { *value = d[3]; *rms = r[3]; }
	else {
		i = 0;
		for (n = 0; n<4; n++) if (d[n]>0.0) { i++; *value += d[n]; *rms += r[n]; }
		if (i == 0) return 0;
		*value /= i; *rms /= i;
	}
	return 1;
}

/* read ionex tec grid file ----------------------------------------------------
* read ionex ionospheric tec grid file
* args   : char   *file       I   ionex tec grid file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*                                 nav->nt, nav->ntmax and nav->tec are modified
*          int    opt         I   read option (1: no clear of tec data,0:clear)
* return : none
* notes  : see ref [1]
*-----------------------------------------------------------------------------*/
/* interpolate tec grid data -------------------------------------------------*/
int IonexModel::interptec(int k, const double *posp, int isession, double *value, double *rms)
{
	double dlat, dlon, a, b, d[4] = { 0 }, r[4] = { 0 };
	int i, j, n, index;

	*value = *rms = 0.0;

	if (this->ionex_tec->lats[2] == 0.0 || this->ionex_tec->lons[2] == 0.0)
		return 0;

	dlat = posp[0] * R2D - this->ionex_tec[isession].lats[0];
	dlon = posp[1] * R2D - this->ionex_tec[isession].lons[0];
	if (this->ionex_tec[isession].lons[2]>0.0)
		dlon -= floor(dlon / 360)*360.0; /*  0<=dlon<360 */
	else
		dlon += floor(-dlon / 360)*360.0; /* -360<dlon<=0 */

	a = dlat / this->ionex_tec[isession].lats[2];
	b = dlon / this->ionex_tec[isession].lons[2];
	i = (int)floor(a); a -= i;
	j = (int)floor(b); b -= j;

	/* get gridded tec data */
	for (n = 0; n<4; n++)
	{
		if ((index = dataindex(i + (n % 2), j + (n<2 ? 0 : 1), k, this->ionex_tec[isession].ndata))<0)
			continue;
		d[n] = this->ionex_tec[isession].data[index];
		r[n] = this->ionex_tec[isession].rms[index];
	}
	if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0)
	{

		/* bilinear interpolation (inside of grid) */
		*value = (1.0 - a)*(1.0 - b)*d[0] + a*(1.0 - b)*d[1] + (1.0 - a)*b*d[2] + a*b*d[3];
		*rms = (1.0 - a)*(1.0 - b)*r[0] + a*(1.0 - b)*r[1] + (1.0 - a)*b*r[2] + a*b*r[3];
	}
	/* nearest-neighbour extrapolation (outside of grid) */
	else if (a <= 0.5&&b <= 0.5&&d[0]>0.0) { *value = d[0]; *rms = r[0]; }
	else if (a> 0.5&&b <= 0.5&&d[1]>0.0) { *value = d[1]; *rms = r[1]; }
	else if (a <= 0.5&&b> 0.5&&d[2]>0.0) { *value = d[2]; *rms = r[2]; }
	else if (a> 0.5&&b> 0.5&&d[3]>0.0) { *value = d[3]; *rms = r[3]; }
	else {
		i = 0;
		for (n = 0; n<4; n++) if (d[n]>0.0) { i++; *value += d[n]; *rms += r[n]; }
		if (i == 0) return 0;
		*value /= i; *rms /= i;
	}
	return 1;
}


/* ionosphere delay by tec grid data -----------------------------------------*/
int IonexModel::ion_tec_one_map(gtime_t time, const ionex_tec_t *tec, double ipp_b, double ipp_l,
	double *tec_sum, double *var)
{
	double posp[3] = { 0 }, vtec, rms, hion;
	int i;


	*tec_sum = *var = 0.0;

	for (i = 0; i<tec->ndata[2]; i++)
	{ /* for a layer */

		hion = tec->hgts[0] + tec->hgts[2] * i;

		posp[0] = ipp_b;
		posp[1] = ipp_l;

		/* earth rotation correction (sun-fixed coordinate) */
		//posp[1]+=2.0*PI*this->pubfun.TimeDiff(time,tec->time)/86400.0;

		/* interpolate tec grid data */
		//if (!interptec(i,posp,&vtec,&rms)) return 0;
		if (!interptec(i, posp, tec, &vtec, &rms)) return 0;

		*tec_sum += vtec;
		*var += rms*rms;
	}

	*var = sqrt(*var);

	return 1;
}

/* ionosphere model by tec grid data -------------------------------------------
* compute ionospheric delay by tec grid data
* args   : gtime_t time     I   time (gpst)
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    opt       I   model option
*                                bit0: 0:earth-fixed,1:sun-fixed
*                                bit1: 0:single-layer,1:modified single-layer
*          double *delay    O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
* notes  : before calling the function, read tec grid data by calling readtec()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
*-----------------------------------------------------------------------------*/
int IonexModel::iontec(gtime_t time, double ipp_b, double ipp_l, double *tec, double *var)
{
	char subroutine[20];
	double dels[2], vars[2], a, tt;
	int i, stat[2];

	memset(subroutine, 0, sizeof(subroutine));
	strcpy(subroutine, "iontec");

	*tec = 0.0;
	*var = 0.0;

	for (i = 0; i<this->nt; i++)
	{
		if (timediff(this->ionex_tec[i].time, time)>0.0)
			break;
	}
	if (i == 0 || i >= this->nt)
	{
		memset(this->reason, 0, sizeof(this->reason));
		//this->pubfun.MSGOUT_STD(PROCESS_INF,subroutine,this->reason,WARNING_INF);
		return 0;
	}
	if ((tt = timediff(this->ionex_tec[i].time, this->ionex_tec[i - 1].time)) == 0.0)
	{
		memset(this->reason, 0, sizeof(this->reason));
		sprintf(this->reason, "%s", "tec grid time interval error");
		//this->pubfun.MSGOUT_STD(PROCESS_INF,subroutine,this->reason,WARNING_INF);
		return 0;
	}
	/* ionospheric delay by tec grid data */
	stat[0] = this->ion_tec_one_map(time, this->ionex_tec + i - 1, ipp_b, ipp_l, dels, vars);
	stat[1] = this->ion_tec_one_map(time, this->ionex_tec + i, ipp_b, ipp_l, dels + 1, vars + 1);

	if (!stat[0] && !stat[1])
	{

		memset(this->reason, 0, sizeof(this->reason));

		//this->pubfun.MSGOUT_STD(PROCESS_INF,subroutine,this->reason,WARNING_INF);
		return 0;
	}
	if (stat[0] && stat[1])
	{
		/* linear interpolation by time */
		a = timediff(time, this->ionex_tec[i - 1].time) / tt;
		*tec = dels[0] * (1.0 - a) + dels[1] * a;
		*var = vars[0] * (1.0 - a) + vars[1] * a;
	}
	else if (stat[0])
	{
		/* nearest-neighbor extrapolation by time */
		*tec = dels[0];
		*var = vars[0];
	}
	else
	{
		*tec = dels[1];
		*var = vars[1];
	}

	return 1;
}

int IonexModel::GetMapCount()
{
	return this->nt;
}

gtime_t IonexModel::GetMapTime(int imap)
{
	gtime_t tmp = { 0 };

	if (imap<0 || imap >= this->nt)
		return tmp;

	return this->ionex_tec[imap].time;
}

int IonexModel::GetGridCount(int imap)
{
	if (imap<0 || imap >= this->nt)
		return 0;
	return this->ionex_tec[imap].ndata[0] * this->ionex_tec[imap].ndata[1] * this->ionex_tec[imap].ndata[2];
}

double IonexModel::GetGridValue(int imap, int m, double &lat, double &lon, double &tec, double &rms)
{
	tec = 0.0;
	rms = 0.0;
	lat = 0.0;
	lon = 0.0;

	if (imap<0 || imap >= this->nt)
		return tec;

	if (m<0 || m >= this->ionex_tec[imap].ndata[0] * this->ionex_tec[imap].ndata[1] * this->ionex_tec[imap].ndata[2])
		return tec;

	tec = this->ionex_tec[imap].data[m];
	rms = this->ionex_tec[imap].rms[m];
	lat = this->ionex_tec[imap].grid_lat[m] * D2R;
	lon = this->ionex_tec[imap].grid_lon[m] * D2R;

	return tec;

}
