#include <ctype.h>

#include "SUPREME_Tide.h"
#include "SUPREME_CommonFunction.h"
#include "SUPREME_SunMoon.h"
#include "SUPREME_Coordinate.h"

/* read blq record -----------------------------------------------------------*/
int readblqrecord(FILE *fp, double *odisp)
{
	double v[11];
	char buff[256];
	int i, n = 0;

	while (fgets(buff, sizeof(buff), fp)) {
		if (!strncmp(buff, "$$", 2)) continue;
		if (sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			v, v + 1, v + 2, v + 3, v + 4, v + 5, v + 6, v + 7, v + 8, v + 9, v + 10)<11) continue;
		for (i = 0; i<11; i++) odisp[n + i * 6] = v[i];
		if (++n == 6) return 1;
	}
	return 0;
}

/// Read BLQ file
int ReadBlqFile(const char *file, const char *sta, double *odisp)
{
	FILE *fp;
	char buff[256], staname[32] = "", name[32], *p;

	/* station name to upper case */
	sscanf(sta, "%16s", staname);
	for (p = staname; (*p = (char)toupper((int)(*p))); p++);

	if (!(fp = fopen(file, "r")))
	{
		printf("Error: Blq file open failed\n", file);
		return 0;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		if (!strncmp(buff, "$$", 2) || strlen(buff)<2) continue;

		if (sscanf(buff + 2, "%16s", name)<1) continue;
		for (p = name; (*p = (char)toupper((int)(*p))); p++);
		if (strcmp(name, staname)) continue;

		/* read blq record */
		if (readblqrecord(fp, odisp)) {
			fclose(fp);
			return 1;
		}
	}
	fclose(fp);

	return 1;
}

/* read earth rotation parameters ----------------------------------------------
* read earth rotation parameters
* args   : char   *file       I   IGS ERP file (IGS ERP ver.2)
*          erp_t  *erp        O   earth rotation parameters
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
int ReadErpFile(const char *file, erp_t *erp)
{
	FILE *fp;
	erpd_t *erp_data;
	double v[14] = { 0 };
	char buff[256];

	if (!(fp = fopen(file, "r"))) {
		printf("erp file open error: file=%s\n", file);
		return 0;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		if (sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			v, v + 1, v + 2, v + 3, v + 4, v + 5, v + 6, v + 7, v + 8, v + 9, v + 10, v + 11, v + 12, v + 13) < 5) {
			continue;
		}
		if (erp->n >= erp->nmax) {
			erp->nmax = erp->nmax <= 0 ? 128 : erp->nmax * 2;
			erp_data = (erpd_t *)realloc(erp->data, sizeof(erpd_t)*erp->nmax);
			if (!erp_data) {
				free(erp->data); erp->data = NULL; erp->n = erp->nmax = 0;
				fclose(fp);
				return 0;
			}
			erp->data = erp_data;
		}
		erp->data[erp->n].mjd = v[0];
		erp->data[erp->n].xp = v[1] * 1E-6*AS2R;
		erp->data[erp->n].yp = v[2] * 1E-6*AS2R;
		erp->data[erp->n].ut1_utc = v[3] * 1E-7;
		erp->data[erp->n].lod = v[4] * 1E-7;
		erp->data[erp->n].xpr = v[12] * 1E-6*AS2R;
		erp->data[erp->n++].ypr = v[13] * 1E-6*AS2R;
	}
	fclose(fp);

	return 1;
}

/* solar/lunar tides (ref [2] 7) ---------------------------------------------*/
static void tide_pl(const double *eu, const double *rp, double GMp,
	const double *pos, double *dr)
{
	const double H3 = 0.292, L3 = 0.015;
	double r, ep[3], latp, lonp, p, K2, K3, a, H2, L2, dp, du, cosp, sinl, cosl;
	int i;

	if ((r = Norm(rp, 3)) <= 0.0) return;

	for (i = 0; i<3; i++) ep[i] = rp[i] / r;

	K2 = GMp / E_GM*SQR(WGS84_A)*SQR(WGS84_A) / (r*r*r);
	K3 = K2*WGS84_A / r;
	latp = asin(ep[2]); lonp = atan2(ep[1], ep[0]);
	cosp = cos(latp); sinl = sin(pos[0]); cosl = cos(pos[0]);

	/* step1 in phase (degree 2) */
	p = (3.0*sinl*sinl - 1.0) / 2.0;
	H2 = 0.6078 - 0.0006*p;
	L2 = 0.0847 + 0.0002*p;
	a = Dot(ep, eu, 3);
	dp = K2*3.0*L2*a;
	du = K2*(H2*(1.5*a*a - 0.5) - 3.0*L2*a*a);

	/* step1 in phase (degree 3) */
	dp += K3*L3*(7.5*a*a - 1.5);
	du += K3*(H3*(2.5*a*a*a - 1.5*a) - L3*(7.5*a*a - 1.5)*a);

	/* step1 out-of-phase (only radial) */
	du += 3.0 / 4.0*0.0025*K2*sin(2.0*latp)*sin(2.0*pos[0])*sin(pos[1] - lonp);
	du += 3.0 / 4.0*0.0022*K2*cosp*cosp*cosl*cosl*sin(2.0*(pos[1] - lonp));

	dr[0] = dp*ep[0] + du*eu[0];
	dr[1] = dp*ep[1] + du*eu[1];
	dr[2] = dp*ep[2] + du*eu[2];
}

/* Displacement by solid earth tide (ref [2] 7) ------------------------------*/
void SolidTideEffect(const double *rsun, const double *rmoon,
	const double *pos, const double *E, double gmst, int opt, double *dr)
{
	double dr1[3], dr2[3], eu[3], du, dn, sinl, sin2l;

	/* step1: time domain */
	eu[0] = E[2]; eu[1] = E[5]; eu[2] = E[8];
	tide_pl(eu, rsun, S_GM, pos, dr1);
	tide_pl(eu, rmoon, M_GM, pos, dr2);

	/* step2: frequency domain, only K1 radial */
	sin2l = sin(2.0*pos[0]);
	du = -0.012*sin2l*sin(gmst + pos[1]);

	dr[0] = dr1[0] + dr2[0] + du*E[2];
	dr[1] = dr1[1] + dr2[1] + du*E[5];
	dr[2] = dr1[2] + dr2[2] + du*E[8];

	/* eliminate permanent deformation */
	if (opt & 8) {
		sinl = sin(pos[0]);
		du = 0.1196*(1.5*sinl*sinl - 0.5);
		dn = 0.0247*sin2l;
		dr[0] += du*E[2] + dn*E[1];
		dr[1] += du*E[5] + dn*E[4];
		dr[2] += du*E[8] + dn*E[7];
	}
}

/* displacement by ocean tide loading (ref [2] 7) ----------------------------*/
void OceanTideEffect(gtime_t tut, const double *odisp, double *denu)
{
	const double args[][5] = {
		{ 1.40519E-4, 2.0, -2.0, 0.0, 0.00 },  /* M2 */
		{ 1.45444E-4, 0.0, 0.0, 0.0, 0.00 },  /* S2 */
		{ 1.37880E-4, 2.0, -3.0, 1.0, 0.00 },  /* N2 */
		{ 1.45842E-4, 2.0, 0.0, 0.0, 0.00 },  /* K2 */
		{ 0.72921E-4, 1.0, 0.0, 0.0, 0.25 },  /* K1 */
		{ 0.67598E-4, 1.0, -2.0, 0.0, -0.25 },  /* O1 */
		{ 0.72523E-4, -1.0, 0.0, 0.0, -0.25 },  /* P1 */
		{ 0.64959E-4, 1.0, -3.0, 1.0, -0.25 },  /* Q1 */
		{ 0.53234E-5, 0.0, 2.0, 0.0, 0.00 },  /* Mf */
		{ 0.26392E-5, 0.0, 1.0, -1.0, 0.00 },  /* Mm */
		{ 0.03982E-5, 2.0, 0.0, 0.0, 0.00 }   /* Ssa */
	};
	const double ep1975[] = { 1975, 1, 1, 0, 0, 0 };
	double ep[6], fday, days, t, t2, t3, a[5], ang, dp[3] = { 0 };
	int i, j;

	/* angular argument: see subroutine arg.f for reference [1] */
	time2epoch(tut, ep);
	fday = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];
	ep[3] = ep[4] = ep[5] = 0.0;
	days = timediff(epoch2time(ep), epoch2time(ep1975)) / 86400.0 + 1.0;
	t = (27392.500528 + 1.000000035*days) / 36525.0;
	t2 = t*t; t3 = t2*t;

	a[0] = fday;
	a[1] = (279.69668 + 36000.768930485*t + 3.03E-4*t2)*D2R; /* H0 */
	a[2] = (270.434358 + 481267.88314137*t - 0.001133*t2 + 1.9E-6*t3)*D2R; /* S0 */
	a[3] = (334.329653 + 4069.0340329577*t - 0.010325*t2 - 1.2E-5*t3)*D2R; /* P0 */
	a[4] = 2.0*PI;

	/* displacements by 11 constituents */
	for (i = 0; i<11; i++) {
		ang = 0.0;
		for (j = 0; j<5; j++) ang += a[j] * args[i][j];
		for (j = 0; j<3; j++) dp[j] += odisp[j + i * 6] * cos(ang - odisp[j + 3 + i * 6] * D2R);
	}
	denu[0] = -dp[1];
	denu[1] = -dp[2];
	denu[2] = dp[0];
}

/* iers mean pole (ref [7] eq.7.25) ------------------------------------------*/
static void iers_mean_pole(gtime_t tut, double *xp_bar, double *yp_bar)
{
	const double ep2000[] = { 2000, 1, 1, 0, 0, 0 };
	double y, y2, y3;

	y = timediff(tut, epoch2time(ep2000)) / 86400.0 / 365.25;

	if (y<3653.0 / 365.25) { /* until 2010.0 */
		y2 = y*y; y3 = y2*y;
		*xp_bar = 55.974 + 1.8243*y + 0.18413*y2 + 0.007024*y3; /* (mas) */
		*yp_bar = 346.346 + 1.7896*y - 0.10729*y2 - 0.000908*y3;
	}
	else { /* after 2010.0 */
		*xp_bar = 23.513 + 7.6141*y; /* (mas) */
		*yp_bar = 358.891 - 0.6287*y;
	}
}

/* get earth rotation parameter values -----------------------------------------
* get earth rotation parameter values
* args   : erp_t  *erp        I   earth rotation parameters
*          gtime_t time       I   time (gpst)
*          double *erpv       O   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
int GetErp(const erp_t *erp, gtime_t time, double *erpv)
{
	const double ep[] = { 2000, 1, 1, 12, 0, 0 };
	double mjd, day, a;
	int i = 0, j, k;

	if (erp->n <= 0) return 0;

	mjd = 51544.5 + (timediff(gpst2utc(time), epoch2time(ep))) / 86400.0;

	if (mjd <= erp->data[0].mjd) {
		day = mjd - erp->data[0].mjd;
		erpv[0] = erp->data[0].xp + erp->data[0].xpr*day;
		erpv[1] = erp->data[0].yp + erp->data[0].ypr*day;
		erpv[2] = erp->data[0].ut1_utc - erp->data[0].lod*day;
		erpv[3] = erp->data[0].lod;
		return 1;
	}
	if (mjd >= erp->data[erp->n - 1].mjd) {
		day = mjd - erp->data[erp->n - 1].mjd;
		erpv[0] = erp->data[erp->n - 1].xp + erp->data[erp->n - 1].xpr*day;
		erpv[1] = erp->data[erp->n - 1].yp + erp->data[erp->n - 1].ypr*day;
		erpv[2] = erp->data[erp->n - 1].ut1_utc - erp->data[erp->n - 1].lod*day;
		erpv[3] = erp->data[erp->n - 1].lod;
		return 1;
	}
	for (j = 0, k = erp->n - 1; j <= k;) {
		i = (j + k) / 2;
		if (mjd<erp->data[i].mjd) k = i - 1;
		else if (mjd>erp->data[i + 1].mjd) j = i + 1;
		else break;
	}
	if (erp->data[i].mjd == mjd - erp->data[i + 1].mjd) {
		a = 0.5;
	}
	else {
		a = (mjd - erp->data[i + 1].mjd) / (erp->data[i].mjd - mjd - erp->data[i + 1].mjd);

		if (i + 1 >= erp->n || i<0) {
			printf("i+1>=erp->n || i<0  %d\n", i);
			getchar();
		}
	}
	erpv[0] = (1.0 - a)*erp->data[i].xp + a*erp->data[i + 1].xp;
	erpv[1] = (1.0 - a)*erp->data[i].yp + a*erp->data[i + 1].yp;
	erpv[2] = (1.0 - a)*erp->data[i].ut1_utc + a*erp->data[i + 1].ut1_utc;
	erpv[3] = (1.0 - a)*erp->data[i].lod + a*erp->data[i + 1].lod;
	return 1;
}

/* displacement by pole tide (ref [7] eq.7.26) --------------------------------*/
void PoleTideEffect(gtime_t tut, const double *pos, const double *erpv, double *denu)
{
	double xp_bar, yp_bar, m1, m2, cosl, sinl;

	/* iers mean pole (mas) */
	iers_mean_pole(tut, &xp_bar, &yp_bar);

	/* ref [7] eq.7.24 */
	m1 = erpv[0] / AS2R - xp_bar*1E-3; /* (as) */
	m2 = -erpv[1] / AS2R + yp_bar*1E-3;

	/* sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi) */
	cosl = cos(pos[1]);
	sinl = sin(pos[1]);
	denu[0] = 9E-3*sin(pos[0])    *(m1*sinl - m2*cosl); /* de= Slambda (m) */
	denu[1] = -9E-3*cos(2.0*pos[0])*(m1*cosl + m2*sinl); /* dn=-Stheta  (m) */
	denu[2] = -33E-3*sin(2.0*pos[0])*(m1*cosl + m2*sinl); /* du= Sr      (m) */
}

/* Tidal displacement ----------------------------------------------------------
* displacements by earth tides
* args   : gtime_t tutc     I   time in utc
*          double *rr       I   site position (ecef) (m)
*          int    opt       I   options (or of the followings)
*                                 1: solid earth tide
*                                 2: ocean tide loading
*                                 4: pole tide
*                                 8: elimate permanent deformation
*          double *erp      I   earth rotation parameters (NULL: not used)
*          double *odisp    I   ocean loading parameters  (NULL: not used)
*                                 odisp[0+i*6]: consituent i amplitude radial(m)
*                                 odisp[1+i*6]: consituent i amplitude west  (m)
*                                 odisp[2+i*6]: consituent i amplitude south (m)
*                                 odisp[3+i*6]: consituent i phase radial  (deg)
*                                 odisp[4+i*6]: consituent i phase west    (deg)
*                                 odisp[5+i*6]: consituent i phase south   (deg)
*                                (i=0:M2,1:S2,2:N2,3:K2,4:K1,5:O1,6:P1,7:Q1,
*                                   8:Mf,9:Mm,10:Ssa)
*          double *dr       O   displacement by earth tides (ecef) (m)
* return : none
* notes  : see ref [1], [2] chap 7
*          see ref [4] 5.2.1, 5.2.2, 5.2.3
*          ver.2.4.0 does not use ocean loading and pole tide corrections
*-----------------------------------------------------------------------------*/
int TideDisp(gnsstime &tutc, const double *rr, int opt, const erp_t *erp, const double *odisp, double *dr)
{
	gtime_t tut;
	double pos[2], E[9], drt[3], denu[3], rs[3], rm[3], gmst, erpv[5] = { 0 };
	int i;

	if (erp)
	{
		// Get erp data
		GetErp(erp, utc2gpst(tutc.m_gtime), erpv);
	}
	tut = timeadd(tutc.m_gtime, erpv[2]);

	dr[0] = dr[1] = dr[2] = 0.0;

	if (Norm(rr, 3) <= 0.0) return 0;

	pos[0] = asin(rr[2] / Norm(rr, 3));
	pos[1] = atan2(rr[1], rr[0]);
	xyz2enu(pos, E);

	// Solid time effect
	GetSunMoonPos(tutc, rs, rm, &gmst);
	SolidTideEffect(rs, rm, pos, E, gmst, opt, drt);
	for (i = 0; i < 3; i++) dr[i] += drt[i];

	// Ocean tide loading
	if (odisp)
	{
		OceanTideEffect(tut, odisp, denu);
		matmul("TN", 3, 1, 3, 1.0, E, denu, 0.0, drt);
		for (i = 0; i<3; i++) dr[i] += drt[i];
	}

	// Pole tide
	if (erp)
	{
		PoleTideEffect(tut, pos, erpv, denu);
		matmul("TN", 3, 1, 3, 1.0, E, denu, 0.0, drt);
		for (i = 0; i<3; i++) dr[i] += drt[i];
	}

	return 1;
}
