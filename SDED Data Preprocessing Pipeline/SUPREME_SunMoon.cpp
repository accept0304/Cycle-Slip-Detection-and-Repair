#include "SUPREME_SunMoon.h"
#include "SUPREME_CommonFunction.h"
#include "SUPREME_Coordinate.h"

/* astronomical arguments: f={l,l',F,D,OMG} (rad) ----------------------------*/
static void ast_args(double t, double *f)
{
	static const double fc[][5] = { /* coefficients for iau 1980 nutation */
		{ 134.96340251, 1717915923.2178, 31.8792, 0.051635, -0.00024470 },
		{ 357.52910918, 129596581.0481, -0.5532, 0.000136, -0.00001149 },
		{ 93.27209062, 1739527262.8478, -12.7512, -0.001037, 0.00000417 },
		{ 297.85019547, 1602961601.2090, -6.3706, 0.006593, -0.00003169 },
		{ 125.04455501, -6962890.2665, 7.4722, 0.007702, -0.00005939 }
	};
	double tt[4];
	int i, j;

	for (tt[0] = t, i = 1; i<4; i++) tt[i] = tt[i - 1] * t;
	for (i = 0; i<5; i++) {
		f[i] = fc[i][0] * 3600.0;
		for (j = 0; j<4; j++) f[i] += fc[i][j + 1] * tt[j];
		f[i] = fmod(f[i] * AS2R, 2.0*PI);
	}
}

/* sun and moon position in eci --------------------------*/
void sunmoonpos_eci(gtime_t tut, double *rsun, double *rmoon)
{
	const double ep2000[] = { 2000, 1, 1, 12, 0, 0 };
	double t, f[5], eps, Ms, ls, rs, lm, pm, rm, sine, cose, sinp, cosp, sinl, cosl;

	t = timediff(tut, epoch2time(ep2000)) / 86400.0 / 36525.0;

	/* astronomical arguments */
	ast_args(t, f);

	/* obliquity of the ecliptic */
	eps = 23.439291 - 0.0130042*t;
	sine = sin(eps*D2R); cose = cos(eps*D2R);

	/* sun position in eci */
	if (rsun) {
		Ms = 357.5277233 + 35999.05034*t;
		ls = 280.460 + 36000.770*t + 1.914666471*sin(Ms*D2R) + 0.019994643*sin(2.0*Ms*D2R);
		rs = AU*(1.000140612 - 0.016708617*cos(Ms*D2R) - 0.000139589*cos(2.0*Ms*D2R));
		sinl = sin(ls*D2R); cosl = cos(ls*D2R);
		rsun[0] = rs*cosl;
		rsun[1] = rs*cose*sinl;
		rsun[2] = rs*sine*sinl;
	}
	/* moon position in eci */
	if (rmoon) {
		lm = 218.32 + 481267.883*t + 6.29*sin(f[0]) - 1.27*sin(f[0] - 2.0*f[3]) +
			0.66*sin(2.0*f[3]) + 0.21*sin(2.0*f[0]) - 0.19*sin(f[1]) - 0.11*sin(2.0*f[2]);
		pm = 5.13*sin(f[2]) + 0.28*sin(f[0] + f[2]) - 0.28*sin(f[2] - f[0]) -
			0.17*sin(f[2] - 2.0*f[3]);
		rm = WGS84_A / sin((0.9508 + 0.0518*cos(f[0]) + 0.0095*cos(f[0] - 2.0*f[3]) +
			0.0078*cos(2.0*f[3]) + 0.0028*cos(2.0*f[0]))*D2R);
		sinl = sin(lm*D2R); cosl = cos(lm*D2R);
		sinp = sin(pm*D2R); cosp = cos(pm*D2R);
		rmoon[0] = rm*cosp*cosl;
		rmoon[1] = rm*(cose*cosp*sinl - sine*sinp);
		rmoon[2] = rm*(sine*cosp*sinl + cose*sinp);
	}
}


/* Get sun moon position by program-self -------------------------------------------------------------
* Parameter:
*            gtime_t      tutc          I                time in ut1
*            double       *erpv         I                erp value {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
*            double       *rsun         O                sun position in ecef
*            double       *rmoon        O                moon position in ecef
*            double       *gmst         O                gmst (rad)
* --------------------------------------------------------------------------------------------------- */
void GetSunMoonPos(gnsstime &obst, double *rsun, double *rmoon, double *gmst)
{
	double erpv[5] = { 0 };
	if (obst.m_Type != UTC)
		obst.Gpst2Utc();

	gtime_t tut;
	double rs[3], rm[3], U[9], gmst_;

	tut = timeadd(obst.m_gtime, erpv[2]); /* utc -> ut1 */

	/* sun and moon position in eci */
	sunmoonpos_eci(tut, rsun ? rs : NULL, rmoon ? rm : NULL);

	/* eci to ecef transformation MatrixT */
	eci2ecef(obst.m_gtime, erpv, U, &gmst_);

	/* sun and moon postion in ecef */
	if (rsun) matmul("NN", 3, 1, 3, 1.0, U, rs, 0.0, rsun);
	if (rmoon) matmul("NN", 3, 1, 3, 1.0, U, rm, 0.0, rmoon);
	if (gmst) *gmst = gmst_;
}
