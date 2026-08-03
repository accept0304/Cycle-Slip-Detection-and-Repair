#include "SUPREME_Matrix.h"
#include "SUPREME_Coordinate.h"

/* coordinate rotation MatrixT ------------------------------------------------*/
#define Rx(t,X) do { \
	(X)[0] = 1.0; (X)[1] = (X)[2] = (X)[3] = (X)[6] = 0.0; \
	(X)[4] = (X)[8] = cos(t); (X)[7] = sin(t); (X)[5] = -(X)[7]; \
} while (0)

#define Ry(t,X) do { \
	(X)[4] = 1.0; (X)[1] = (X)[3] = (X)[5] = (X)[7] = 0.0; \
	(X)[0] = (X)[8] = cos(t); (X)[2] = sin(t); (X)[6] = -(X)[2]; \
} while (0)

#define Rz(t,X) do { \
	(X)[8] = 1.0; (X)[2] = (X)[5] = (X)[6] = (X)[7] = 0.0; \
	(X)[0] = (X)[4] = cos(t); (X)[3] = sin(t); (X)[1] = -(X)[3]; \
} while (0)

/* iau 1980 nutation ---------------------------------------------------------*/
static void nut_iau1980(double t, const double *f, double *dpsi, double *deps)
{
	static const double nut[106][10] =
	{
		{ 0, 0, 0, 0, 1, -6798.4, -171996, -174.2, 92025, 8.9 },
		{ 0, 0, 2, -2, 2, 182.6, -13187, -1.6, 5736, -3.1 },
		{ 0, 0, 2, 0, 2, 13.7, -2274, -0.2, 977, -0.5 },
		{ 0, 0, 0, 0, 2, -3399.2, 2062, 0.2, -895, 0.5 },
		{ 0, -1, 0, 0, 0, -365.3, -1426, 3.4, 54, -0.1 },
		{ 1, 0, 0, 0, 0, 27.6, 712, 0.1, -7, 0.0 },
		{ 0, 1, 2, -2, 2, 121.7, -517, 1.2, 224, -0.6 },
		{ 0, 0, 2, 0, 1, 13.6, -386, -0.4, 200, 0.0 },
		{ 1, 0, 2, 0, 2, 9.1, -301, 0.0, 129, -0.1 },
		{ 0, -1, 2, -2, 2, 365.2, 217, -0.5, -95, 0.3 },
		{ -1, 0, 0, 2, 0, 31.8, 158, 0.0, -1, 0.0 },
		{ 0, 0, 2, -2, 1, 177.8, 129, 0.1, -70, 0.0 },
		{ -1, 0, 2, 0, 2, 27.1, 123, 0.0, -53, 0.0 },
		{ 1, 0, 0, 0, 1, 27.7, 63, 0.1, -33, 0.0 },
		{ 0, 0, 0, 2, 0, 14.8, 63, 0.0, -2, 0.0 },
		{ -1, 0, 2, 2, 2, 9.6, -59, 0.0, 26, 0.0 },
		{ -1, 0, 0, 0, 1, -27.4, -58, -0.1, 32, 0.0 },
		{ 1, 0, 2, 0, 1, 9.1, -51, 0.0, 27, 0.0 },
		{ -2, 0, 0, 2, 0, -205.9, -48, 0.0, 1, 0.0 },
		{ -2, 0, 2, 0, 1, 1305.5, 46, 0.0, -24, 0.0 },
		{ 0, 0, 2, 2, 2, 7.1, -38, 0.0, 16, 0.0 },
		{ 2, 0, 2, 0, 2, 6.9, -31, 0.0, 13, 0.0 },
		{ 2, 0, 0, 0, 0, 13.8, 29, 0.0, -1, 0.0 },
		{ 1, 0, 2, -2, 2, 23.9, 29, 0.0, -12, 0.0 },
		{ 0, 0, 2, 0, 0, 13.6, 26, 0.0, -1, 0.0 },
		{ 0, 0, 2, -2, 0, 173.3, -22, 0.0, 0, 0.0 },
		{ -1, 0, 2, 0, 1, 27.0, 21, 0.0, -10, 0.0 },
		{ 0, 2, 0, 0, 0, 182.6, 17, -0.1, 0, 0.0 },
		{ 0, 2, 2, -2, 2, 91.3, -16, 0.1, 7, 0.0 },
		{ -1, 0, 0, 2, 1, 32.0, 16, 0.0, -8, 0.0 },
		{ 0, 1, 0, 0, 1, 386.0, -15, 0.0, 9, 0.0 },
		{ 1, 0, 0, -2, 1, -31.7, -13, 0.0, 7, 0.0 },
		{ 0, -1, 0, 0, 1, -346.6, -12, 0.0, 6, 0.0 },
		{ 2, 0, -2, 0, 0, -1095.2, 11, 0.0, 0, 0.0 },
		{ -1, 0, 2, 2, 1, 9.5, -10, 0.0, 5, 0.0 },
		{ 1, 0, 2, 2, 2, 5.6, -8, 0.0, 3, 0.0 },
		{ 0, -1, 2, 0, 2, 14.2, -7, 0.0, 3, 0.0 },
		{ 0, 0, 2, 2, 1, 7.1, -7, 0.0, 3, 0.0 },
		{ 1, 1, 0, -2, 0, -34.8, -7, 0.0, 0, 0.0 },
		{ 0, 1, 2, 0, 2, 13.2, 7, 0.0, -3, 0.0 },
		{ -2, 0, 0, 2, 1, -199.8, -6, 0.0, 3, 0.0 },
		{ 0, 0, 0, 2, 1, 14.8, -6, 0.0, 3, 0.0 },
		{ 2, 0, 2, -2, 2, 12.8, 6, 0.0, -3, 0.0 },
		{ 1, 0, 0, 2, 0, 9.6, 6, 0.0, 0, 0.0 },
		{ 1, 0, 2, -2, 1, 23.9, 6, 0.0, -3, 0.0 },
		{ 0, 0, 0, -2, 1, -14.7, -5, 0.0, 3, 0.0 },
		{ 0, -1, 2, -2, 1, 346.6, -5, 0.0, 3, 0.0 },
		{ 2, 0, 2, 0, 1, 6.9, -5, 0.0, 3, 0.0 },
		{ 1, -1, 0, 0, 0, 29.8, 5, 0.0, 0, 0.0 },
		{ 1, 0, 0, -1, 0, 411.8, -4, 0.0, 0, 0.0 },
		{ 0, 0, 0, 1, 0, 29.5, -4, 0.0, 0, 0.0 },
		{ 0, 1, 0, -2, 0, -15.4, -4, 0.0, 0, 0.0 },
		{ 1, 0, -2, 0, 0, -26.9, 4, 0.0, 0, 0.0 },
		{ 2, 0, 0, -2, 1, 212.3, 4, 0.0, -2, 0.0 },
		{ 0, 1, 2, -2, 1, 119.6, 4, 0.0, -2, 0.0 },
		{ 1, 1, 0, 0, 0, 25.6, -3, 0.0, 0, 0.0 },
		{ 1, -1, 0, -1, 0, -3232.9, -3, 0.0, 0, 0.0 },
		{ -1, -1, 2, 2, 2, 9.8, -3, 0.0, 1, 0.0 },
		{ 0, -1, 2, 2, 2, 7.2, -3, 0.0, 1, 0.0 },
		{ 1, -1, 2, 0, 2, 9.4, -3, 0.0, 1, 0.0 },
		{ 3, 0, 2, 0, 2, 5.5, -3, 0.0, 1, 0.0 },
		{ -2, 0, 2, 0, 2, 1615.7, -3, 0.0, 1, 0.0 },
		{ 1, 0, 2, 0, 0, 9.1, 3, 0.0, 0, 0.0 },
		{ -1, 0, 2, 4, 2, 5.8, -2, 0.0, 1, 0.0 },
		{ 1, 0, 0, 0, 2, 27.8, -2, 0.0, 1, 0.0 },
		{ -1, 0, 2, -2, 1, -32.6, -2, 0.0, 1, 0.0 },
		{ 0, -2, 2, -2, 1, 6786.3, -2, 0.0, 1, 0.0 },
		{ -2, 0, 0, 0, 1, -13.7, -2, 0.0, 1, 0.0 },
		{ 2, 0, 0, 0, 1, 13.8, 2, 0.0, -1, 0.0 },
		{ 3, 0, 0, 0, 0, 9.2, 2, 0.0, 0, 0.0 },
		{ 1, 1, 2, 0, 2, 8.9, 2, 0.0, -1, 0.0 },
		{ 0, 0, 2, 1, 2, 9.3, 2, 0.0, -1, 0.0 },
		{ 1, 0, 0, 2, 1, 9.6, -1, 0.0, 0, 0.0 },
		{ 1, 0, 2, 2, 1, 5.6, -1, 0.0, 1, 0.0 },
		{ 1, 1, 0, -2, 1, -34.7, -1, 0.0, 0, 0.0 },
		{ 0, 1, 0, 2, 0, 14.2, -1, 0.0, 0, 0.0 },
		{ 0, 1, 2, -2, 0, 117.5, -1, 0.0, 0, 0.0 },
		{ 0, 1, -2, 2, 0, -329.8, -1, 0.0, 0, 0.0 },
		{ 1, 0, -2, 2, 0, 23.8, -1, 0.0, 0, 0.0 },
		{ 1, 0, -2, -2, 0, -9.5, -1, 0.0, 0, 0.0 },
		{ 1, 0, 2, -2, 0, 32.8, -1, 0.0, 0, 0.0 },
		{ 1, 0, 0, -4, 0, -10.1, -1, 0.0, 0, 0.0 },
		{ 2, 0, 0, -4, 0, -15.9, -1, 0.0, 0, 0.0 },
		{ 0, 0, 2, 4, 2, 4.8, -1, 0.0, 0, 0.0 },
		{ 0, 0, 2, -1, 2, 25.4, -1, 0.0, 0, 0.0 },
		{ -2, 0, 2, 4, 2, 7.3, -1, 0.0, 1, 0.0 },
		{ 2, 0, 2, 2, 2, 4.7, -1, 0.0, 0, 0.0 },
		{ 0, -1, 2, 0, 1, 14.2, -1, 0.0, 0, 0.0 },
		{ 0, 0, -2, 0, 1, -13.6, -1, 0.0, 0, 0.0 },
		{ 0, 0, 4, -2, 2, 12.7, 1, 0.0, 0, 0.0 },
		{ 0, 1, 0, 0, 2, 409.2, 1, 0.0, 0, 0.0 },
		{ 1, 1, 2, -2, 2, 22.5, 1, 0.0, -1, 0.0 },
		{ 3, 0, 2, -2, 2, 8.7, 1, 0.0, 0, 0.0 },
		{ -2, 0, 2, 2, 2, 14.6, 1, 0.0, -1, 0.0 },
		{ -1, 0, 0, 0, 2, -27.3, 1, 0.0, -1, 0.0 },
		{ 0, 0, -2, 2, 1, -169.0, 1, 0.0, 0, 0.0 },
		{ 0, 1, 2, 0, 1, 13.1, 1, 0.0, 0, 0.0 },
		{ -1, 0, 4, 0, 2, 9.1, 1, 0.0, 0, 0.0 },
		{ 2, 1, 0, -2, 0, 131.7, 1, 0.0, 0, 0.0 },
		{ 2, 0, 0, 2, 0, 7.1, 1, 0.0, 0, 0.0 },
		{ 2, 0, 2, -2, 1, 12.8, 1, 0.0, -1, 0.0 },
		{ 2, 0, -2, 0, 1, -943.2, 1, 0.0, 0, 0.0 },
		{ 1, -1, 0, -2, 0, -29.3, 1, 0.0, 0, 0.0 },
		{ -1, 0, 0, 1, 1, -388.3, 1, 0.0, 0, 0.0 },
		{ -1, -1, 0, 2, 1, 35.0, 1, 0.0, 0, 0.0 },
		{ 0, 1, 0, 1, 0, 27.3, 1, 0.0, 0, 0.0 }
	};
	double ang;
	int i, j;

	*dpsi = *deps = 0.0;

	for (i = 0; i<106; i++) {
		ang = 0.0;
		for (j = 0; j<5; j++) ang += nut[i][j] * f[j];
		*dpsi += (nut[i][6] + nut[i][7] * t)*sin(ang);
		*deps += (nut[i][8] + nut[i][9] * t)*cos(ang);
	}
	*dpsi *= 1E-4*AS2R; /* 0.1 mas -> rad */
	*deps *= 1E-4*AS2R;
}

/// ecef to position
void ecef2pos(const double *r, double *pos)
{
	double e2 = WGS84_F*(2.0 - WGS84_F), r2 = Dot(r, r, 2), z, zk, v = WGS84_A, sinp;

	for (z = r[2], zk = 0.0; fabs(z - zk) >= 1E-4;) {
		zk = z;
		sinp = z / sqrt(r2 + z*z);
		v = WGS84_A / sqrt(1.0 - e2*sinp*sinp);
		z = r[2] + v*e2*sinp;
	}
	pos[0] = r2>1E-12 ? atan(z / sqrt(r2)) : (r[2]>0.0 ? PI / 2.0 : -PI / 2.0);
	pos[1] = r2>1E-12 ? atan2(r[1], r[0]) : 0.0;
	pos[2] = sqrt(r2 + z*z) - v;
}

/// xyz to neu
void xyz2enu(double *pos, double *E)
{
	double sinp = sin(pos[0]), cosp = cos(pos[0]), sinl = sin(pos[1]), cosl = cos(pos[1]);

	E[0] = -sinl;      E[3] = cosl;       E[6] = 0.0;
	E[1] = -sinp*cosl; E[4] = -sinp*sinl; E[7] = cosp;
	E[2] = cosp*cosl;  E[5] = cosp*sinl;  E[8] = sinp;
}

/* transform ecef vector to local tangental coordinate -------------------------
* transform ecef vector to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *r        I   vector in ecef coordinate {x,y,z}
*          double *e        O   vector in local tangental coordinate {e,n,u}
* return : none
*-----------------------------------------------------------------------------*/
void ecef2enu(double *pos, double *r, double *e)
{
	double E[9];

	xyz2enu(pos, E);
	matmul("NN", 3, 1, 3, 1.0, E, r, 0.0, e);
}

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

/* eci to ecef transformation MatrixT -------------------------------------------
* compute eci to ecef transformation MatrixT
* args   : gtime_t tutc     I   time in utc
*          double *erpv     I   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
*          double *U        O   eci to ecef transformation MatrixT (3 x 3)
*          double *gmst     IO  greenwich mean sidereal time (rad)
*                               (NULL: no output)
* return : none
* note   : see ref [3] chap 5
*          not thread-safe
*-----------------------------------------------------------------------------*/
void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst)
{
	const double ep2000[] = { 2000, 1, 1, 12, 0, 0 };
	static gtime_t tutc_;
	static double U_[9], gmst_;
	gtime_t tgps;
	double eps, ze, th, z, t, t2, t3, dpsi, deps, gast, f[5];
	double R1[9], R2[9], R3[9], R[9], W[9], N[9], P[9], NP[9];
	int i;

	if (fabs(timediff(tutc, tutc_))<0.01) { /* read cache */
		for (i = 0; i<9; i++) U[i] = U_[i];
		if (gmst) *gmst = gmst_;
		return;
	}
	tutc_ = tutc;

	/* terrestrial time */
	tgps = utc2gpst(tutc_);
	t = (timediff(tgps, epoch2time(ep2000)) + 19.0 + 32.184) / 86400.0 / 36525.0;
	t2 = t*t; t3 = t2*t;

	/* astronomical arguments */
	ast_args(t, f);

	/* iau 1976 precession */
	ze = (2306.2181*t + 0.30188*t2 + 0.017998*t3)*AS2R;
	th = (2004.3109*t - 0.42665*t2 - 0.041833*t3)*AS2R;
	z = (2306.2181*t + 1.09468*t2 + 0.018203*t3)*AS2R;
	eps = (84381.448 - 46.8150*t - 0.00059*t2 + 0.001813*t3)*AS2R;
	Rz(-z, R1); Ry(th, R2); Rz(-ze, R3);
	matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, R);
	matmul("NN", 3, 3, 3, 1.0, R, R3, 0.0, P); /* P=Rz(-z)*Ry(th)*Rz(-ze) */

	/* iau 1980 nutation */
	nut_iau1980(t, f, &dpsi, &deps);
	Rx(-eps - deps, R1); Rz(-dpsi, R2); Rx(eps, R3);
	matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, R);
	matmul("NN", 3, 3, 3, 1.0, R, R3, 0.0, N); /* N=Rx(-eps)*Rz(-dspi)*Rx(eps) */

	/* greenwich aparent sidereal time (rad) */
	gmst_ = utc2gmst(tutc_, erpv[2]);
	gast = gmst_ + dpsi*cos(eps);
	gast += (0.00264*sin(f[4]) + 0.000063*sin(2.0*f[4]))*AS2R;

	/* eci to ecef transformation MatrixT */
	Ry(-erpv[0], R1); Rx(-erpv[1], R2); Rz(gast, R3);
	matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, W);
	matmul("NN", 3, 3, 3, 1.0, W, R3, 0.0, R); /* W=Ry(-xp)*Rx(-yp) */
	matmul("NN", 3, 3, 3, 1.0, N, P, 0.0, NP);
	matmul("NN", 3, 3, 3, 1.0, R, NP, 0.0, U_); /* U=W*Rz(gast)*N*P */

	for (i = 0; i<9; i++) U[i] = U_[i];
	if (gmst) *gmst = gmst_;
}

void Coordinate::_xyz2blh()
{
	double a = WGS84_A, e = WGS84_E;
	double pi = 4 * atan(1.0);
	double rr = sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y);

	if (rr<1.0e-15)
	{
		BLH.setB(pi / 2.0);
		BLH.setL(0.0);
		BLH.setH(XYZ._Z);
		return;
	}
	else
	{
		BLH.setB(atan2(XYZ._Z, rr));
	}
	if (fabs(XYZ._X)<1.0e-8)
	{
		BLH.setL(0.0);
	}
	else
	{
		BLH.setL(atan2(XYZ._Y, XYZ._X));
	}

	double temp1 = 0.0, temp2 = 0.0, N = 0.0;
	do
	{
		N = a / (sqrt(1 - e*sin(BLH._B)*sin(BLH._B)));
		temp1 = atan2(XYZ._Z + N*e*sin(BLH._B), sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y));
		temp2 = temp1 - BLH._B;
		BLH.setB(temp1);
	} while (fabs(temp2)>1.0e-10);

	BLH.setH(sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y) / cos(BLH._B) - N);

	if (fabs(fabs(BLH._B) - pi / 2.0) <= 1.0e-9)
	{
		double dz = 0.0;
		BLH.setB(1.535);
		do
		{
			N = a / (sqrt(1 - e*sin(BLH._B)*sin(BLH._B)));
			dz = e*N*sin(BLH._B);
			temp1 = atan2(XYZ._Z + dz, sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y));
			temp2 = temp1 - BLH._B;
			BLH.setB(temp1);
		} while (fabs(temp2)>1.0e-10);

		BLH.setH(sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y + (XYZ._Z + dz)*(XYZ._Z + dz)) - N);
	}
}
void Coordinate::_xyz2blh(double a, double e){
	double pi = 4 * atan(1.0);
	double rr = sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y);

	if (rr<1.0e-15){
		BLH.setB(pi / 2.0);
		BLH.setL(0.0);
		BLH.setH(XYZ._Z);
		return;
	}
	else{
		BLH.setB(atan2(XYZ._Z, rr));
	}
	if (fabs(XYZ._X)<1.0e-8)	{
		BLH.setL(0.0);
	}
	else{
		BLH.setL(atan2(XYZ._Y, XYZ._X));
	}

	double temp1 = 0.0, temp2 = 0.0, N = 0.0;
	do{
		N = a / (sqrt(1 - e*sin(BLH._B)*sin(BLH._B)));
		temp1 = atan2(XYZ._Z + N*e*sin(BLH._B), sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y));
		temp2 = temp1 - BLH._B;
		BLH.setB(temp1);
	} while (fabs(temp2)>1.0e-10);

	BLH.setH(sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y) / cos(BLH._B) - N);

	if (fabs(fabs(BLH._B) - pi / 2.0) <= 1.0e-9)	{
		double dz = 0.0;
		BLH.setB(1.535);
		do{
			N = a / (sqrt(1 - e*sin(BLH._B)*sin(BLH._B)));
			dz = e*N*sin(BLH._B);
			temp1 = atan2(XYZ._Z + dz, sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y));
			temp2 = temp1 - BLH._B;
			BLH.setB(temp1);
		} while (fabs(temp2)>1.0e-10);

		BLH.setH(sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y + (XYZ._Z + dz)*(XYZ._Z + dz)) - N);
	}
}

/* Geodetic Coordinate to Cartesian Coordinate */
void Coordinate::_blh2xyz(){
	double a = WGS84_A, e = WGS84_E;
	double n = a / (sqrt(1 - e*sin(BLH._B)*sin(BLH._B)));
	XYZ.setX((n + BLH._H)*cos(BLH._B)*cos(BLH._L));
	XYZ.setY((n + BLH._H)*cos(BLH._B)*sin(BLH._L));
	XYZ.setZ((n*(1 - e) + BLH._H)*sin(BLH._B));
}

void Coordinate::_blh2xyz(double a, double e){
	double n = a / (sqrt(1 - e*sin(BLH._B)*sin(BLH._B)));
	XYZ.setX((n + BLH._H)*cos(BLH._B)*cos(BLH._L));
	XYZ.setY((n + BLH._H)*cos(BLH._B)*sin(BLH._L));
	XYZ.setZ((n*(1 - e) + BLH._H)*sin(BLH._B));
}

/* Cartesian Coordinate to Topcentric Cooridnate */
void Coordinate::_xyz2neu(Cart_Crd& StaCent)
{
	double a = WGS84_A, e = WGS84_E;
	double tmp_x = 0, tmp_y = 0, tmp_z = 0;
	_xyz2blh();
	tmp_x = XYZ._X - StaCent._X;
	tmp_y = XYZ._Y - StaCent._Y;
	tmp_z = XYZ._Z - StaCent._Z;

	NEU.setN(-1 * sin(BLH._B)*cos(BLH._L)*tmp_x - sin(BLH._B)*sin(BLH._L)*tmp_y + cos(BLH._B)*tmp_z);
	NEU.setE(-1 * sin(BLH._L)*tmp_x + cos(BLH._L)*tmp_y);
	NEU.setU(cos(BLH._B)*cos(BLH._L)*tmp_x + cos(BLH._B)*sin(BLH._L)*tmp_y + sin(BLH._B)*tmp_z);

	/*double E[9] = { 0.0 };
	double sinp = sin(BLH._B), cosp = cos(BLH._B), sinl = sin(BLH._H), cosl = cos(BLH._H);
	E[0] = -sinl;      E[3] = cosl;       E[6] = 0.0;
	E[1] = -sinp*cosl; E[4] = -sinp*sinl; E[7] = cosp;
	E[2] = cosp*cosl;  E[5] = cosp*sinl;  E[8] = sinp;*/
}
void Coordinate::_xyz2neu(double a, double e, Cart_Crd& StaCent)
{
	double tmp_x = 0, tmp_y = 0, tmp_z = 0;

	Coordinate blh(StaCent);
	blh._xyz2blh();
	tmp_x = XYZ._X - StaCent._X;
	tmp_y = XYZ._Y - StaCent._Y;
	tmp_z = XYZ._Z - StaCent._Z;

	NEU.setN(-1 * sin(blh.BLH._B)*cos(blh.BLH._L)*tmp_x - sin(blh.BLH._B)*sin(blh.BLH._L)*tmp_y + cos(blh.BLH._B)*tmp_z);
	NEU.setE(-1 * sin(blh.BLH._L)*tmp_x + cos(blh.BLH._L)*tmp_y);
	NEU.setU(cos(blh.BLH._B)*cos(blh.BLH._L)*tmp_x + cos(blh.BLH._B)*sin(blh.BLH._L)*tmp_y + sin(blh.BLH._B)*tmp_z);
}

void Coordinate::_neu2xyz(Cart_Crd& StaCent)
{
	Coordinate tmp(StaCent);
	Geo_Crd blh;
	tmp._xyz2blh();
	blh = tmp.BLH;
	double B = blh._B, L = blh._L, H = blh._H;

	MatrixT R(3, 3), m_neu(3, 1), m_xyz(3, 1);
	R(1, 1) = -sin(B)*cos(L);
	R(1, 2) = -sin(B)*sin(L);
	R(1, 3) = cos(B);
	R(2, 1) = -sin(L);
	R(2, 2) = cos(L);
	R(2, 3) = 0;
	R(3, 1) = cos(B)*cos(L);
	R(3, 2) = cos(B)*sin(L);
	R(3, 3) = sin(B);

	m_neu(1, 1) = NEU._N;
	m_neu(2, 1) = NEU._E;
	m_neu(3, 1) = NEU._U;

	m_xyz = (R.Inv())*m_neu;
	//m_xyz(1, 1) += StaCent._X;
	//m_xyz(2, 1) += StaCent._Y;
	//m_xyz(3, 1) += StaCent._Z;

	XYZ.setX(m_xyz(1, 1));
	XYZ.setY(m_xyz(2, 1));
	XYZ.setZ(m_xyz(3, 1));
}

void Coordinate::_xyz2pol()
{
	double pi = 4 * atan(1.0);
	double dis = sqrt(XYZ._X*XYZ._X + XYZ._Y*XYZ._Y);

	Pol.setR(sqrt(dis*dis + XYZ._Z*XYZ._Z));
	Pol.setE(atan2(XYZ._Z, dis));

	if (fabs(Pol._E - pi / 2.0)<1.0e-9)
	{
		if (XYZ._Z>0.0){ Pol.setE(pi / 2.0); }
		if (XYZ._Z < 0.0){ Pol.setE(pi / -2.0); }
		Pol.setA(0.0);
	}
	else
	{
		Pol.setA(atan2(fabs(XYZ._Y), fabs(XYZ._X)));
		if (XYZ._X == 0.0&&XYZ._Y>0.0)
		{
			Pol.setA(pi / 2.0);
		}
		if (XYZ._X == 0.0&&XYZ._Y<0.0)
		{
			Pol.setA(3 * pi / 2.0);
		}
		if (XYZ._X>0.0&&XYZ._Y == 0.0)
		{
			Pol.setA(0.0);
		}
		if (XYZ._X<0.0&&XYZ._Y == 0.0)
		{
			Pol.setA(pi);
		}
		if (XYZ._X>0.0&&XYZ._Y>0.0)
		{
			Pol.setA(Pol._A);
		}
		if (XYZ._X>0.0&&XYZ._Y<0.0)
		{
			Pol.setA(2.0*pi - Pol._A);
		}
		if (XYZ._X<0.0&&XYZ._Y<0.0)
		{
			Pol.setA(pi + Pol._A);
		}
		if (XYZ._X<0.0&&XYZ._Y>0.0)
		{
			Pol.setA(pi - Pol._A);
		}
	}
}

/// Note: pos and vel are used to define the rtn coordinate
/// Parameter:
///            double  *rac       O          ratial along cross
///            double  *vrac      O          v-ratial v-along v-cross
void Coordinate::_xyz2rtn(double* pos, double* vel, double* dxyz, double* rac, double* dvxyz, double* vrac)
{
	double eA[3] = { 0.0 }, eC[3] = { 0.0 }, eR[3] = { 0.0 };

	eA[0] = vel[0];  eA[1] = vel[1];  eA[2] = vel[2];
	Norm3(eA, eA);

	Cross3(pos, vel, eC);
	Norm3(eC, eC);

	Cross3(eA, eC, eR);
	Norm3(eR, eR);

	// Radial Along Cross
	rac[0] = dxyz[0] * eR[0] + dxyz[1] * eR[1] + dxyz[2] * eR[2];
	rac[1] = dxyz[0] * eA[0] + dxyz[1] * eA[1] + dxyz[2] * eA[2];
	rac[2] = dxyz[0] * eC[0] + dxyz[1] * eC[1] + dxyz[2] * eC[2];

	// Velocity
	vrac[0] = dvxyz[0] * eR[0] + dvxyz[1] * eR[1] + dvxyz[2] * eR[2];
	vrac[1] = dvxyz[0] * eA[0] + dvxyz[1] * eA[1] + dvxyz[2] * eA[2];
	vrac[2] = dvxyz[0] * eC[0] + dvxyz[1] * eC[1] + dvxyz[2] * eC[2];

}

void Coordinate::_geo2local(double a, double e, Cart_Crd& StaCent)
{
	Cart_Crd temp;
	Coordinate station(StaCent);
	station._xyz2blh();
	temp.setXYZ((XYZ._X - StaCent._X), (XYZ._Y - StaCent._Y), (XYZ._Z - StaCent._Z));

	NEU.setN(-1 * sin(station.BLH._B)*cos(station.BLH._L)*temp._X - sin(station.BLH._B)*sin(station.BLH._L)*temp._Y + cos(station.BLH._B)*temp._Z);
	NEU.setE(-1 * sin(station.BLH._L)*temp._X + cos(station.BLH._L)*temp._Y);
	NEU.setU(cos(station.BLH._B)*cos(station.BLH._L)*temp._X + cos(station.BLH._B)*sin(station.BLH._L)*temp._Y + sin(station.BLH._B)*temp._Z);

}

