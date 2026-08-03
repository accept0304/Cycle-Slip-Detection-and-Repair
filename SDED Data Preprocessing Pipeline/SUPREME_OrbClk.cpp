#include <vector>

#include "SUPREME_OrbClk.h"
#include "SUPREME_Coordinate.h"
#include "SUPREME_CommonFunction.h"

int glo_Prn[32] = { 0 };

///% Satellite information class
SatInfo::SatInfo()
{
	m_Size = 0;
	release = 1;
	flag = NULL;
	m_Prn = NULL;
	m_Xs = NULL; m_Vs = NULL;
	m_Clk = NULL; m_Ele = NULL; m_Azi = NULL;
}
SatInfo::SatInfo(unsigned int satn)
{
	this->Initialize(satn);
}

SatInfo::~SatInfo()
{
	this->Clear();
}

int SatInfo::DeleSat(int index)
{
	//
	return 1;
}

int SatInfo::Initialize(unsigned int satn)
{
	if (!release) Clear();

	m_Size = satn;
	if (m_Size > 0)
	{
		flag = (unsigned int*)malloc(sizeof(unsigned int)*satn);
		m_Prn = (unsigned int*)malloc(sizeof(unsigned int)*satn);
		m_Xs = (double*)malloc(sizeof(double)*satn * 3);
		m_Vs = (double*)malloc(sizeof(double)*satn * 3);
		m_Clk = (double*)malloc(sizeof(double)*satn);
		m_Ele = (double*)malloc(sizeof(double)*satn);
		m_Azi = (double*)malloc(sizeof(double)*satn);

		memset(flag, 0, sizeof(unsigned int)*satn);
		memset(m_Prn, 0, sizeof(unsigned int)*satn);
		memset(m_Xs, 0.0, sizeof(double)*satn * 3);
		memset(m_Vs, 0.0, sizeof(double)*satn * 3);
		memset(m_Clk, 0.0, sizeof(double)*satn);
		memset(m_Ele, 0.0, sizeof(double)*satn);
		memset(m_Azi, 0.0, sizeof(double)*satn);

		release = 0;
	}

	return 1;
}

int SatInfo::Clear()
{
	m_Size = 0;
	if (!release)
	{
		if (flag) free(this->flag);
		if (m_Prn) free(this->m_Prn);
		if (m_Xs) free(this->m_Xs); 
		if (m_Vs) free(this->m_Vs);
		if (m_Clk) free(this->m_Clk);
		if (m_Ele) free(this->m_Ele);
		if (m_Azi) free(this->m_Azi);

		release = 1;
	}
	return 1;
}

/// Find Gps navigation ephemeris by prn, time
GpsEph* Find_GPS_NavEph(unsigned int prn, gnsstime& t, NavData* navdata)
{
	Nav_GpsEph* p = Nav_GPS_FindPrn(prn, navdata);
	if (p == NULL) return NULL;

	GpsEph* ptr = p->m_GpsEph;
	GpsEph* des = NULL;

	double t_diff = 9000;       // 9000s = 2.5h

	while (ptr != NULL)
	{
		double dt = fabs(t - ptr->m_Toc_Gps);
		if (dt < t_diff)
		{
			des = ptr;
			t_diff = dt;
		}
		ptr = ptr->m_Next;
	}

	return des;
}

/// Find Glo navigation ephemeris by prn, time
GloEph* Find_GLO_NavEph(unsigned int prn, gnsstime& t, NavData* navdata)
{
	Nav_GloEph* p = Nav_GLO_FindPrn(prn, navdata);
	if (p == NULL) return NULL;

	GloEph* ptr = p->m_GloEph;
	GloEph* des = NULL;

	double t_diff = 9000;          // 9000s = 2.5h

	while (ptr != NULL)
	{
		double dt = fabs(t - ptr->m_Toc_Gps);
		if (dt < t_diff)
		{
			des = ptr;
			t_diff = dt;
		}
		ptr = ptr->m_Next;
	}

	return des;
}

/// Find Gal Ephemeris by prn and time
GalEph* Find_GAL_NavEph(unsigned int prn, gnsstime& t, NavData* navdata)
{
	Nav_GalEph* p = Nav_GAL_FindPrn(prn, navdata);
	if (p == NULL) return NULL;
	GalEph* ptr = p->m_GalEph;
	GalEph* des = NULL;

	double t_diff = 9000;         // 9000s = 2.5h

	while (ptr != NULL)
	{
		double dt = fabs(t - ptr->m_Toc_Gps);
		if (dt < t_diff)
		{
			des = ptr;
			t_diff = dt;
		}
		ptr = ptr->m_Next;
	}
	return des;
}

/// Find Bds navigation ephememris by prn and time
BdsEph* Find_BDS_NavEph(unsigned int prn, gnsstime& t, NavData* navdata)
{
	Nav_BdsEph* p = Nav_BDS_FindPrn(prn, navdata);
	if (p == NULL) return NULL;

	BdsEph* ptr = p->m_BdsEph;
	BdsEph* des = NULL;

	double t_diff = 9000;       // 9000s = 2.5h

	while (ptr != NULL)
	{
		double dt = fabs(t - ptr->m_Toc_Gps);
		if (dt < t_diff)
		{
			des = ptr;
			t_diff = dt;
		}
		ptr = ptr->m_Next;
	}

	return des;
}

/// Find Navigation Ephemeris by prn and time.(Gps,Glo,Gal,Bds)
void* Find_NavEph(unsigned int prn, gnsstime& t, NavData* navdata)
{
	void* re = NULL;

	if (GetSystem(prn) == 'G')
	{
		re = (void*)(Find_GPS_NavEph(prn, t, navdata));
	}
	else if (GetSystem(prn) == 'R')
	{
		re = (void*)(Find_GLO_NavEph(prn, t, navdata));
	}
	else if (GetSystem(prn) == 'E')
	{
		re = (void*)(Find_GAL_NavEph(prn, t, navdata));
	}
	else if (GetSystem(prn) == 'C')
	{
		re = (void*)(Find_BDS_NavEph(prn, t, navdata));
	}
	else
		re = NULL;

	return re;
}

double motion_equationX(double x, double y, double z, double Vy, double Ax)
{
	double r = sqrt(x*x + y*y + z*z);
	return -1 * E_GM*x / (r*r*r) - 1.5*GLO_J20*E_GM*PZ90_A*PZ90_A*x / (pow(r, 5))*(1 - 5 * z*z / r / r)
		+ OMEGA_GLO*OMEGA_GLO*x + 2 * OMEGA_GLO*Vy + Ax;
}
double motion_equationY(double x, double y, double z, double Vx, double Ay)
{
	double r = sqrt(x*x + y*y + z*z);
	return -1 * E_GM*y / (r*r*r) - 1.5*GLO_J20*E_GM*PZ90_A*PZ90_A*y / (pow(r, 5))*(1 - 5 * z*z / r / r)
		+ OMEGA_GLO*OMEGA_GLO*y - 2 * OMEGA_GLO*Vx + Ay;
}
double motion_equationZ(double x, double y, double z, double Az)
{
	double r = sqrt(x*x + y*y + z*z);
	return -1 * E_GM*z / (r*r*r) - 1.5*GLO_J20*E_GM*PZ90_A*PZ90_A*z / (pow(r, 5))*(3 - 5 * z*z / r / r)
		+ Az;
}

void get_klm(double P[3], double V[3], double A[3], double klm[3])
{
	klm[0] = motion_equationX(P[0], P[1], P[2], V[1], A[0]);
	klm[1] = motion_equationY(P[0], P[1], P[2], V[0], A[1]);
	klm[2] = motion_equationZ(P[0], P[1], P[2], A[2]);
}
void get_npq(double k, double l, double m, double V[3], double npq[3])
{
	npq[0] = V[0] + k;
	npq[1] = V[1] + l;
	npq[2] = V[2] + m;
}

/// Runge Kutta Integral
void RungeKutta(double P[3], double V[3], double A[3], double stride)
{
	double klm1[3] = { 0 }, klm2[3] = { 0 }, klm3[3] = { 0 }, klm4[3] = { 0 };
	double npq1[3] = { 0 }, npq2[3] = { 0 }, npq3[3] = { 0 }, npq4[3] = { 0 };

	get_klm(P, V, A, klm1);
	get_npq(0, 0, 0, V, npq1);

	double P1[3] = { P[0] + stride*npq1[0] / 2, P[1] + stride*npq1[1] / 2, P[2] + stride*npq1[2] / 2 };
	double V1[3] = { V[0] + stride*klm1[0] / 2, V[1] + stride*klm1[1] / 2, V[2] + stride*klm1[2] / 2 };
	get_klm(P1, V1, A, klm2);
	get_npq(stride*klm1[0] / 2, stride*klm1[1] / 2, stride*klm1[2] / 2, V, npq2);

	double P2[3] = { P[0] + stride*npq2[0] / 2, P[1] + stride*npq2[1] / 2, P[2] + stride*npq2[2] / 2 };
	double V2[3] = { V[0] + stride*klm2[0] / 2, V[1] + stride*klm2[1] / 2, V[2] + stride*klm2[2] / 2 };
	get_klm(P2, V2, A, klm3);
	get_npq(stride*klm2[0] / 2, stride*klm2[1] / 2, stride*klm2[2] / 2, V, npq3);

	double P3[3] = { P[0] + stride*npq3[0], P[1] + stride*npq3[1], P[2] + stride*npq3[2] };
	double V3[3] = { V[0] + stride*klm3[0], V[1] + stride*klm3[1], V[2] + stride*klm3[2] };
	get_klm(P3, V3, A, klm4);
	get_npq(stride*klm3[0], stride*klm3[1], stride*klm3[2], V, npq4);

	V[0] += stride*(klm1[0] + 2 * klm2[0] + 2 * klm3[0] + klm4[0]) / 6;
	V[1] += stride*(klm1[1] + 2 * klm2[1] + 2 * klm3[1] + klm4[1]) / 6;
	V[2] += stride*(klm1[2] + 2 * klm2[2] + 2 * klm3[2] + klm4[2]) / 6;
	P[0] += stride*(npq1[0] + 2 * npq2[0] + 2 * npq3[0] + npq4[0]) / 6;
	P[1] += stride*(npq1[1] + 2 * npq2[1] + 2 * npq3[1] + npq4[1]) / 6;
	P[2] += stride*(npq1[2] + 2 * npq2[2] + 2 * npq3[2] + npq4[2]) / 6;
}

/// Calculate GPS Satellite Position and Velocity by Keplerian elements
int Cal_GPS_PosVel(gnsstime& t, GpsEph* gpseph, double pos[3], double vel[3])
{
	if (t.m_Week == 0)
		t._date2gpst();
	double tk = t.m_SecondofWeek - gpseph->m_TOE;
	if (tk > 302400)
		tk -= 604800;
	else
	if (tk < -302400)
		tk += 604800;

	double n = sqrt(E_GM_GPS) / pow(gpseph->m_sqrt_A, 3) + gpseph->m_Delta_n;
	double Mk = gpseph->m_M0 + n*tk;

	double Ek = Mk, dek = 0.0, temp = 0.0;
	do
	{
		temp = Mk + gpseph->m_e*sin(Ek);
		dek = temp - Ek;
		Ek = temp;
	} while (fabs(dek)>1.0e-8);

	double fk = 2 * atan(sqrt((1 + gpseph->m_e) / (1 - gpseph->m_e))*tan(Ek / 2.0));
	//double fk = atan2(sqrt(1.0 - gpseph->m_e*gpseph->m_e)*sin(Ek), cos(Ek) - gpseph->m_e);
	if (fk < 0.0)
		fk += 2 * PI;

	if (fk<0.0)
	{
		fk += 2 * PI;
	}

	double FFk = fk + gpseph->m_omega;
	double du = gpseph->m_Cuc*cos(2 * FFk) + gpseph->m_Cus*sin(2 * FFk);
	double dr = gpseph->m_Crc*cos(2 * FFk) + gpseph->m_Crs*sin(2 * FFk);
	double di = gpseph->m_Cic*cos(2 * FFk) + gpseph->m_Cis*sin(2 * FFk);
	double uk = FFk + du;
	double rk = pow(gpseph->m_sqrt_A, 2)*(1 - gpseph->m_e*cos(Ek)) + dr;
	double ik = gpseph->m_i0 + di + gpseph->m_IDOT*tk;

	double xk = rk*cos(uk);
	double yk = rk*sin(uk);
	double wk = gpseph->m_OMEGA0 + (gpseph->m_OMEGADOT - OMEGA_GPS)*tk - OMEGA_GPS*gpseph->m_TOE;

	pos[0] = xk*cos(wk) - yk*cos(ik)*sin(wk);
	pos[1] = xk*sin(wk) + yk*cos(ik)*cos(wk);
	pos[2] = yk*sin(ik);

	double  vEk = n / (1 - gpseph->m_e*cos(Ek));
	double vFFk = sqrt((1 + gpseph->m_e) / (1 - gpseph->m_e))*cos(fk / 2.0)*cos(fk / 2.0)*vEk / cos(Ek / 2.0) / cos(Ek / 2.0);
	double vdu = (1 + 2 * (gpseph->m_Cus*cos(2 * FFk) - gpseph->m_Cuc*sin(2 * FFk)))*vFFk;
	double vrk = gpseph->m_sqrt_A*gpseph->m_sqrt_A*gpseph->m_e*sin(Ek)*vEk + 2 * (gpseph->m_Crs*cos(2 * FFk) - gpseph->m_Crc*sin(2 * FFk))*vFFk;
	double vik = 2 * (gpseph->m_Cis*cos(2 * FFk) - gpseph->m_Cic*sin(2 * FFk))*vFFk + gpseph->m_IDOT;
	double vwk = gpseph->m_OMEGADOT - OMEGA_GPS;

	double vxk = vrk*cos(uk) - rk*sin(uk)*vdu;
	double vyk = vrk*sin(uk) + rk*cos(uk)*vdu;

	vel[0] = vxk*cos(wk) - vyk*cos(ik)*sin(wk) + yk*sin(wk)*sin(ik)*vik - (xk*sin(wk) + yk*cos(wk)*cos(ik))*vwk;
	vel[1] = vxk*sin(wk) + vyk*cos(ik)*cos(wk) - yk*cos(wk)*sin(ik)*vik + (xk*cos(wk) - yk*sin(wk)*cos(ik))*vwk;
	vel[2] = vyk*sin(ik) + yk*cos(ik)*vik;

	return 1;
}

/// Calclulate GLO Satellite Position and Velocity
int Cal_GLO_PosVel(gnsstime& t, GloEph* gloeph, double pos[3], double vel[3])
{
	if (t.m_Week == 0)
		t._date2gpst();

	double timediff = t - gloeph->m_Toc_Gps;

	int Init = 6;
	double stride = timediff;
	double x1[3] = { 0.0 }, x2[3] = { 0.0 };

	do
	{
		double P1[3] = { gloeph->m_X, gloeph->m_Y, gloeph->m_Z };
		double V1[3] = { gloeph->m_Vx, gloeph->m_Vy, gloeph->m_Vz };
		double A1[3] = { gloeph->m_Ax, gloeph->m_Ay, gloeph->m_Az };

		stride = timediff / Init;
		RungeKutta(P1, V1, A1, stride);
		x1[0] = P1[0];
		x1[1] = P1[1];
		x1[2] = P1[2];

		double P2[3] = { gloeph->m_X, gloeph->m_Y, gloeph->m_Z };
		double V2[3] = { gloeph->m_Vx, gloeph->m_Vy, gloeph->m_Vz };
		double A2[3] = { gloeph->m_Ax, gloeph->m_Ay, gloeph->m_Az };
		stride /= 2.0;

		RungeKutta(P2, V2, A2, stride);
		RungeKutta(P2, V2, A2, stride);
		x2[0] = P2[0];
		x2[1] = P2[1];
		x2[2] = P2[2];
		Init *= 2;

	} while (fabs(x2[0] - x1[0]) > 0.001 || fabs(x2[1] - x1[1]) > 0.001 || fabs(x2[2] - x1[2]) > 0.001);

	int n = (int)(timediff / stride), i = 0;
	double m = timediff - n*stride;

	n = abs(n);

	double A[3] = { gloeph->m_Ax, gloeph->m_Ay, gloeph->m_Az };

	pos[0] = gloeph->m_X; pos[1] = gloeph->m_Y; pos[2] = gloeph->m_Z;
	vel[0] = gloeph->m_Vx; vel[1] = gloeph->m_Vy; vel[2] = gloeph->m_Vz;

	for (i = 0; i < n; ++i)
	{
		RungeKutta(pos, vel, A, stride);
	}

	if (fabs(m)>1.0e-8)
		RungeKutta(pos, vel, A, m);

	return 1;
}

/// Calculate GAL Satellite Position and Velocity
int Cal_GAL_PosVel(gnsstime& t, GalEph* galeph, double pos[3], double vel[3])
{
	if (t.m_SecondofWeek == 0)
		t._date2gpst();

	double tk = t.m_SecondofWeek - galeph->m_TOE;
	if (tk > 302400)
		tk -= 604800;
	else
	if (tk < -302400)
		tk += 604800;

	double n = sqrt(E_GM) / pow(galeph->m_sqrt_A, 3) + galeph->m_Delta_n;
	double Mk = galeph->m_M0 + n*tk;

	double Ek = Mk, dek = 0.0, temp = 0.0;
	do
	{
		temp = Mk + galeph->m_e*sin(Ek);
		dek = temp - Ek;
		Ek = temp;
	} while (fabs(dek) > 1.0e-8);

	double  fk = 2 * atan(sqrt((1 + galeph->m_e) / (1 - galeph->m_e))*tan(Ek / 2.0));
	//double fk = atan2(sqrt(1.0 - galeph->m_e*galeph->m_e)*sin(Ek), cos(Ek) - galeph->m_e);
	if (fk < 0.0)
		fk += 2 * PI;

	double FFk = fk + galeph->m_omega;
	double du = galeph->m_Cuc*cos(2 * FFk) + galeph->m_Cus*sin(2 * FFk);
	double dr = galeph->m_Crc*cos(2 * FFk) + galeph->m_Crs*sin(2 * FFk);
	double di = galeph->m_Cic*cos(2 * FFk) + galeph->m_Cis*sin(2 * FFk);
	double uk = FFk + du;
	double rk = pow(galeph->m_sqrt_A, 2)*(1 - galeph->m_e*cos(Ek)) + dr;
	double ik = galeph->m_i0 + di + galeph->m_IDOT*tk;

	double xk = rk*cos(uk);
	double yk = rk*sin(uk);
	double wk = galeph->m_OMEGA0 + (galeph->m_OMEGADOT - OMEGA_GPS)*tk - OMEGA_GPS*galeph->m_TOE;

	pos[0] = xk*cos(wk) - yk*cos(ik)*sin(wk);
	pos[1] = xk*sin(wk) + yk*cos(ik)*cos(wk);
	pos[2] = yk*sin(ik);

	double  vEk = n / (1 - galeph->m_e*cos(Ek));
	double vFFk = sqrt((1 + galeph->m_e) / (1 - galeph->m_e))*cos(fk / 2.0)*cos(fk / 2.0)*vEk / cos(Ek / 2.0) / cos(Ek / 2.0);
	double  vdu = (1 + 2 * (galeph->m_Cus*cos(2 * FFk) - galeph->m_Cuc*sin(2 * FFk)))*vFFk;
	double vrk = galeph->m_sqrt_A*galeph->m_sqrt_A*galeph->m_e*sin(Ek)*vEk + 2 * (galeph->m_Crs*cos(2 * FFk) - galeph->m_Crc*sin(2 * FFk))*vFFk;
	double vik = 2 * (galeph->m_Cis*cos(2 * FFk) - galeph->m_Cic*sin(2 * FFk))*vFFk + galeph->m_IDOT;
	double vwk = galeph->m_OMEGADOT - OMEGA_GPS;

	double vxk = vrk*cos(uk) - rk*sin(uk)*vdu;
	double vyk = vrk*sin(uk) + rk*cos(uk)*vdu;

	vel[0] = vxk*cos(wk) - vyk*cos(ik)*sin(wk) + yk*sin(wk)*sin(ik)*vik - (xk*sin(wk) + yk*cos(wk)*cos(ik))*vwk;
	vel[1] = vxk*sin(wk) + vyk*cos(ik)*cos(wk) - yk*cos(wk)*sin(ik)*vik + (xk*cos(wk) - yk*sin(wk)*cos(ik))*vwk;
	vel[2] = vyk*sin(ik) + yk*cos(ik)*vik;

	return 1;
}

/// Calculate BDS Satellite Position and Velocity
int Cal_BDS_PosVel(gnsstime& t, BdsEph* bdseph, double pos[3], double vel[3])
{
	if (t.m_Week == 0)
		t._date2gpst();

	double tow = t.m_SecondofWeek - 14.0;
	if (tow<0)
	{
		tow += 604800.0;
	}
	double tk = tow - bdseph->m_TOE;
	if (tk > 302400)
		tk -= 604800;
	else
	if (tk < -302400)
		tk += 604800;

	double n = sqrt(E_GM) / pow(bdseph->m_sqrt_A, 3) + bdseph->m_Delta_n;
	double Mk = bdseph->m_M0 + n*tk;

	double Ek = Mk, dek = 0.0, temp = 0.0;
	do
	{
		temp = Mk + bdseph->m_e*sin(Ek);
		dek = temp - Ek;
		Ek = temp;
	} while (fabs(dek) > 1.0e-8);

	double fk = 2 * atan(sqrt((1 + bdseph->m_e) / (1 - bdseph->m_e))*tan(Ek / 2.0));
	//double fk = atan2(sqrt(1.0 - bdseph->m_e*bdseph->m_e)*sin(Ek), cos(Ek) - bdseph->m_e);

	if (fk < 0.0)
		fk += 2 * PI;

	double FFk = fk + bdseph->m_omega;
	double du = bdseph->m_Cuc*cos(2 * FFk) + bdseph->m_Cus*sin(2 * FFk);
	double dr = bdseph->m_Crc*cos(2 * FFk) + bdseph->m_Crs*sin(2 * FFk);
	double di = bdseph->m_Cic*cos(2 * FFk) + bdseph->m_Cis*sin(2 * FFk);
	double uk = FFk + du;
	double rk = pow(bdseph->m_sqrt_A, 2)*(1 - bdseph->m_e*cos(Ek)) + dr;
	double ik = bdseph->m_i0 + di + bdseph->m_IDOT*tk;

	double xk = rk*cos(uk);
	double yk = rk*sin(uk);
	double wk = bdseph->m_OMEGA0 + (bdseph->m_OMEGADOT - OMEGA_BDS)*tk - OMEGA_BDS*bdseph->m_TOE;

	int GEO_Flag = 0;
	if (bdseph->m_prn >= MIN_BDS_PRN && bdseph->m_prn < MIN_BDS_PRN + 5)
	{
		wk = bdseph->m_OMEGA0 + bdseph->m_OMEGADOT*tk - OMEGA_BDS*bdseph->m_TOE;
		GEO_Flag = 1;
	}

	pos[0] = xk*cos(wk) - yk*cos(ik)*sin(wk);
	pos[1] = xk*sin(wk) + yk*cos(ik)*cos(wk);
	pos[2] = yk*sin(ik);

	double  vEk = n / (1 - bdseph->m_e*cos(Ek));
	double vFFk = sqrt((1 + bdseph->m_e) / (1 - bdseph->m_e))*cos(fk / 2.0)*cos(fk / 2.0)*vEk / cos(Ek / 2.0) / cos(Ek / 2.0);
	double  vdu = (1 + 2 * (bdseph->m_Cus*cos(2 * FFk) - bdseph->m_Cuc*sin(2 * FFk)))*vFFk;
	double vrk = bdseph->m_sqrt_A*bdseph->m_sqrt_A*bdseph->m_e*sin(Ek)*vEk + 2 * (bdseph->m_Crs*cos(2 * FFk) - bdseph->m_Crc*sin(2 * FFk))*vFFk;
	double vik = 2 * (bdseph->m_Cis*cos(2 * FFk) - bdseph->m_Cic*sin(2 * FFk))*vFFk + bdseph->m_IDOT;
	double vwk = bdseph->m_OMEGADOT - OMEGA_BDS;

	if (GEO_Flag)
		vwk = bdseph->m_OMEGADOT;

	double vxk = vrk*cos(uk) - rk*sin(uk)*vdu;
	double vyk = vrk*sin(uk) + rk*cos(uk)*vdu;

	double coswk = cos(wk), sinwk = sin(wk);
	double cosik = cos(ik), sinik = sin(ik);

	vel[0] = vxk*coswk - vyk*cosik*sinwk + yk*sinwk*sinik*vik - (xk*sinwk + yk*coswk*cosik)*vwk;
	vel[1] = vxk*sinwk + vyk*cosik*coswk - yk*coswk*sinik*vik + (xk*coswk - yk*sinwk*cosik)*vwk;
	vel[2] = vyk*sinik + yk*cosik*vik;

	if (GEO_Flag)
	{
		double tp1 = -5.0*PI / 180.0;
		double tp2 = OMEGA_BDS*tk;

		double tx = pos[0], ty = pos[1], tz = pos[2];
		double vtx = vel[0], vty = vel[1], vtz = vel[2];

		double cost1 = cos(tp1), sint1 = sin(tp1);
		double cost2 = cos(tp2), sint2 = sin(tp2);

		pos[0] = cost2*tx + sint2*cost1*ty + sint2*sint1*tz;
		pos[1] = -1 * sint2*tx + cost2*cost1*ty + cost2*sint1*tz;
		pos[2] = -1 * sint1*ty + cost1*tz;

		vel[0] = cost2*vtx + sint2*cost1*vty + sint2*sint1*vtz;
		vel[1] = -1 * sint2*vtx + cost2*cost1*vty + cost2*sint1*vtz;
		vel[2] = -1 * sint1*vty + cost1*vtz;
	}

	return 1;
}

/// Calculate GPS Clock Offsets
double Cal_GPS_Clk(gnsstime& t, GpsEph* gpseph)
{
	//double dt = t.m_SecondofWeek - gpseph->m_Toc_Gps.m_SecondofWeek;

	//if (dt > 302400)
	//	dt -= 604800;
	//else
	//if (dt < -302400)
	//	dt += 604800;

	double dt = t - gpseph->m_Toc_Gps;

	return gpseph->m_clock_bias + gpseph->m_clock_drift*dt + gpseph->m_clock_driftrate*dt*dt;
}

/// Calculate Glonass Clock Offset
double Cal_GLO_Clk(gnsstime& t, GloEph *gloeph)
{
	double dt = t.m_SecondofWeek - gloeph->m_Toc_Gps.m_SecondofWeek;
	if (dt>302400)
		dt -= 604800;
	else
	if (dt < -302400)
		dt += 604800;

	return gloeph->m_tau - gloeph->m_gamma*dt;
}

/// Calculate Galileo Clock Offset
double Cal_GAL_Clk(gnsstime& t, GalEph* galeph)
{
	double dt = t.m_SecondofWeek - galeph->m_Toc_Gps.m_SecondofWeek;
	if (dt>302400)
		dt -= 604800;
	else
	if (dt < -302400)
		dt += 604800;

	return galeph->m_clock_bias + galeph->m_clock_drift*dt + galeph->m_clock_driftrate*dt*dt;
}

/// Calculate BDS Clock Offset
double Cal_BDS_Clk(gnsstime& t, BdsEph* bdseph)
{
	double dt = t.m_SecondofWeek - bdseph->m_Toc_Gps.m_SecondofWeek;
	if (dt>302400)
		dt -= 604800;
	else
	if (dt<-302400)
		dt += 604800;

	return bdseph->m_clock_bias + bdseph->m_clock_drift*dt + bdseph->m_clock_driftrate*dt*dt;
}

double Precise_Interp_Clock(gnsstime& time, Clk_EpochData *ptr, int len, int n)
{
	if (ptr == NULL || len == 0 || len != (n + 1))
	{
		return 999999;
	}
	if (time.m_Jd<ptr[0].m_time.m_Jd || time.m_Jd>ptr[len - 1].m_time.m_Jd)
	{
		return 999999;
	}

	double satclk = 0.0;

	int i = 0;

	for (i = 0; i != len; ++i)
	{
		double temp = 1.0;

		int j = 0;

		for (j = 0; j != len; ++j)
		{
			if (i == j)
			{
				continue;
			}
			temp *= ((time - ptr[j].m_time) / (ptr[i].m_time - ptr[j].m_time));
		}
		satclk += (temp*ptr[i].m_clockcor);
	}

	return satclk;
}

/// Convert navigation positon and velocity  to precise position and velocity by realtime product
int Nav2Pre(double pos[3], double vel[3], double dr, double da, double dc, double dx[3])
{
	double eA[3] = { 0.0 }, eC[3] = { 0.0 }, eR[3] = { 0.0 };

	eA[0] = vel[0];  eA[1] = vel[1];  eA[2] = vel[2];

	Norm3(eA);

	Cross3(pos, vel, eC);

	Norm3(eC);

	Cross3(eA, eC, eR);

	////Norm3(eR);

	dx[0] = eR[0] * dr + eA[0] * da + eC[0] * dc;
	dx[1] = eR[1] * dr + eA[1] * da + eC[1] * dc;
	dx[2] = eR[2] * dr + eA[2] * da + eC[2] * dc;

	return 1;
}

/// Calculate Satellite Position and Velocity by Time ---------------------------
/* Parameters:
* void           *eph               I              Navigation ephemeris
* gnsstime       &obst              I              Gnss time
* double         pos[3]             O              Satellite positions
* double         vel[3]             O              Satellite velocity
* double         &SatClk            O              Satellite clock offset
* char           Sys                I              GNSS system
* unsigned int   flag               I              Use orb clk flag
* ------------------------------------------------------------------------------ */
double Cal_Satellite_PosVel_byTime(void* eph, gnsstime& obst, double pos[3], double vel[3],
	double& SatClk, char Sys, unsigned int flag, double& TGD1, double& TGD2)
{
	double Omega = OMEGA_GPS;

	switch (Sys)
	{
		case 'G':
		{
					GpsEph *ptr = (GpsEph*)eph;

					Cal_GPS_PosVel(obst, ptr, pos, vel);     // Get GPS satellite Pos and Vel by navigation
					SatClk = Cal_GPS_Clk(obst, ptr);         // Get GPS satellite clock by navigation
					TGD1 = ptr->m_TGD;
					Omega = OMEGA_GPS;

					break;
		}
		case 'R':
		{
					GloEph *ptr = (GloEph*)(eph);

					Cal_GLO_PosVel(obst, ptr, pos, vel);
					SatClk = Cal_GLO_Clk(obst, ptr);

					Omega = OMEGA_GLO;

					break;
		}
		case 'C':
		{
					BdsEph *ptr = (BdsEph*)(eph);

					Cal_BDS_PosVel(obst, ptr, pos, vel);
					SatClk = Cal_BDS_Clk(obst, ptr);
					TGD1 = ptr->m_TGD1;
					TGD2 = ptr->m_TGD2;
					Omega = OMEGA_BDS;

					break;
		}
		case 'E':
		{
					GalEph *ptr = (GalEph*)(eph);

					Cal_GAL_PosVel(obst, ptr, pos, vel);
					SatClk = Cal_GAL_Clk(obst, ptr);

					Omega = OMEGA_GPS;

					break;
		}
		default:break;
	}

	return Omega;
}

/// Realtime Calculate Satellite Position and Velocity --------------------------------------
/* Parameter:
* unsigned int   prn
* gnsstime       t                    observation time
* double         pr                   pesudo range
* double         Sta[3]               statioin position
* void          *eph                  ephemeris
* double         pos[3]               satellite position
* double         vel[3]               satellite velocity
* double         SatClk               satellite clock offset
* int            flag                 1:use orbclk corrections   0:not use orbclk corrections
* ------------------------------------------------------------------------------------------- */
int Cal_Sat_PosVel_Real(unsigned int prn, gnsstime& t, double pr, double Sta[3], void* eph,
	double pos[3], double vel[3], double& SatClk, int flag, double& TGD1, double& TGD2)
{
	double Rclk = 0.0, dif_t = 100.0, Mjd = 0.0, old_dt = 0.0;
	int Initialize = 0, loopCount = 0;
	double Omega = OMEGA_GPS;

	gnsstime obst;

	char Sys = GetSystem(prn);

	while (fabs(dif_t) > 1.0e-11)
	{
		obst = t;

		double dt = pr / LIGHTSPEED + Rclk;
		double dt_ro = pr / LIGHTSPEED;
		Mjd = obst.m_Jd - 2400000.5 - dt / 86400.0;
		obst.m_Jd = Mjd + 2400000.5;

		dif_t = dt - old_dt;

		obst.m_SecondofWeek -= dt;
		if (obst.m_SecondofWeek < 0)
		{
			obst.m_SecondofWeek += 604800;
			--obst.m_Week;
		}

		Omega = Cal_Satellite_PosVel_byTime(eph, obst, pos, vel, SatClk, Sys, flag, TGD1, TGD2);
		if (Omega < 0) return 0;

		// Sagnac Influence of earth rotation
		//double mg = Omega*dt_ro;
		//double x = pos[0], y = pos[1];
		//pos[0] = x*cos(mg) + y*sin(mg);
		//pos[1] = y*cos(mg) - x*sin(mg);

		// Distance between station and satellite
		double Receiver_Satellite = sqrt((pos[0] - Sta[0])*(pos[0] - Sta[0])
			+ (pos[1] - Sta[1])*(pos[1] - Sta[1])
			+ (pos[2] - Sta[2])*(pos[2] - Sta[2]));

		if (Initialize == 0)
		{
			Rclk = (pr - Receiver_Satellite) / LIGHTSPEED + (SatClk);
			Initialize = 1;
		}

		pr = Receiver_Satellite;

		old_dt = dt;

		if (loopCount > 15)	return 0;
		loopCount++;
	}

	// return signal launch time
	//obst.m_Jd = Mjd + 2400000.5;
	t = obst;

	return 1;
}

/// Post Calculate Satellite Position and Velocity --------------------------------------
/* Parameter:
* unsigned int      prn                 I                satellite prn
* gnsstime          &t                  I                observation time
* double            pr                  I                satellite to receiver range
* double            Sta[3]              I                station position
* PreEpochData      *ptr_pre            I                precise ephemeris data
* Clk_EpochData     *ptr_clk            I                precise clock data
* double            satpos[3]           O                satellite position
* double            satvel[3]           O                satellite velocity
* double            &SatClk             O                satellite clock offset
* ------------------------------------------------------------------------------------------- */
int Cal_Sat_PosVel_Post(unsigned int prn, gnsstime& t, double pr, double Sta[3],
	PreEpochData* ptr_pre, int len_pre, Clk_EpochData* ptr_clk, int len_clk,
	double satpos[3], double satvel[3], double& SatClk)
{
	double Rclk = 0.0, dif_t = 100.0, Mjd = 0.0, old_dt = 0.0;
	int Initialize = 0, loopCount = 0;
	double Omega = OMEGA_GPS;

	gnsstime obst;

	char Sys = 0;
	unsigned int prn_tmp = 0;
	Sys = GetSysPrn(prn, prn_tmp);

	while (fabs(dif_t) > 1.0e-8)
	{
		obst = t;

		double dt = pr / LIGHTSPEED + Rclk;
		double dt_ro = pr / LIGHTSPEED;
		Mjd = obst.m_Jd - 2400000.5 - dt / 86400.0;

		dif_t = dt - old_dt;

		obst.m_SecondofWeek -= dt;
		if (obst.m_SecondofWeek < 0)
		{
			obst.m_SecondofWeek += 604800;
			--obst.m_Week;
		}

		if (Precise_Interp_Pos(obst, ptr_pre, len_pre, satpos, len_pre - 1) == 0) return 0;
		if (Precise_Interp_Vel(obst, ptr_pre, len_pre, satvel, len_pre - 1) == 0) return 0;
		SatClk = Precise_Interp_Clock(obst, ptr_clk, len_clk, len_clk - 1);
		if (SatClk == 999999) return 0;

		// Sagnac Influence of earth rotation
		/*double mg = Omega*dt_ro;
		double x = satpos[0], y = satpos[1];
		satpos[0] = x*cos(mg) + y*sin(mg);
		satpos[1] = y*cos(mg) - x*sin(mg);*/

		// Distance between station and satellite
		double Receiver_Satellite = sqrt((satpos[0] - Sta[0])*(satpos[0] - Sta[0])
			+ (satpos[1] - Sta[1])*(satpos[1] - Sta[1])
			+ (satpos[2] - Sta[2])*(satpos[2] - Sta[2]));

		if (Initialize == 0)
		{
			Rclk = (pr - Receiver_Satellite) / LIGHTSPEED + (SatClk);
			Initialize = 1;
		}

		pr = Receiver_Satellite;

		old_dt = dt;

		loopCount++;
		if (loopCount > 15)
			break;
	}

	// return signal launch time
	obst.m_Jd = Mjd + 2400000.5;
	t = obst;

	return 1;
}

/// Calculate Satellite Elevation and Azimuth -----------------------------
/* Parameter:
* double   sta[3]             I              station position
* double   sat[3]             I              satellite position
* double   el                 O              satellite elevation
* double   azimuth            O              satellite azimuth
* ------------------------------------------------------------------------- */
int Cal_Sat_EleAzimuth(double sta[3], double sat[3], double& el, double& azimuth)
{
	Cart_Crd sta_cor(sta[0], sta[1], sta[2]);
	Cart_Crd sat_cor(sat[0], sat[1], sat[2]);
	Coordinate local(sat_cor);
	Coordinate polar;

	local._xyz2neu(WGS84_A, WGS84_E, sta_cor);
	polar.XYZ.setX(local.NEU._N);
	polar.XYZ.setY(local.NEU._E);
	polar.XYZ.setZ(local.NEU._U);

	polar._xyz2pol();

	el = polar.Pol._E;
	azimuth = polar.Pol._A;

	return 1;
}

/// Get earth shadow or not -------------------------------
/* Parameter:
* double    sun[3]              sun position
* double    sat[3]              satellite position
* char     *sat_type            satellite type
* --------------------------------------------------------- */
int EarthShadow(double sun[3], double sat[3], char *sat_type)
{
	int i = 0;

	double rs = Norm(sun, 3);
	double ss = Norm(sat, 3);
	double ca = 0.0;

	if (sat_type == NULL || strcmp(sat_type, "IIA") != 0)
	{
		return 1;
	}

	double es[3] = { 0.0 };

	for (i = 0; i<3; ++i)
	{
		es[i] = sun[i] / rs;
		ca += es[i] * sat[i];
	}

	ca /= ss;

	if (ca<-1.0)
	{
		ca = -1.0;
	}
	if (ca>1.0)
	{
		ca = 1.0;
	}

	double aca = acos(ca);

	if (aca<PI / 2 || ss*sin(aca)>WGS84_A)
	{
		return 1;
	}
	else
	{
		return 0;    // Shadow
	}
}
