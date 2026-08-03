/* -------------------------------------------------------------------------
* SUPREME_Coordinate.h : Coordinate convention
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2019.01.08
* ------------------------------------------------------------------------- */
#include "SUPREME_CommonFunction.h"
#include "SUPREME_Constant.h"
#include "SUPREME_GnssTime.h"

/* Cartesin Coordinate */
class Cart_Crd
{
public:
	Cart_Crd() :_X(0), _Y(0), _Z(0){}
	Cart_Crd(double x, double y, double z) :_X(x), _Y(y), _Z(z){}
	Cart_Crd(const Cart_Crd& C){ _X = C._X; _Y = C._Y; _Z = C._Z; }

	double _X, _Y, _Z;

	void setX(const double x) { _X = x; }
	void setY(const double y) { _Y = y; }
	void setZ(const double z) { _Z = z; }
	void setXYZ(const double x, const double y, const double z)
	{
		_X = x; _Y = y; _Z = z;
	}
	double getNorm() const{ return sqrt(_X*_X + _Y*_Y + _Z*_Z); }
};

/* Geodetic Coordinate */
class Geo_Crd
{
public:
	Geo_Crd() :_B(0), _L(0), _H(0){}
	Geo_Crd(double b, double l, double h) :_B(b), _L(l), _H(h){}
	Geo_Crd(const Geo_Crd& G){ _B = G._B; _L = G._L; _H = G._H; }

	double _B, _L, _H;

	void setB(const double b) { _B = b; }
	void setL(const double l) { _L = l; }
	void setH(const double h) { _H = h; }
	void setBLH(const double b, const double l, const double h)
	{
		_B = b; _L = l; _H = h;
	}
};

/* Topocentric Terrestrial Coordinate */
class NEU_Crd
{
public:
	NEU_Crd() :_N(0), _E(0), _U(0){}
	NEU_Crd(double n, double e, double u) :_N(n), _E(e), _U(u){}
	NEU_Crd(const NEU_Crd& NEU){ _N = NEU._N; _E = NEU._E; _U = NEU._U; }

	double _N, _E, _U;

	void setN(const double n) { _N = n; }
	void setE(const double e) { _E = e; }
	void setU(const double u) { _U = u; }
	void setNEU(const double n, const double e, const double u)
	{
		_N = n; _E = e; _U = u;
	}
};

/* Polar Coordinate */
class Pol_Crd
{
public:
	Pol_Crd() :_E(0), _A(0), _R(0){}
	Pol_Crd(double e, double a, double r) :_E(e), _A(a), _R(r){}
	Pol_Crd(const Pol_Crd& Pol){ _E = Pol._E; _A = Pol._A; _R = Pol._R; }

	double _E, _A, _R; // Elevation, Azimuth, Radius

	void setE(const double e) { _E = e; }
	void setA(const double a) { _A = a; }
	void setR(const double r) { _R = r; }
	void setEAR(const double e, const double a, const double r)
	{
		_E = e; _A = a; _R = r;
	}
};

/* RTN Coordinate */
class RTN_Crd
{
public:
	RTN_Crd() :_R(0), _A(0), _C(0){}
	RTN_Crd(double r, double a, double c) :_R(r), _A(a), _C(c){}
	RTN_Crd(const RTN_Crd& Pol){ _R = Pol._R; _A = Pol._A; _C = Pol._C; }

	double _R, _A, _C; // Radial Along Cross

	void setR(const double r) { _R = r; }
	void setA(const double a) { _A = a; }
	void setC(const double c) { _C = c; }
	void setRAC(const double r, const double a, const double c)
	{
		_R = r; _A = a; _C = c;
	}
};

class Coordinate
{
public:
	Coordinate(){};
	Coordinate(Cart_Crd xyz) :XYZ(xyz){};       // Cartesin Coordinate
	Coordinate(Geo_Crd blh) :BLH(blh){};        // Geodesy Coordinate
	Coordinate(NEU_Crd neu) :NEU(neu){};        // Topocentric terrestrial Coordinate
	Coordinate(Pol_Crd pol) :Pol(pol){};        // Polar Coordinate
	Coordinate(RTN_Crd arc) :ARC(arc){};        // RTN Coordinate

	Cart_Crd XYZ;
	Geo_Crd BLH;
	NEU_Crd NEU;
	Pol_Crd Pol;
	RTN_Crd ARC;

	void _xyz2blh();
	void _xyz2blh(double a, double e);
	void _blh2xyz();
	void _blh2xyz(double a, double e);
	void _xyz2neu(Cart_Crd& StaCent);
	void _neu2xyz(Cart_Crd& StaCent);
	void _xyz2neu(double a, double e, Cart_Crd& StaCent);
	void _xyz2pol();
	void _xyz2rtn(double* pos, double* vel, double* dxyz, double* rac, double* dvxyz, double* vrac);
	void _rtn2xyz();
	void _geo2local(double a, double e, Cart_Crd& StaCent);
};

void ecef2pos(const double *r, double *pos);

void ecef2enu(double *pos, double *r, double *e);

void xyz2enu(double *pos, double *E);

void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst);
