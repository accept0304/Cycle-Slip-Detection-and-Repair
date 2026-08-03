#include "SUPREME_CommonFunction.h"
#include "SUPREME_Troposphere.h"
#include "SUPREME_Matrix.h"

int ipow(int base, int exp)
{
	int result = 1;

	while (exp) {
		if (exp & 1) result *= base;
		exp >>= 1;
		base *= base;
	}
	return result;
}

/* hydro-ave-a,b,c, hydro-amp-a,b,c, wet-a,b,c at latitude 15,30,45,60,75 */
static const double coef[][5] =
{
	{ 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3 },
	{ 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3 },
	{ 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3 },

	{ 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5 },
	{ 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5 },
	{ 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5 },

	{ 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4 },
	{ 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3 },
	{ 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2 }
};
/* height correction */
static const double aht[] = { 2.53E-5, 5.49E-3, 1.14E-3 };

static double interpc(const double coef[], double lat)
{
	int i = (int)(lat / 15.0);
	if (i<1) return coef[0]; else if (i>4) return coef[4];
	return coef[i - 1] * (1.0 - lat / 15.0 + i) + coef[i] * (lat / 15.0 - i);
}

static double mapf(double el, double a, double b, double c)
{
	double sinel = sin(el);
	return (1.0 + a / (1.0 + b / (1.0 + c))) / (sinel + (a / (sinel + b / (sinel + c))));
}


/// NMF Mapping Function -----------------------------------------------------------
/* Parameter:
*            double    h             I      height
*            double    elev          I      satellite elevation
*            double    latitude      I      receiver latitude
*            int       doy           I      day of year
*            double    dry_me        O      trop dry mapping function
*            double    wet_me        O      trop wet mapping function
* Ref. RTKLIB 2.4.3 : rtkcmn.c nmf()
* --------------------------------------------------------------------------------- */
void NMF(int doy, double latitude, double h, double elev, double& dry_me, double& wet_me)
{
	int i;
	double y, cosy, ah[3], aw[3], dm, lat = latitude * R2D;

	if (elev <= 0.0)
		return;

	/* year from doy 28, added half a year for southern latitudes */
	y = (doy - 28.0) / 365.25 + (lat<0.0 ? 0.5 : 0.0);

	cosy = cos(2.0*PI*y);
	lat = fabs(lat);

	for (i = 0; i<3; i++)
	{
		ah[i] = interpc(coef[i], lat) - interpc(coef[i + 3], lat)*cosy;
		aw[i] = interpc(coef[i + 6], lat);
	}
	/* ellipsoidal height is used instead of height above sea level */
	dm = (1.0 / sin(elev) - mapf(elev, aht[0], aht[1], aht[2]))*h / 1E3;

	wet_me = mapf(elev, aw[0], aw[1], aw[2]);        // Trop wet mapping function
	dry_me = mapf(elev, ah[0], ah[1], ah[2]) + dm;   // Trop dry mapping function
}

static double ah_mean[] =
{
	+1.2517e+02, +8.503e-01, +6.936e-02, -6.760e+00, +1.771e-01,
	+1.130e-02, +5.963e-01, +1.808e-02, +2.801e-03, -1.414e-03,
	-1.212e+00, +9.300e-02, +3.683e-03, +1.095e-03, +4.671e-05,
	+3.959e-01, -3.867e-02, +5.413e-03, -5.289e-04, +3.229e-04,
	+2.067e-05, +3.000e-01, +2.031e-02, +5.900e-03, +4.573e-04,
	-7.619e-05, +2.327e-06, +3.845e-06, +1.182e-01, +1.158e-02,
	+5.445e-03, +6.219e-05, +4.204e-06, -2.093e-06, +1.540e-07,
	-4.280e-08, -4.751e-01, -3.490e-02, +1.758e-03, +4.019e-04,
	-2.799e-06, -1.287e-06, +5.468e-07, +7.580e-08, -6.300e-09,
	-1.160e-01, +8.301e-03, +8.771e-04, +9.955e-05, -1.718e-06,
	-2.012e-06, +1.170e-08, +1.790e-08, -1.300e-09, +1.000e-10
};

static double bh_mean[] =
{
	+0.000e+00, +0.000e+00, +3.249e-02, +0.000e+00, +3.324e-02,
	+1.850e-02, +0.000e+00, -1.115e-01, +2.519e-02, +4.923e-03,
	+0.000e+00, +2.737e-02, +1.595e-02, -7.332e-04, +1.933e-04,
	+0.000e+00, -4.796e-02, +6.381e-03, -1.599e-04, -3.685e-04,
	+1.815e-05, +0.000e+00, +7.033e-02, +2.426e-03, -1.111e-03,
	-1.357e-04, -7.828e-06, +2.547e-06, +0.000e+00, +5.779e-03,
	+3.133e-03, -5.312e-04, -2.028e-05, +2.323e-07, -9.100e-08,
	-1.650e-08, +0.000e+00, +3.688e-02, -8.638e-04, -8.514e-05,
	-2.828e-05, +5.403e-07, +4.390e-07, +1.350e-08, +1.800e-09,
	+0.000e+00, -2.736e-02, -2.977e-04, +8.113e-05, +2.329e-07,
	+8.451e-07, +4.490e-08, -8.100e-09, -1.500e-09, +2.000e-10
};

static double ah_amp[] =
{
	-2.738e-01, -2.837e+00, +1.298e-02, -3.588e-01, +2.413e-02,
	+3.427e-02, -7.624e-01, +7.272e-02, +2.160e-02, -3.385e-03,
	+4.424e-01, +3.722e-02, +2.195e-02, -1.503e-03, +2.426e-04,
	+3.013e-01, +5.762e-02, +1.019e-02, -4.476e-04, +6.790e-05,
	+3.227e-05, +3.123e-01, -3.535e-02, +4.840e-03, +3.025e-06,
	-4.363e-05, +2.854e-07, -1.286e-06, -6.725e-01, -3.730e-02,
	+8.964e-04, +1.399e-04, -3.990e-06, +7.431e-06, -2.796e-07,
	-1.601e-07, +4.068e-02, -1.352e-02, +7.282e-04, +9.594e-05,
	+2.070e-06, -9.620e-08, -2.742e-07, -6.370e-08, -6.300e-09,
	+8.625e-02, -5.971e-03, +4.705e-04, +2.335e-05, +4.226e-06,
	+2.475e-07, -8.850e-08, -3.600e-08, -2.900e-09, +0.000e+00
};

static double bh_amp[] =
{
	+0.000e+00, +0.000e+00, -1.136e-01, +0.000e+00, -1.868e-01,
	-1.399e-02, +0.000e+00, -1.043e-01, +1.175e-02, -2.240e-03,
	+0.000e+00, -3.222e-02, +1.333e-02, -2.647e-03, -2.316e-05,
	+0.000e+00, +5.339e-02, +1.107e-02, -3.116e-03, -1.079e-04,
	-1.299e-05, +0.000e+00, +4.861e-03, +8.891e-03, -6.448e-04,
	-1.279e-05, +6.358e-06, -1.417e-07, +0.000e+00, +3.041e-02,
	+1.150e-03, -8.743e-04, -2.781e-05, +6.367e-07, -1.140e-08,
	-4.200e-08, +0.000e+00, -2.982e-02, -3.000e-03, +1.394e-05,
	-3.290e-05, -1.705e-07, +7.440e-08, +2.720e-08, -6.600e-09,
	+0.000e+00, +1.236e-02, -9.981e-04, -3.792e-05, -1.355e-05,
	+1.162e-06, -1.789e-07, +1.470e-08, -2.400e-09, -4.000e-10
};

static double aw_mean[] =
{
	+5.640e+01, +1.555e+00, -1.011e+00, -3.975e+00, +3.171e-02,
	+1.065e-01, +6.175e-01, +1.376e-01, +4.229e-02, +3.028e-03,
	+1.688e+00, -1.692e-01, +5.478e-02, +2.473e-02, +6.059e-04,
	+2.278e+00, +6.614e-03, -3.505e-04, -6.697e-03, +8.402e-04,
	+7.033e-04, -3.236e+00, +2.184e-01, -4.611e-02, -1.613e-02,
	-1.604e-03, +5.420e-05, +7.922e-05, -2.711e-01, -4.406e-01,
	-3.376e-02, -2.801e-03, -4.090e-04, -2.056e-05, +6.894e-06,
	+2.317e-06, +1.941e+00, -2.562e-01, +1.598e-02, +5.449e-03,
	+3.544e-04, +1.148e-05, +7.503e-06, -5.667e-07, -3.660e-08,
	+8.683e-01, -5.931e-02, -1.864e-03, -1.277e-04, +2.029e-04,
	+1.269e-05, +1.629e-06, +9.660e-08, -1.015e-07, -5.000e-10
};

static double bw_mean[] =
{
	+0.000e+00, +0.000e+00, +2.592e-01, +0.000e+00, +2.974e-02,
	-5.471e-01, +0.000e+00, -5.926e-01, -1.030e-01, -1.567e-02,
	+0.000e+00, +1.710e-01, +9.025e-02, +2.689e-02, +2.243e-03,
	+0.000e+00, +3.439e-01, +2.402e-02, +5.410e-03, +1.601e-03,
	+9.669e-05, +0.000e+00, +9.502e-02, -3.063e-02, -1.055e-03,
	-1.067e-04, -1.130e-04, +2.124e-05, +0.000e+00, -3.129e-01,
	+8.463e-03, +2.253e-04, +7.413e-05, -9.376e-05, -1.606e-06,
	+2.060e-06, +0.000e+00, +2.739e-01, +1.167e-03, -2.246e-05,
	-1.287e-04, -2.438e-05, -7.561e-07, +1.158e-06, +4.950e-08,
	+0.000e+00, -1.344e-01, +5.342e-03, +3.775e-04, -6.756e-05,
	-1.686e-06, -1.184e-06, +2.768e-07, +2.730e-08, +5.700e-09
};

static double aw_amp[] =
{
	+1.023e-01, -2.695e+00, +3.417e-01, -1.405e-01, +3.175e-01,
	+2.116e-01, +3.536e+00, -1.505e-01, -1.660e-02, +2.967e-02,
	+3.819e-01, -1.695e-01, -7.444e-02, +7.409e-03, -6.262e-03,
	-1.836e+00, -1.759e-02, -6.256e-02, -2.371e-03, +7.947e-04,
	+1.501e-04, -8.603e-01, -1.360e-01, -3.629e-02, -3.706e-03,
	-2.976e-04, +1.857e-05, +3.021e-05, +2.248e+00, -1.178e-01,
	+1.255e-02, +1.134e-03, -2.161e-04, -5.817e-06, +8.836e-07,
	-1.769e-07, +7.313e-01, -1.188e-01, +1.145e-02, +1.011e-03,
	+1.083e-04, +2.570e-06, -2.140e-06, -5.710e-08, +2.000e-08,
	-1.632e+00, -6.948e-03, -3.893e-03, +8.592e-04, +7.577e-05,
	+4.539e-06, -3.852e-07, -2.213e-07, -1.370e-08, +5.800e-09
};

static double bw_amp[] =
{
	+0.000e+00, +0.000e+00, -8.865e-02, +0.000e+00, -4.309e-01,
	+6.340e-02, +0.000e+00, +1.162e-01, +6.176e-02, -4.234e-03,
	+0.000e+00, +2.530e-01, +4.017e-02, -6.204e-03, +4.977e-03,
	+0.000e+00, -1.737e-01, -5.638e-03, +1.488e-04, +4.857e-04,
	-1.809e-04, +0.000e+00, -1.514e-01, -1.685e-02, +5.333e-03,
	-7.611e-05, +2.394e-05, +8.195e-06, +0.000e+00, +9.326e-02,
	-1.275e-02, -3.071e-04, +5.374e-05, -3.391e-05, -7.436e-06,
	+6.747e-07, +0.000e+00, -8.637e-02, -3.807e-03, -6.833e-04,
	-3.861e-05, -2.268e-05, +1.454e-06, +3.860e-07, -1.068e-07,
	+0.000e+00, -2.658e-02, -1.947e-03, +7.131e-04, -3.506e-05,
	+1.885e-07, +5.792e-07, +3.990e-08, +2.000e-08, -5.700e-09
};

/// GMF mapping function
void GMF(double mjd, double lat, double lon, double h, double elev, double& dry_me, double& wet_me)
{
	double pi = 4 * atan(1.0);
	double zenith = pi / 2.0 - elev;

	double doy = mjd - 44239.0 + 1 - 28;

	double x = cos(lat)*cos(lon);
	double y = cos(lat)*sin(lon);
	double z = sin(lat);

	int nmax = 9;

	double v[20][20] = { 0.0 };
	double w[20][20] = { 0.0 };

	v[0][0] = 1.0;
	w[0][0] = 0.0;
	v[1][0] = z*v[0][0];
	w[1][0] = 0.0;

	int i = 0, j = 0;

	for (i = 2; i <= nmax; ++i)
	{
		v[i][0] = ((2 * i - 1)*z*v[i - 1][0] - (i - 1)*v[i - 2][0]) / i;
		w[i][0] = 0.0;
	}

	for (j = 1; j <= nmax; ++j)
	{
		v[j][j] = (2 * j - 1)*(x*v[j - 1][j - 1] - y*w[j - 1][j - 1]);
		w[j][j] = (2 * j - 1)*(x*w[j - 1][j - 1] + y*v[j - 1][j - 1]);
		if (j<nmax)
		{
			v[j + 1][j] = (2 * j + 1)*z*v[j][j];
			w[j + 1][j] = (2 * j + 1)*z*w[j][j];
		}

		for (i = j + 2; i <= nmax; ++i)
		{
			v[i][j] = ((2 * i - 1)*z*v[i - 1][j] - (i + j - 1)*v[i - 2][j]) / (i - j);
			w[i][j] = ((2 * i - 1)*z*w[i - 1][j] - (i + j - 1)*w[i - 2][j]) / (i - j);
		}
	}

	double bh = 0.0029;
	double c0h = 0.062;
	double phh = 0.0, c11h = 0.0, c10h = 0.0;

	if (lat<0)
	{
		phh = pi;
		c11h = 0.007;
		c10h = 0.002;
	}
	else
	{
		phh = 0;
		c11h = 0.005;
		c10h = 0.001;
	}

	double ch = c0h + ((cos(doy / 365.250 * 2 * pi + phh) + 1)*c11h / 2 + c10h)*(1 - cos(lat));

	double ahm = 0.0, aha = 0.0;
	int k = 0;

	for (i = 0; i <= nmax; ++i)
	{
		for (j = 0; j <= i; ++j)
		{
			k = k + 1;
			ahm = ahm + (ah_mean[k - 1] * v[i][j] + bh_mean[k - 1] * w[i][j]);
			aha = aha + (ah_amp[k - 1] * v[i][j] + bh_amp[k - 1] * w[i][j]);
		}
	}

	double ah = (ahm + aha*cos(doy / 365.25*2.0*pi))*1e-5;

	double sine = sin(pi / 2 - zenith);
	double cose = cos(pi / 2 - zenith);
	double beta = bh / (sine + ch);
	double gamma = ah / (sine + beta);
	double topcon = (1.0 + ah / (1.0 + bh / (1.0 + ch)));

	dry_me = topcon / (sine + gamma);

	//height correction for hydrostatic mapping function from Niell (1996)
	double a_ht = 2.53e-5;
	double b_ht = 5.49e-3;
	double c_ht = 1.14e-3;
	double hs_km = h / 1000.0;

	beta = b_ht / (sine + c_ht);
	gamma = a_ht / (sine + beta);
	topcon = (1.0 + a_ht / (1.0 + b_ht / (1.0 + c_ht)));
	double ht_corr_coef = 1 / sine - topcon / (sine + gamma);
	double ht_corr = ht_corr_coef*hs_km;

	dry_me = (dry_me)+ht_corr;

	double bw = 0.00146;
	double cw = 0.04391;

	double awm = 0.0;
	double awa = 0.0;

	k = 0;
	for (i = 0; i <= nmax; ++i)
	{
		for (j = 0; j <= i; ++j)
		{
			k = k + 1;
			awm = awm + (aw_mean[k - 1] * v[i][j] + bw_mean[k - 1] * w[i][j]);
			awa = awa + (aw_amp[k - 1] * v[i][j] + bw_amp[k - 1] * w[i][j]);
		}
	}

	double aw = (awm + awa*cos(doy / 365.25 * 2 * pi))*1e-5;

	beta = bw / (sine + cw);
	gamma = aw / (sine + beta);
	topcon = (1.0 + aw / (1.0 + bw / (1.0 + cw)));

	wet_me = topcon / (sine + gamma);
}

void VMF(double ah, double aw, double E, double Lat, double H, int doy, double &dry_me, double &wet_me)
{
	double bd = 0.0029, cd = 0;//干量投影系数
	double bw = 0.00146, cw = 0.04391;
	double aht = 2.53e-5, bht = 5.49e-3, cht = 1.14e-3;
	if (Lat < 0)
	{
		double c0 = 0.062, c11 = 0.001, c10 = 0.006, ph = PI;//PI
		cd = c0 + ((cos(2 * PI*(doy - 28) / 365 + ph) + 1)*c11 / 2 + c10)*(1 - cos(-Lat));
	}
	else
	{
		double c0 = 0.062, c11 = 0.0, c10 = 0.006, ph = 0;
		cd = c0 + ((cos(2 * PI*(doy - 28) / 365 + ph) + 1)*c11 / 2 + c10)*(1 - cos(Lat));
	}
	//计算投影函数
	dry_me = (1 + ah / (1 + bd / (1 + cd))) / (sin(E) + ah / (sin(E) + bd / (sin(E) + cd)))
		+ (1 / sin(E) - (1 + aht / (1 + bht / (1 + cht))) / (sin(E) + (aht / (sin(E) + bht / (sin(E) + cht)))))*H / 1000;
	wet_me = (1 + aw / (1 + bw / (1 + cw))) / (sin(E) + aw / (sin(E) + bw / (sin(E) + cw)));

}

/// Hopfield Model ------------------------------------------------------------------------------
/* Parameter:
*            double     h             I         height
*            double     el            I         satellite elevation
*            double     dry_delay     O         trop dry delay
*            double     wet_delay     O         trop wet delay
*            int        flag          I         0: zenith delay    1: slant delay(line of sight)
*  return  total delay(line of sight)
* ---------------------------------------------------------------------------------------------- */
double Trop_Hopfield(double h, double elev, double& dry_delay, double& wet_delay, int flag)
{
	if (fabs(h)>11000)
	{
		dry_delay = 0.0;
		wet_delay = 0.0;
		return 0.0;
	}
	double pi = 4 * atan(1.0);

	double p0 = 1013.25;
	double t0 = 18.0 + 273.16;
	double rh0 = 0.5;
	double h0 = 0.0;

	double p = p0*pow((1 - 0.0000226*(h - h0)), 5.225);
	double t = t0 - 0.0065*(h - h0);
	double rh = rh0*exp(-0.0006396*(h - h0));

	double es = rh*exp(-37.2465 + 0.213166*t - 0.000256908*t*t);

	double hw = 11000.0;
	double hd = 40136.0 + 148.72*(t0 - 273.16);
	double kw = 0.0;
	double me0 = rh0*exp(-37.2465 + 0.213166*t0 - 0.000256908*t0*t0);
	if (h<11000.0)
	{
		kw = (7.46512e-2)*me0 / (t0*t0)*(pow(hw - h, 5) / pow(hw, 4));
	}
	double kd = (155.2e-7)*p0 / t0*(pow(hd - h, 5) / pow(hd, 4));
	double EE = elev*180.0 / pi;

	double temp1 = sqrt(EE*EE + 6.25)*pi / 180.0;
	double temp2 = sqrt(EE*EE + 2.25)*pi / 180.0;

	dry_delay = kd;
	wet_delay = kw;

	double delay = kd / sin(temp1) + kw / sin(temp2);

	if (flag == 1)
	{
		dry_delay = (dry_delay) / sin(temp1);
		wet_delay = (wet_delay) / sin(temp2);
	}

	return delay;
}

/// Saastamoinen Model
// el-unit is radian. h-unit is m.
// flag=0,dry_delay and wet_delay are zenith delays.
// flag=1,dry_delay and wet_delay are transfer route delays.
double Trop_Saastamoinen(double h, double elev, double& dry_delay, double& wet_delay, int flag)
{
	double p0 = 1013.25;
	double t0 = 18.0 + 273.16;
	double rh0 = 0.5;
	double h0 = 0.0;

	if (fabs(h)>11000)
	{
		dry_delay = 0.0;
		wet_delay = 0.0;
		return 0.0;
	}

	double p = p0*pow((1 - 0.0000226*(h - h0)), 5.225);
	double t = t0 - 0.0065*(h - h0);
	double rh = rh0*exp(-0.0006396*(h - h0));

	double es = rh*exp(-37.2465 + 0.213166*t - 0.000256908*t*t);

	double de = 0.0, e = 0.0;
	double hh = h / 1000.0;
	if (fabs(elev) < 1.0e-7)
		de = 0.0;
	else
		de = 16.0 / 206265 / t*(p + 4810 / t*es) / tan(elev);

	e = elev + de;
	double a = 1.16 - 0.15e-3*hh + 0.716e-8*h*hh;
	double b = 0.0;
	if (fabs(fabs(elev) - 3.141592653 / 2.0)>1.0e-9)
		b = a / tan(elev) / tan(elev);

	dry_delay = 0.002277*p;
	wet_delay = 0.002277*(1255.0 / t + 0.05)*es;

	double delay = ((dry_delay)+(wet_delay)-0.002277*b) / sin(elev);

	if (flag == 1)
	{
		dry_delay = (dry_delay) / sin(elev);
		wet_delay = (wet_delay) / sin(elev);
	}

	return delay;
}

double Trop_UNB3(double blh[3], double doy, double elev, double& dry_delay, double& wet_delay, int flag)
{
	double LATRAD = blh[0];
	double HEIGHTM = blh[2];

	MatrixT AVG(5, 6);

	AVG.m_mat << 15.0, 1013.25, 299.65, 75.00, 6.30, 2.77,
		30.0, 1017.25, 294.15, 80.00, 6.05, 3.15,
		45.0, 1015.75, 283.15, 76.00, 5.58, 2.57,
		60.0, 1011.75, 272.15, 77.50, 5.39, 1.81,
		75.0, 1013.00, 263.65, 82.50, 4.53, 1.55;
	MatrixT AMP(5, 6);
	AMP.m_mat << 15.0, 0.00, 0.00, 0.00, 0.00, 0.00,
		30.0, -3.75, 7.00, 0.00, 0.25, 0.33,
		45.0, -2.25, 11.00, -1.00, 0.32, 0.46,
		60.0, -1.75, 15.00, -2.50, 0.81, 0.74,
		75.0, -0.50, 14.50, 2.50, 0.62, 0.30;

	double EXCEN2 = 6.6943799901413e-03;
	double MD = 28.9644;
	double MW = 18.0152;
	double K1 = 77.604;
	double K2 = 64.79;
	double K3 = 3.776e5;
	double R = 8314.34;
	double C1 = 2.2768e-03;
	double K2PRIM = K2 - K1*(MW / MD);
	double RD = R / MD;

	double DOY2RAD = (0.31415926535897935601e01) * 2 / 365.25;
	MatrixT ABC_AVG(5, 4);
	ABC_AVG.m_mat << 15.0, 1.2769934e-3, 2.9153695e-3, 62.610505e-3,
		30.0, 1.2683230e-3, 2.9152299e-3, 62.837393e-3,
		45.0, 1.2465397e-3, 2.9288445e-3, 63.721774e-3,
		60.0, 1.2196049e-3, 2.9022565e-3, 63.824265e-3,
		75.0, 1.2045996e-3, 2.9024912e-3, 64.258455e-3;
	MatrixT ABC_AMP(5, 4);
	ABC_AMP.m_mat << 15.0, 0.0, 0.0, 0.0,
		30.0, 1.2709626e-5, 2.1414979e-5, 9.0128400e-5,
		45.0, 2.6523662e-5, 3.0160779e-5, 4.3497037e-5,
		60.0, 3.4000452e-5, 7.2562722e-5, 84.795348e-5,
		75.0, 4.1202191e-5, 11.723375e-5, 170.37206e-5;

	double A_HT = 2.53e-5;
	double B_HT = 5.49e-3;
	double C_HT = 1.14e-3;
	double HT_TOPCON = 1 + A_HT / (1 + B_HT / (1 + C_HT));

	MatrixT ABC_W2P0(5, 4);
	ABC_W2P0.m_mat << 15.0, 5.8021897e-4, 1.4275268e-3, 4.3472961e-2,
		30.0, 5.6794847e-4, 1.5138625e-3, 4.6729510e-2,
		45.0, 5.8118019e-4, 1.4572752e-3, 4.3908931e-2,
		60.0, 5.9727542e-4, 1.5007428e-3, 4.4626982e-2,
		75.0, 6.1641693e-4, 1.7599082e-3, 5.4736038e-2;

	double LATDEG = LATRAD * 180.0 / PI;
	double  TD_O_Y = doy;
	if (LATDEG<0)
	{
		TD_O_Y = TD_O_Y + 182.625;
	}
	double COSPHS = cos((TD_O_Y - 28) * DOY2RAD);
	double LAT = abs(LATDEG);
	double P1 = 1, P2 = 1, M = 0;
	if (LAT >= 75)
	{
		P1 = 5;
		P2 = 5;
		M = 0;
	}
	else if (LAT <= 15)
	{
		P1 = 1;
		P2 = 1;
		M = 0;
	}
	else
	{
		P1 = int((LAT - 15) / 15) + 1;
		P2 = P1 + 1;
		double aa = LAT - AVG(P1, 1);
		double bb = AVG(P2, 1) - AVG(P1, 1);
		M = (aa) / (bb);
	}

	double PAVG = M * (AVG(P2, 2) - AVG(P1, 2)) + AVG(P1, 2);
	double TAVG = M * (AVG(P2, 3) - AVG(P1, 3)) + AVG(P1, 3);
	double EAVG = M * (AVG(P2, 4) - AVG(P1, 4)) + AVG(P1, 4);
	double BETAAVG = M * (AVG(P2, 5) - AVG(P1, 5)) + AVG(P1, 5);
	double LAMBDAAVG = M * (AVG(P2, 6) - AVG(P1, 6)) + AVG(P1, 6);

	double PAMP = M * (AMP(P2, 2) - AMP(P1, 2)) + AMP(P1, 2);
	double TAMP = M * (AMP(P2, 3) - AMP(P1, 3)) + AMP(P1, 3);
	double EAMP = M * (AMP(P2, 4) - AMP(P1, 4)) + AMP(P1, 4);
	double BETAAMP = M * (AMP(P2, 5) - AMP(P1, 5)) + AMP(P1, 5);
	double LAMBDAAMP = M * (AMP(P2, 6) - AMP(P1, 6)) + AMP(P1, 6);

	double P0 = PAVG - PAMP * COSPHS;
	double T0 = TAVG - TAMP * COSPHS;
	double E0 = EAVG - EAMP * COSPHS;
	double BETA = BETAAVG - BETAAMP * COSPHS;
	BETA = BETA / 1000;
	double LAMBDA = LAMBDAAVG - LAMBDAAMP * COSPHS;

	double ES = 0.01 * exp(1.2378847e-5 * (pow(T0, 2)) - 1.9121316e-2 * T0 + 3.393711047e1 - 6.3431645e3 * (pow(T0, -1)));
	double FW = 1.00062 + 3.14e-6 * P0 + 5.6e-7 * (pow((T0 - 273.15), 2));
	E0 = (E0 / 1.00e2) * ES * FW;

	double EP = 9.80665 / 287.054 / BETA;

	double T = T0 - BETA * HEIGHTM;
	double P = P0 * pow((T / T0), EP);
	double E = E0 *pow((T / T0), (EP * (LAMBDA + 1)));

	double GEOLAT = atan((1.0 - EXCEN2)*tan(LATRAD));
	double DGREF = 1.0 - 2.66e-03*cos(2.0*GEOLAT) - 2.8e-07*HEIGHTM;
	double GM = 9.784 * DGREF;
	double DEN = (LAMBDA + 1.0) * GM;

	double TM = T * (1 - BETA * RD / DEN);

	double HZD = C1 / DGREF * P;

	double WZD = 1.0e-6 * (K2PRIM + K3 / TM) * RD * E / DEN;

	double A_AVG = M * (ABC_AVG(P2, 2) - ABC_AVG(P1, 2)) + ABC_AVG(P1, 2);
	double B_AVG = M * (ABC_AVG(P2, 3) - ABC_AVG(P1, 3)) + ABC_AVG(P1, 3);
	double C_AVG = M * (ABC_AVG(P2, 4) - ABC_AVG(P1, 4)) + ABC_AVG(P1, 4);

	double A_AMP = M * (ABC_AMP(P2, 2) - ABC_AMP(P1, 2)) + ABC_AMP(P1, 2);
	double B_AMP = M * (ABC_AMP(P2, 3) - ABC_AMP(P1, 3)) + ABC_AMP(P1, 3);
	double C_AMP = M * (ABC_AMP(P2, 4) - ABC_AMP(P1, 4)) + ABC_AMP(P1, 4);

	double A = A_AVG - A_AMP * COSPHS;
	double B = B_AVG - B_AMP * COSPHS;
	double C = C_AVG - C_AMP * COSPHS;

	double SINE = sin(elev);

	double ALPHA = B / (SINE + C);
	double GAMMA = A / (SINE + ALPHA);
	double TOPCON = (1 + A / (1 + B / (1 + C)));
	double HMF = TOPCON / (SINE + GAMMA);

	ALPHA = B_HT / (SINE + C_HT);
	GAMMA = A_HT / (SINE + ALPHA);
	double HT_CORR_COEF = 1 / SINE - HT_TOPCON / (SINE + GAMMA);
	double HT_CORR = HT_CORR_COEF * HEIGHTM / 1000;
	HMF = HMF + HT_CORR;

	A = M * (ABC_W2P0(P2, 2) - ABC_W2P0(P1, 2)) + ABC_W2P0(P1, 2);
	B = M * (ABC_W2P0(P2, 3) - ABC_W2P0(P1, 3)) + ABC_W2P0(P1, 3);
	C = M * (ABC_W2P0(P2, 4) - ABC_W2P0(P1, 4)) + ABC_W2P0(P1, 4);

	ALPHA = B / (SINE + C);
	GAMMA = A / (SINE + ALPHA);
	TOPCON = (1 + A / (1 + B / (1 + C)));
	double WMF = TOPCON / (SINE + GAMMA);
	double RTROP = HZD*HMF + WZD*WMF;

	dry_delay = HZD;
	wet_delay = WZD;

	return (HZD + WZD);
}

// latitude(radian)
double Trop_NMF_Saastamoinen(int doy, double latitude, double h, double elev, double& dry_delay, double& wet_delay, int flag)
{
	double p0 = 1013.25;
	double t0 = 18.0 + 273.16;
	double rh0 = 0.5;
	double h0 = 0.0;

	if (h > 15000.0)
		h = 15000.0;
	else if (h < 0)
		h = 0;

	double p = p0*pow((1 - 0.0000226*(h - h0)), 5.225);
	double t = t0 - 0.0065*(h - h0);
	double rh = rh0*exp(-0.0006396*(h - h0));

	double es = rh*exp(-37.2465 + 0.213166*t - 0.000256908*t*t);

	dry_delay = 0.002277*p;
	wet_delay = 0.002277*es*(1255.0 / t + 0.05);

	double dry_me = 1.0, wet_me = 1.0;

	NMF(doy, latitude, h, elev, dry_me, wet_me);

	if (flag == 1)
	{
		(dry_delay) = (dry_delay)*dry_me;
		(wet_delay) = (wet_delay)*wet_me;
	}

	return (dry_delay)*dry_me + (wet_delay)*wet_me;
}

// latitude(radian)
double Trop_GMF_Saastamoinen( double mjd, double latitude,double lontitude,double h,double elev, double& dry_delay, double& wet_delay, int flag)
{
	double p0 = 1013.25;
	double t0 = 18.0 + 273.16;
	double rh0 = 0.5;
	double h0 = 0.0;

	if (h > 15000.0)
		h = 15000.0;

	double p = p0*pow((1 - 0.0000226*(h - h0)), 5.225);
	double t = t0 - 0.0065*(h - h0);
	double rh = rh0*exp(-0.0006396*(h - h0));

	double es = rh*exp(-37.2465 + 0.213166*t - 0.000256908*t*t);

	dry_delay = 0.002277*p;
	wet_delay = 0.002277*es*(1255.0 / t + 0.05);

	double dry_me = 1.0, wet_me = 1.0;

	GMF(mjd, latitude, lontitude, h, elev, dry_me, wet_me);

	if (flag == 1)
	{
		(dry_delay) = (dry_delay)*dry_me;
		(wet_delay) = (wet_delay)*wet_me;
	}

	return (dry_delay)*dry_me + (wet_delay)*wet_me;
}
