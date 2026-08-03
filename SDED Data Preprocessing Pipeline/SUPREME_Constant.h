/* -------------------------------------------------------------------------
* SUPREME_Constant.h : Constant in GNSS
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2018.12.07
* ------------------------------------------------------------------------- */
#ifndef SUPREME_CONSTANT_H_HH
#define SUPREME_CONSTANT_H_HH

#define VERSION	     7.0

#define LEAPSECOND   18                                                   // Current leap second

#define LIGHTSPEED   299792458.0                                          // Light speed in vacuum (m/sec)
#define DCB_Parameter         1E-9 * LIGHTSPEED
#define MAXFREQ      7                                                    // Maximun number of frequency

#define GPS_FREQU_L1          1575420000.0                                // Hz               
#define GPS_FREQU_L2          1227600000.0                                // Hz 
#define GPS_FREQU_L5          1176450000.0								  // Hz

#define G_B_12    -(SQR(GPS_FREQU_L2)/(SQR(GPS_FREQU_L1)-SQR(GPS_FREQU_L2)))
#define G_B_13    -(SQR(GPS_FREQU_L5)/(SQR(GPS_FREQU_L1)-SQR(GPS_FREQU_L5)))

#define GPS_WAVELENGTH_L1     (LIGHTSPEED / GPS_FREQU_L1)                 // m 
#define GPS_WAVELENGTH_L2     (LIGHTSPEED / GPS_FREQU_L2)                 // m
#define GPS_WAVELENGTH_L5     (LIGHTSPEED / GPS_FREQU_L5)                 // m
#define BDS_L1_L2_gamma		 (SQR(GPS_WAVELENGTH_L1))/(SQR(GPS_WAVELENGTH_L2)-SQR(GPS_WAVELENGTH_L1))
#define BDS_FREQU_B1          1561098000.0                                // Hz
#define BDS_FREQU_B2          1207140000.0                                // Hz
#define BDS_FREQU_B3          1268520000.0                                // Hz
#define BDS_FREQU_B1C         1575420000.0                                // Hz
#define BDS_FREQU_B2a         1176450000.0                                // Hz

#define C_B_12    -(SQR(BDS_FREQU_B3)/(SQR(BDS_FREQU_B1)-SQR(BDS_FREQU_B3)))
#define C_B_13    -(SQR(BDS_FREQU_B2)/(SQR(BDS_FREQU_B1)-SQR(BDS_FREQU_B2)))
#define C_B_14    -(SQR(BDS_FREQU_B1C)/(SQR(BDS_FREQU_B1)-SQR(BDS_FREQU_B1C)))
#define C_B_15    -(SQR(BDS_FREQU_B2a)/(SQR(BDS_FREQU_B1)-SQR(BDS_FREQU_B2a)))

#define BDS_FREQU_B2b         1207140000.0                                // Hz
#define BDS_WAVELENGTH_B1     (LIGHTSPEED / BDS_FREQU_B1)                 // m
#define BDS_WAVELENGTH_B2     (LIGHTSPEED / BDS_FREQU_B2)                 // m
#define BDS_WAVELENGTH_B3     (LIGHTSPEED / BDS_FREQU_B3)                 // m

#define BDS_B1_B3_gamma      (SQR(BDS_WAVELENGTH_B1))/(SQR(BDS_WAVELENGTH_B3)-SQR(BDS_WAVELENGTH_B1))
#define BDS_WAVELENGTH_B1C    (LIGHTSPEED / BDS_FREQU_B1C)                // m
#define BDS_WAVELENGTH_B2a    (LIGHTSPEED / BDS_FREQU_B2a)                // m
#define BDS_WAVELENGTH_B2b    (LIGHTSPEED / BDS_FREQU_B2b)                // m

#define GLO_FREQU_L1_BASE     1602000000.0                                // Hz 
#define GLO_FREQU_L2_BASE     1246000000.0                                // Hz 
#define GLO_FREQU_L1a         1600995000.0
#define GLO_FREQU_L2a         1248060000.0
#define GLO_FREQU_L3          1202025000.0
#define GLO_FREQU_L1_STEP     562500.0                                    // Hz 
#define GLO_FREQU_L2_STEP     437500.0                                    // Hz 
#define GLO_FREQU_L1(a)       (GLO_FREQU_L1_BASE+(a)*GLO_FREQU_L1_STEP)
#define GLO_FREQU_L2(a)       (GLO_FREQU_L2_BASE+(a)*GLO_FREQU_L2_STEP)
#define GLO_WAVELENGTH_L1(a)  (LIGHTSPEED / GLO_FREQU_L1(a))              // m
#define GLO_WAVELENGTH_L2(a)  (LIGHTSPEED / GLO_FREQU_L2(a))              // m
#define GLO_WAVELENGTH_L1a    (LIGHTSPEED / GLO_FREQU_L1a)    
#define GLO_WAVELENGTH_L2a    (LIGHTSPEED / GLO_FREQU_L2a) 
#define GLO_WAVELENGTH_L3     (LIGHTSPEED / GLO_FREQU_L3) 

#define GAL_FREQU_E1          1575420000.0                                // HZ
#define GAL_FREQU_E5a         1176450000.0                                // HZ
#define GAL_FREQU_E5b         1207140000.0                                // HZ
#define GAL_FREQU_E5          1191795000.0   
#define GAL_FREQU_E6          1278750000.0  

#define E_B_12    -(SQR(GAL_FREQU_E5a)/(SQR(GAL_FREQU_E1)-SQR(GAL_FREQU_E5a)))
#define E_B_13    -(SQR(GAL_FREQU_E5b)/(SQR(GAL_FREQU_E1)-SQR(GAL_FREQU_E5b)))
#define E_B_14    -(SQR(GAL_FREQU_E5)/(SQR(GAL_FREQU_E1)-SQR(GAL_FREQU_E5)))
#define E_B_15    -(SQR(GAL_FREQU_E6)/(SQR(GAL_FREQU_E1)-SQR(GAL_FREQU_E6)))

#define GAL_WAVELENGTH_E1     (LIGHTSPEED / GAL_FREQU_E1)                 // m
#define GAL_WAVELENGTH_E5a    (LIGHTSPEED / GAL_FREQU_E5a)                // m
#define GAL_WAVELENGTH_E5b    (LIGHTSPEED / GAL_FREQU_E5b)                // m
#define GAL_WAVELENGTH_E5     (LIGHTSPEED / GAL_FREQU_E5) 
#define GAL_WAVELENGTH_E6     (LIGHTSPEED / GAL_FREQU_E6) 

#define PI                    (4*atan(1.0))
#define SC2RAD                3.1415926535898                             // semi-circle to radian (IS-GPS)
#define AU                    149597870691.0                              // 1 AU (m)

// WGS84 Ellipsoid
#define WGS84_A 6378137.0
#define WGS84_E 0.00669437999013
#define WGS84_F 0.003352810664747
// PZ90 Ellipsoid
#define PZ90_A  6378136.0
#define PZ90_E  0.006694366193099
#define PZ90_F  0.003352803743019

#define E_GM        3986004.415e+8                                       // Earth gravitational constant
#define S_GM        1.327124E+20                                         // Sun gravitational constant
#define M_GM        4.902801E+12                                         // Moon gravitational constant
#define E_GM_GPS    3986005.0e+8

#define OMEGA_GPS   7.2921151467e-5                                      // Earth angular velocity (IS-GPS) (rad/s)
#define OMEGA_GLO   7.2921151467e-5
#define OMEGA_BDS   7.292115e-5
#define GLO_J20     1082625.7e-9

#define HION        350000.0                                             // Ionosphere height (m)
#define ION_FACTEC  40.28e16                                             // Ionospheric Factor

enum NavSystem { SYS_G, SYS_R, SYS_E, SYS_C, SYS_Q, SYS_S, SYS_NUM };    // 7 Navigation system
const char SYSTEM_STR[SYS_NUM][5]{"GPS", "GLO", "GAL", "BDS", "QZSS", "SBAS"};
const char SYSTEM_FLAG[SYS_NUM]{'G', 'R', 'E', 'C', 'J', 'S'};

#define SIGEMA_P    0.30                                                 // (m) default sigema of psuedo-range
#define SIGEMA_L    0.003                                                // (m) default sigema of carrier phase

static char CODE_PRIS[7][MAXFREQ][16] = {                                 // Code Priority Table

	/// L1/E1      L2/B1        L5/E5a/L3 L6/LEX/B3 E5b/B2    E5(a+b)  S
	{ "CPYWMNSL", "PYWCMNDSLX", "IQX", "", "", "", "" },                  // GPS
	{ "PC", "PC", "IQX", "", "", "", "" },                                // GLONASS
	{ "CABXZ", "", "IQX", "ABCXZ", "IQX", "IQX", "" },                    // GALILEO
	{ "CSLXZ", "SLX", "IQX", "SLX", "", "", "" },                         // QZSS
	{ "C", "", "IQX", "", "", "", "" },                                   // SBAS
	{ "IQX", "IQX", "IQX", "IQX", "IQX", "", "" },                        // BDS
	{ "", "", "ABCX", "", "", "", "ABCX" }                                // IRNSS
};

class GNSS_CONSTANT
{
	GNSS_CONSTANT();

};

#endif
