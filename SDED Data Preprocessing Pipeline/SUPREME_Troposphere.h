/* -------------------------------------------------------------------------
* SUPREME_Troposphere.h : Troposphere models in SUPREME
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2019.01.08
* ------------------------------------------------------------------------- */
/// Troposphere delay strategy
#define TROP_ESTIMATE              0      // Estimate trop wet delay
#define TROP_GRADIENT              1      // Estimate trop gradient
#define TROP_GRID                  2      // Correct by grid

/// Troposphere model
#define TROP_MODEL_HOPFIELd        1      // Holfield model
#define TROP_MODEL_SAASTAMOINEn    2      // Saastamoinen model
#define TROP_UNB3                  3      // UNB3 model
#define TROP_GPT2W                 4      // GPT2W model
#define TROP_IGG                   5      // IGGTrop Model

/// Mapping function
void NMF(int doy, double latitude, double h, double elev, double& dry_me, double& wet_me);
void GMF(double mjd, double lat, double lon, double h, double elev, double& dry_me, double& wet_me);
void VMF(double ah, double aw, double E, double Lat, double H, int doy, double &md, double &mw);

/// Troposphere model
double Trop_Hopfield(double h, double elev, double& dry_delay, double& wet_delay, int flag);
double Trop_Saastamoinen(double h, double elev, double& dry_delay, double& wet_delay, int flag);
double Trop_UNB3(double blh[3], double doy, double elev, double& dry_delay, double& wet_delay, int flag);

double Trop_NMF_Saastamoinen(int doy, double latitude, double h, double elev, double& dry_delay, double& wet_delay, int flag);
double Trop_GMF_Saastamoinen(double mjd, double latitude, double lontitude, double h, double elev, double& dry_delay, double& wet_delay, int flag);
