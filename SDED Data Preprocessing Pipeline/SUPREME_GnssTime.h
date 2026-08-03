/* -------------------------------------------------------------------------
* SUPREME_GnssTime.h : GNSS Time in Multi-GNSS Solutions
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2018.12.07
* ------------------------------------------------------------------------- */
#ifndef SUPREME_GNSS_TIME_H_HH
#define SUPREME_GNSS_TIME_H_HH

#include <stdio.h>
#include <stdlib.h>

#define MAX_LEAPS    64                                 /* Max number of leap seconds table   */

const static double gpst0[] = { 1980, 1, 6, 0, 0, 0 };  /* GPS start time */
const static double gst0[] = { 1999, 8, 22, 0, 0, 0 };  /* GAL start time */
const static double bdt0[] = { 2006, 1, 1, 0, 0, 0 };   /* BDS start time */

static double leaps[MAX_LEAPS + 1][7] =                 /* Leap Seconds (Y,M,D,h,m,s,utc-gpst) */
{
	{ 2017, 1, 1, 0, 0, 0, -18 },
	{ 2015, 7, 1, 0, 0, 0, -17 },
	{ 2012, 7, 1, 0, 0, 0, -16 },
	{ 2009, 1, 1, 0, 0, 0, -15 },
	{ 2006, 1, 1, 0, 0, 0, -14 },
	{ 1999, 1, 1, 0, 0, 0, -13 },
	{ 1997, 7, 1, 0, 0, 0, -12 },
	{ 1996, 1, 1, 0, 0, 0, -11 },
	{ 1994, 7, 1, 0, 0, 0, -10 },
	{ 1993, 7, 1, 0, 0, 0, -9 },
	{ 1992, 7, 1, 0, 0, 0, -8 },
	{ 1991, 1, 1, 0, 0, 0, -7 },
	{ 1990, 1, 1, 0, 0, 0, -6 },
	{ 1988, 1, 1, 0, 0, 0, -5 },
	{ 1985, 7, 1, 0, 0, 0, -4 },
	{ 1983, 7, 1, 0, 0, 0, -3 },
	{ 1982, 7, 1, 0, 0, 0, -2 },
	{ 1981, 7, 1, 0, 0, 0, -1 },
	{ 0 }
};

enum GNSS_TimeYype{ NONET, UTC, GPST, BDST, GALT, GLOT, TimeTypeNUM };                         // GNSS Time Type

typedef struct
{
	time_t time;        /* time (s) expressed by standard time_t */
	double sec;         /* fraction of second under 1 s          */
} gtime_t;

typedef struct {        /* second of day                         */
	long   sn;
	double tos;
} sod_t;

typedef struct {        /* modified Julian date                  */
	long  day;
	sod_t ds;
} mjd_t;

gtime_t timeadd(gtime_t t, double sec);
double  timediff(gtime_t t1, gtime_t t2);

gtime_t epoch2time(const double ep[6]);
void    time2epoch(gtime_t t, double ep[6]);

gtime_t gpst2time(int week, double sec);
double  time2gpst(gtime_t t, int *week);

gtime_t gst2time(int week, double sec);
double  time2gst(gtime_t t, int *week);

gtime_t bdt2time(int week, double sec);
double time2bdt(gtime_t t, int *week);

gtime_t gpst2utc(gtime_t t);
gtime_t utc2gpst(gtime_t t);

gtime_t gpst2bdt(gtime_t t);
gtime_t bdt2gpst(gtime_t t);

double  utc2gmst(gtime_t t, double ut1_utc);

double mjd2jd(double mjd);
double jd2mjd(double jd);

double  time2doy(gtime_t t);

int tow2dow(double secofweek);

gtime_t timeget(void);
void    timeset(gtime_t t);

int LongYear(int yy, int mm);


/// SUPREME GNSS Time Class
class gnsstime
{
public:
	int m_Type, m_Year, m_Month, m_Day, m_Hour, m_Min, m_Doy, m_Week, m_Dow;
	double m_Sec, m_SecondofWeek, m_Jd, m_Mjd;
	double m_Epoch[6];

	gtime_t m_gtime;

	gnsstime();
	gnsstime(double jd);
	gnsstime(double ep[6]);
	gnsstime(int ww, double ss);
	gnsstime(int YY, int MM, int DD, int hh, int mm, double ss);

	void Init();

	double _jd2mjd();
	double _mjd2jd();
	double _date2jd();
	double _date2mjd();

	int _GetDoy();
	int _date2utc();
	int _date2gpst();
	int _gpst2date();
	int _date2bdst();
	int _bdst2date();
	int _date2ep();
	int _ep2date();
	int _jd2date();
	int _mjd2date();

	int Utc2Gpst();
	int Gpst2Utc();
	int Utc2Bdst();
	int Bdst2Utc();
	int Gpst2Bdst();
	int Bdst2Gpst();

	void _GetLocalTime();
	void _GetLocalTime_GPST();
	void _GetLocalTime_GPST(double offset);

	void _GetDaySecond();
	const int _GetDay() const{ return _day; }
	const double _GetSecond() const{ return _second; }

	int _LongYear();
	char *_Time2Str();

	friend double operator-(const gnsstime &T1, const gnsstime &T2);
	friend gnsstime operator+(const gnsstime &T, const double &sec);

private:
	int _day;
	double _second;
	char _TimeStr[30];
	//gtime_t m_gtime;
};

inline double operator-(const gnsstime &T1, const gnsstime &T2);
gnsstime operator+(const gnsstime &T, const double &sec);

/// Operator Overloading for GNSS Time
bool operator>(const gnsstime &T1, const gnsstime &T2);
bool operator<(const gnsstime &T1, const gnsstime &T2);
bool operator>=(const gnsstime &T1, const gnsstime &T2);
bool operator<=(const gnsstime &T1, const gnsstime &T2);
bool operator==(const gnsstime &T1, const gnsstime &T2);

#endif
