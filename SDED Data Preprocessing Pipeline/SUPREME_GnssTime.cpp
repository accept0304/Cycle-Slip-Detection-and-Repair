
#include <time.h>
#include <math.h>

#ifdef _WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <pthread.h>
#include <sys/time.h>    
#include <unistd.h>
#endif

#include "SUPREME_GnssTime.h"

#ifndef PI
#define PI       (4*atan(1.0))
#endif

/* time offset (s) */
static double timeoffset_ = 0.0;

/* Add time -------------------------------------
 * add time to gtime_t struct
 * args   : gtime_t t        I   gtime_t struct
 *          double sec       I   time to add (s)
 * return : gtime_t struct (t+sec)
 * ---------------------------------------------- */
gtime_t timeadd(gtime_t t, double sec)
{
	double tt;

	t.sec += sec; tt = floor(t.sec); t.time += (int)tt; t.sec -= tt;
	return t;
}

/* Time difference ------------------------------
 * difference between gtime_t structs
 * args   : gtime_t t1,t2    I   gtime_t structs
 * return : time difference (t1-t2) (s)
 * ---------------------------------------------- */
double timediff(gtime_t t1, gtime_t t2)
{
	return difftime(t1.time, t2.time) + t1.sec - t2.sec;
}

/* convert calendar day/time to time -----------------------------------
 * convert calendar day/time to gtime_t struct
 * args   : double *ep       I   day/time {year,month,day,hour,min,sec}
 * return : gtime_t struct
 * notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
 * --------------------------------------------------------------------- */
gtime_t epoch2time(const double ep[6])
{
	const int doy[] = { 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 };
	gtime_t time = { 0 };
	int days, sec, year = (int)ep[0], mon = (int)ep[1], day = (int)ep[2];

	if (year<1970 || 2099<year || mon<1 || 12<mon) return time;

	/* leap year if year%4==0 in 1901-2099 */
	days = (year - 1970) * 365 + (year - 1969) / 4 + doy[mon - 1] + day - 2 + (year % 4 == 0 && mon >= 3 ? 1 : 0);
	sec = (int)floor(ep[5]);
	time.time = (time_t)days * 86400 + (int)ep[3] * 3600 + (int)ep[4] * 60 + sec;
	time.sec = ep[5] - sec;
	return time;
}

/* time to calendar day/time -------------------------------------------
 * convert gtime_t struct to calendar day/time
 * args   : gtime_t t        I   gtime_t struct
 *          double *ep       O   day/time {year,month,day,hour,min,sec}
 * return : none
 * notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
 * ---------------------------------------------------------------------- */
void time2epoch(gtime_t t, double ep[6])
{
	const int mday[] =
	{ /* # of days in a month */
		31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
		31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
	};
	int days, sec, mon, day;

	/* leap year if year%4==0 in 1901-2099 */
	days = (int)(t.time / 86400);
	sec = (int)(t.time - (time_t)days * 86400);
	for (day = days % 1461, mon = 0; mon<48; mon++)
	{
		if (day >= mday[mon]) day -= mday[mon]; else break;
	}
	ep[0] = 1970 + days / 1461 * 4 + mon / 12; ep[1] = mon % 12 + 1; ep[2] = day + 1;
	ep[3] = sec / 3600; ep[4] = sec % 3600 / 60; ep[5] = sec % 60 + t.sec;
}

/* gps time to time ------------------------------------------
 * convert week and tow in gps time to gtime_t struct
 * args   : int    week      I   week number in gps time
 *          double sec       I   time of week in gps time (tow) (s)
 * return : gtime_t struct
 * ----------------------------------------------------------- */
gtime_t gpst2time(int week, double sec)
{
	gtime_t t = epoch2time(gpst0);

	if (sec<-1E9 || 1E9<sec) sec = 0.0;
	t.time += 86400 * 7 * week + (int)sec;
	t.sec = sec - (int)sec;
	return t;
}

/* time to gps time -------------------------------------------------------
 * convert gtime_t struct to week and tow in gps time
 * args   : gtime_t t        I   gtime_t struct
 *          int    *week     IO  week number in gps time (NULL: no output)
 * return : time of week in gps time (tow) (s)
 * ------------------------------------------------------------------------ */
double time2gpst(gtime_t t, int *week)
{
	gtime_t t0 = epoch2time(gpst0);
	time_t sec = t.time - t0.time;
	int w = (int)(sec / (86400 * 7));

	if (week) *week = w;
	return (double)(sec - w * 86400 * 7) + t.sec;
}

/* galileo system time to time ----------------------------------------
 * convert week and tow in galileo system time (gst) to gtime_t struct
 * args   : int    week      I   week number in gst
 *          double sec       I   time of week in gst (s)
 * return : gtime_t struct
 * -------------------------------------------------------------------- */
gtime_t gst2time(int week, double sec)
{
	gtime_t t = epoch2time(gst0);

	if (sec<-1E9 || 1E9<sec) sec = 0.0;
	t.time += 86400 * 7 * week + (int)sec;
	t.sec = sec - (int)sec;
	return t;
}

/* time to galileo system time ----------------------------------------
 * convert gtime_t struct to week and tow in galileo system time (gst)
 * args   : gtime_t t        I   gtime_t struct
 *          int    *week     IO  week number in gst (NULL: no output)
 * return : time of week in gst (s)
 * -------------------------------------------------------------------- */
double time2gst(gtime_t t, int *week)
{
	gtime_t t0 = epoch2time(gst0);
	time_t sec = t.time - t0.time;
	int w = (int)(sec / (86400 * 7));

	if (week) *week = w;
	return (double)(sec - w * 86400 * 7) + t.sec;
}

/* beidou time (bdt) to time ----------------------------------
 * convert week and tow in beidou time (bdt) to gtime_t struct
 * args   : int    week      I   week number in bdt
 *          double sec       I   time of week in bdt (s)
 * return : gtime_t struct
 * ------------------------------------------------------------ */
gtime_t bdt2time(int week, double sec)
{
	gtime_t t = epoch2time(bdt0);

	if (sec<-1E9 || 1E9<sec) sec = 0.0;
	t.time += 86400 * 7 * week + (int)sec;
	t.sec = sec - (int)sec;
	return t;
}

/* time to beidou time (bdt) -----------------------------------------
 * convert gtime_t struct to week and tow in beidou time (bdt)
 * args   : gtime_t t        I   gtime_t struct
 *          int    *week     IO  week number in bdt (NULL: no output)
 * return : time of week in bdt (s)
 * ------------------------------------------------------------------- */
double time2bdt(gtime_t t, int *week)
{
	gtime_t t0 = epoch2time(bdt0);
	time_t sec = t.time - t0.time;
	int w = (int)(sec / (86400 * 7));

	if (week) *week = w;
	return (double)(sec - w * 86400 * 7) + t.sec;
}

/* gpstime to utc -----------------------------------------
 * convert gpstime to utc considering leap seconds
 * args   : gtime_t t        I   time expressed in gpstime
 * return : time expressed in utc
 * notes  : ignore slight time offset under 100 ns
 * -------------------------------------------------------- */
gtime_t gpst2utc(gtime_t t)
{
	gtime_t tu;
	int i;

	for (i = 0; leaps[i][0]>0; i++) {
		tu = timeadd(t, leaps[i][6]);
		if (timediff(tu, epoch2time(leaps[i])) >= 0.0) return tu;
	}
	return t;
}

/* utc to gpstime -------------------------------------
 * convert utc to gpstime considering leap seconds
 * args   : gtime_t t        I   time expressed in utc
 * return : time expressed in gpstime
 * notes  : ignore slight time offset under 100 ns
 * ---------------------------------------------------- */
gtime_t utc2gpst(gtime_t t)
{
	int i;

	for (i = 0; leaps[i][0]>0; i++) {
		if (timediff(t, epoch2time(leaps[i])) >= 0.0) return timeadd(t, -leaps[i][6]);
	}
	return t;
}

/* gpst to bdt */
gtime_t gpst2bdt(gtime_t t)
{
	return timeadd(t, -14.0);
}

/* bdt to gpstime --------------------------------------------------
 * convert bdt (beidou navigation satellite system time) to gpstime
 * args   : gtime_t t        I   time expressed in bdt
 * return : time expressed in gpstime
 * notes  : see gpst2bdt()
 * ----------------------------------------------------------------- */
gtime_t bdt2gpst(gtime_t t)
{
	return timeadd(t, 14.0);
}

/* time to day and sec ------------------------------------
 * convert gtime_t struct to day and second of day
 * args   : gtime_t  time     I      gtime_t struct
 *          gtime_t  *day     O      gtime_t struct of day
 * return : second of day (s)
 * -------------------------------------------------------- */
static double time2sec(gtime_t time, gtime_t *day)
{
	double ep[6], sec;
	time2epoch(time, ep);
	sec = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];
	ep[3] = ep[4] = ep[5] = 0.0;
	*day = epoch2time(ep);
	return sec;
}

/* utc to gmst ----------------------------------------
 * convert utc to gmst (Greenwich mean sidereal time)
 * args   : gtime_t t        I   time expressed in utc
 *          double ut1_utc   I   UT1-UTC (s)
 * return : gmst (rad)
 * ---------------------------------------------------- */
double utc2gmst(gtime_t t, double ut1_utc)
{
	const double ep2000[] = { 2000, 1, 1, 12, 0, 0 };
	gtime_t tut, tut0;
	double ut, t1, t2, t3, gmst0, gmst;

	tut = timeadd(t, ut1_utc);
	ut = time2sec(tut, &tut0);
	t1 = timediff(tut0, epoch2time(ep2000)) / 86400.0 / 36525.0;
	t2 = t1*t1; t3 = t2*t1;
	gmst0 = 24110.54841 + 8640184.812866*t1 + 0.093104*t2 - 6.2E-6*t3;
	gmst = gmst0 + 1.002737909350795*ut;

	return fmod(gmst, 86400.0)*PI / 43200.0; /* 0 <= gmst <= 2*PI */
}

double mjd2jd(double mjd)
{
	return (mjd + 2400000.5);
}

double jd2mjd(double jd)
{
	return (jd - 2400000.5);
}

/* time to day of year -------------------------
 * convert time to day of year
 * args   : gtime_t t        I   gtime_t struct
 * return : day of year (days)
 * --------------------------------------------- */
double time2doy(gtime_t t)
{
	double ep[6];

	time2epoch(t, ep);
	ep[1] = ep[2] = 1.0; ep[3] = ep[4] = ep[5] = 0.0;
	return timediff(t, epoch2time(ep)) / 86400.0 + 1.0;
}

/* get day of week by second of week ----
 * return : day of week
 * -------------------------------------- */
int tow2dow(double secofweek)
{
	int dow = 0;
	dow = (int)(secofweek / 86400);
	return dow;
}

/* get current time in utc ------
 * get current time in utc
 * args   : none
 * return : current time in utc
 * ------------------------------ */
gtime_t timeget(void)
{
	gtime_t time;
	double ep[6] = { 0 };
#ifdef _WIN32
	SYSTEMTIME ts;

	GetSystemTime(&ts); /* utc */
	ep[0] = ts.wYear; ep[1] = ts.wMonth;  ep[2] = ts.wDay;
	ep[3] = ts.wHour; ep[4] = ts.wMinute; ep[5] = ts.wSecond + ts.wMilliseconds*1E-3;
#else
	struct timeval tv;
	struct tm *tt;

	if (!gettimeofday(&tv, NULL) && (tt = gmtime(&tv.tv_sec))) {
		ep[0] = tt->tm_year + 1900; ep[1] = tt->tm_mon + 1; ep[2] = tt->tm_mday;
		ep[3] = tt->tm_hour; ep[4] = tt->tm_min; ep[5] = tt->tm_sec + tv.tv_usec*1E-6;
	}
#endif
	time = epoch2time(ep);

#ifdef CPUTIME_IN_GPST /* cputime operated in gpst */
	time = gpst2utc(time);
#endif
	return timeadd(time, timeoffset_);
}

/* set current time in utc -----------------------------------------
 * set current time in utc
 * args   : gtime_t          I   current time in utc
 * return : none
 * notes  : just set time offset between cpu time and current time
 *          the time offset is reflected to only timeget()
 *          not reentrant
 * ------------------------------------------------------------------ */
void timeset(gtime_t t)
{
	timeoffset_ += timediff(t, timeget());
}

int LongYear(int yy, int mm)
{
	if (!(yy % 4) && (!(yy % 400) || (yy % 100)))
	{
		if (!mm || mm == 2)	return 1;
	}
	return 0;
}


gnsstime::gnsstime()
{
	this->Init();
}

gnsstime::gnsstime(double jd)
{
	m_Jd = jd;
	_jd2mjd();

	_mjd2date();
	_date2gpst();

	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];
}

gnsstime::gnsstime(double ep[6])
{
	m_Epoch[0] = ep[0]; m_Epoch[1] = ep[1]; m_Epoch[2] = ep[2];
	m_Epoch[3] = ep[3]; m_Epoch[4] = ep[4]; m_Epoch[5] = ep[5];

	m_Year = (int)ep[0]; m_Month = (int)ep[1]; m_Day = (int)ep[2];
	m_Hour = (int)ep[3]; m_Min = (int)ep[4]; m_Sec = ep[5];

	_date2gpst();
}

gnsstime::gnsstime(int ww, double ss)
{
	m_Week = ww; m_SecondofWeek = ss;

	_gpst2date();

	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	_date2jd();
}

gnsstime::gnsstime(int YY, int MM, int DD, int hh, int mm, double ss)
{
	m_Epoch[0] = YY; m_Epoch[1] = MM; m_Epoch[2] = DD;
	m_Epoch[3] = hh; m_Epoch[4] = mm; m_Epoch[5] = ss;
	m_Year = YY; m_Month = MM; m_Day = DD;
	m_Hour = hh; m_Min = mm; m_Sec = ss;

	m_gtime = epoch2time(m_Epoch);
	_date2gpst(); _date2jd();
	_GetDoy();
}

void gnsstime::Init()
{
	m_Epoch[0] = 0.0; m_Epoch[1] = 0.0; m_Epoch[2] = 0.0;
	m_Epoch[3] = 0.0; m_Epoch[4] = 0.0; m_Epoch[5] = 0.0;

	m_Year = 0; m_Month = 0; m_Day = 0;
	m_Hour = 0; m_Min = 0; m_Sec = 0.0;

	m_Doy = 0; m_Week = 0; m_Dow = 0;
	m_Jd = 0.0; m_Mjd = 0.0;
	m_SecondofWeek = 0.0;

	m_gtime = epoch2time(m_Epoch);
}

double gnsstime::_jd2mjd()
{
	m_Mjd = jd2mjd(m_Jd);

	return m_Mjd;
}

double gnsstime::_mjd2jd()
{
	m_Jd = mjd2jd(m_Mjd);

	return m_Jd;
}

double gnsstime::_date2jd()
{
	int m_year = m_Year, m_month = m_Month, m_day = m_Day;
	double m_hour = m_Hour + m_Min / 60.0 + m_Sec / 3600.0;
	if (m_month <= 2)
	{
		m_year -= 1;
		m_month += 12;
	}
	m_Jd = (int)(365.25*m_year) + (int)(30.6001*(m_month + 1)) + m_day + m_hour / 24.0 + 1720981.5;
	this->_jd2mjd();

	return m_Jd;
}

double gnsstime::_date2mjd()
{
	this->_date2jd();

	return m_Mjd;
}

int gnsstime::_GetDoy()
{
	m_Doy = (int)time2doy(m_gtime);

	return m_Doy;
}

int gnsstime::_date2utc()
{
	m_Epoch[0] = m_Year; m_Epoch[1] = m_Month; m_Epoch[2] = m_Day;
	m_Epoch[3] = m_Hour; m_Epoch[4] = m_Min; m_Epoch[5] = m_Sec;

	m_gtime = epoch2time(m_Epoch);

	return 1;
}

int gnsstime::_date2gpst()
{
	m_Epoch[0] = m_Year; m_Epoch[1] = m_Month; m_Epoch[2] = m_Day;
	m_Epoch[3] = m_Hour; m_Epoch[4] = m_Min; m_Epoch[5] = m_Sec;

	m_gtime = epoch2time(m_Epoch);
	m_SecondofWeek = time2gpst(m_gtime, &m_Week);

	m_Doy = _GetDoy();
	m_Dow = tow2dow(m_SecondofWeek);

	return 1;
}

int gnsstime::_gpst2date()
{
	m_gtime = gpst2time(m_Week, m_SecondofWeek);
	time2epoch(m_gtime, m_Epoch);
	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	return 1;
}

int gnsstime::_date2bdst()
{
	m_Epoch[0] = m_Year; m_Epoch[1] = m_Month; m_Epoch[2] = m_Day;
	m_Epoch[3] = m_Hour; m_Epoch[4] = m_Min; m_Epoch[5] = m_Sec;
	m_gtime = epoch2time(m_Epoch);
	m_SecondofWeek = time2bdt(m_gtime, &m_Week);

	return 1;
}

int gnsstime::_bdst2date()
{
	m_gtime = bdt2time(m_Week, m_SecondofWeek);
	time2epoch(m_gtime, m_Epoch);
	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	return 1;
}

int gnsstime::_date2ep()
{
	m_Epoch[0] = m_Year; m_Epoch[1] = m_Month; m_Epoch[2] = m_Day;
	m_Epoch[3] = m_Hour;	m_Epoch[4] = m_Min;	m_Epoch[5] = m_Sec;

	return 1;
}

int gnsstime::_ep2date()
{
	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	return 1;
}

int gnsstime::_jd2date()
{
	int a = (int)(m_Jd + 0.5);
	int b = a + 1537;
	int c = (int)((b - 122.1) / 365.25);
	int d = (int)(365.25*c);
	int e = (int)((b - d) / 30.6001);

	double day_t = m_Jd + 0.5 + b - (int)(30.6001*e) - d - a;
	m_Day = (int)(day_t);
	m_Month = e - 1 - 12 * (int)(e / 14.0);
	m_Year = c - 4715 - (int)((m_Month + 7) / 10.0);

	double h = (day_t - m_Day) * 24;

	if (fabs(h - (int)(h + 0.5))<1.0e-3)
	{
		h = (double)((int)(h + 0.5));
	}

	m_Hour = (int)(h);

	double min_t = (h - (int)(h)) * 60;

	if (fabs(min_t - (int)(min_t + 0.5))<1.0e-3)
	{
		min_t = (double)((int)(min_t + 0.5));
	}

	m_Min = (int)(min_t);
	m_Sec = (min_t - (int)(min_t))*60.0;

	this->_date2ep();

	return 1;
}

int gnsstime::_mjd2date()
{
	_mjd2jd();
	_jd2date();

	return 1;
}

/// Convert UTC to GPS Time
int gnsstime::Utc2Gpst()
{
	m_gtime = utc2gpst(m_gtime);
	m_SecondofWeek = time2gpst(m_gtime, &m_Week);

	time2epoch(m_gtime, m_Epoch);
	m_Year = (int)m_Epoch[0]; int m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	_GetDoy(); _date2jd();

	m_Dow = tow2dow(m_SecondofWeek);

	return 1;
}

/// Convert GPS to UTC Time
int gnsstime::Gpst2Utc()
{
	m_gtime = gpst2utc(m_gtime);

	time2epoch(m_gtime, m_Epoch);
	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	_GetDoy(); _date2jd();

	m_Type = UTC;

	return 1;
}

/// Convert UTC to BDS Time
int gnsstime::Utc2Bdst()
{
	m_gtime = utc2gpst(m_gtime);
	m_gtime = gpst2bdt(m_gtime);
	m_SecondofWeek = time2bdt(m_gtime, &m_Week);

	time2epoch(m_gtime, m_Epoch);
	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	_GetDoy(); _date2jd();

	return 1;
}

/// Convert BDS to UTC Time
int gnsstime::Bdst2Utc()
{
	m_gtime = bdt2gpst(m_gtime);
	m_gtime = gpst2utc(m_gtime);

	time2epoch(m_gtime, m_Epoch);
	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	_GetDoy(); _date2jd();

	return 1;
}

/// Convert GPS to BDS Time
int gnsstime::Gpst2Bdst()
{
	m_gtime = gpst2bdt(m_gtime);
	m_SecondofWeek = time2bdt(m_gtime, &m_Week);

	time2epoch(m_gtime, m_Epoch);
	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	_GetDoy(); _date2jd();

	return 1;
}

/// BDS to GPS Time
int gnsstime::Bdst2Gpst()
{
	m_gtime = bdt2gpst(m_gtime);
	m_SecondofWeek = time2gpst(m_gtime, &m_Week);

	time2epoch(m_gtime, m_Epoch);
	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	_GetDoy(); _date2jd();

	return 1;
}

/// Get local time (UTC)
void gnsstime::_GetLocalTime()
{
	m_gtime = timeget();

	time2epoch(m_gtime, m_Epoch);
	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];

	_GetDoy(); _date2jd();
}

/// Get local gps time (GPST)
void gnsstime::_GetLocalTime_GPST()
{
	_GetLocalTime();

	Utc2Gpst();
}

void gnsstime::_GetLocalTime_GPST(double offset)
{
	_GetLocalTime_GPST();

	m_SecondofWeek -= offset;
	if (m_SecondofWeek < 0)
	{
		m_SecondofWeek += 604800;
		m_Week--;
	}
	else if (m_SecondofWeek >= 604800)
	{
		m_SecondofWeek -= 604800;
		m_Week++;
	}

	m_gtime = gpst2time(m_Week, m_SecondofWeek);

	time2epoch(m_gtime, m_Epoch);
	m_Year = (int)m_Epoch[0]; m_Month = (int)m_Epoch[1]; m_Day = (int)m_Epoch[2];
	m_Hour = (int)m_Epoch[3]; m_Min = (int)m_Epoch[4]; m_Sec = m_Epoch[5];
}

void gnsstime::_GetDaySecond()
{
	_day = (int)(367 * (m_Year - 1950) - (7 * (m_Year + (m_Month + 9) / 12)) / 4 + (275 * m_Month) / 9 + m_Day + 3381);
	_second = 3600.0*(double)m_Hour + 60.0*(double)m_Min + m_Sec;
}

int gnsstime::_LongYear()
{
	if (LongYear(m_Year, m_Month))
		return 1;

	return 0;
}

char* gnsstime::_Time2Str()
{
	memset(_TimeStr, 0, sizeof(_TimeStr));
	sprintf_s(_TimeStr, "%04d/%02d/%02d %02d:%02d:%02d", 
		m_Year, m_Month, m_Day, m_Hour, m_Min, (int)m_Sec);
	return _TimeStr;
}

double operator-(const gnsstime &T1, const gnsstime &T2)
{
	return (T1.m_Week - T2.m_Week)*(double)(604800.0) + T1.m_SecondofWeek - T2.m_SecondofWeek;
}

gnsstime operator+(const gnsstime &T, const double &sec)
{
	gnsstime re_time;
	re_time.m_gtime.time = T.m_gtime.time + (int)sec;
	re_time.m_gtime.sec = T.m_gtime.sec + sec - (int)sec;

	time2epoch(re_time.m_gtime, re_time.m_Epoch);

	re_time.m_Year  = (int)re_time.m_Epoch[0];
	re_time.m_Month = (int)re_time.m_Epoch[1];
	re_time.m_Day   = (int)re_time.m_Epoch[2];
	re_time.m_Hour  = (int)re_time.m_Epoch[3];
	re_time.m_Min   = (int)re_time.m_Epoch[4];
	re_time.m_Sec   = re_time.m_Epoch[5];

	re_time._date2gpst();
	time2epoch(re_time.m_gtime, re_time.m_Epoch);

	re_time.m_Dow = tow2dow(re_time.m_SecondofWeek);
	re_time._GetDoy(); re_time._date2jd();

	return re_time;
}

bool operator>(const gnsstime &T1, const gnsstime &T2)
{
	if (T1 - T2 > 0) return true;
	else return false;
}
bool operator<(const gnsstime &T1, const gnsstime &T2)
{
	if (T1 - T2 < 0) return true;
	else return false;
}
bool operator>=(const gnsstime &T1, const gnsstime &T2)
{
	if (T1 - T2 >= 0) return true;
	else return false;
}
bool operator<=(const gnsstime &T1, const gnsstime &T2)
{
	if (T1 - T2 <= 0) return true;
	else return false;
}
bool operator==(const gnsstime &T1, const gnsstime &T2)
{
	if (T1 - T2 == 0) return true;
	else return false;
}
