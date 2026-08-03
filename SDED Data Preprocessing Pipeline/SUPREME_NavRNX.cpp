#include "SUPREME_NavRNX.h"
#include "SUPREME_CommonFunction.h"

#define  NAV_FILE_LINE                   100
#define  PRECISE_LINE_CHARACTER_NUM      100

void Nav_ana_Time(char* strline, unsigned int len, int flag, gnsstime& t)
{
	char s_year[5] = "", s_month[5] = "", s_day[5] = "", s_hour[5] = "", s_minu[5] = "", s_sec[5] = "";

	if (flag == 2)
	{
		GetSubStr(strline, len, s_year, 5, 3, 2);
		GetSubStr(strline, len, s_month, 5, 6, 2);
		GetSubStr(strline, len, s_day, 5, 9, 2);
		GetSubStr(strline, len, s_hour, 5, 12, 2);
		GetSubStr(strline, len, s_minu, 5, 15, 2);
		GetSubStr(strline, len, s_sec, 5, 18, 4);
	}
	else
	{
		GetSubStr(strline, len, s_year, 5, 4, 4);
		GetSubStr(strline, len, s_month, 5, 9, 2);
		GetSubStr(strline, len, s_day, 5, 12, 2);
		GetSubStr(strline, len, s_hour, 5, 15, 2);
		GetSubStr(strline, len, s_minu, 5, 18, 2);
		GetSubStr(strline, len, s_sec, 5, 21, 2);
	}

	t.m_Year = str_to_i(s_year);
	if (flag == 2)
	{
		if (t.m_Year>80){ t.m_Year += 1900; }
		else
		{
			t.m_Year += 2000;
		}
	}
	t.m_Month = str_to_i(s_month);
	t.m_Day = str_to_i(s_day);
	t.m_Min = str_to_i(s_minu);
	t.m_Hour = str_to_i(s_hour);
	t.m_Sec = str_to_f(s_sec);
}

/// Get epoch time in precise ephemeris string line
void Pre_ana_Epochtime(gnsstime& time, char *strline, unsigned int len)
{
	char year[5] = "", month[3] = "", day[3] = "";
	char hour[3] = "", minute[3] = "", sec[12] = "";

	GetSubStr(strline, len, year, 5, 3, 4);
	time.m_Year = str_to_i(year);
	GetSubStr(strline, len, month, 3, 8, 2);
	time.m_Month = str_to_i(month);
	GetSubStr(strline, len, day, 3, 11, 2);
	time.m_Day = str_to_i(day);
	GetSubStr(strline, len, hour, 3, 14, 2);
	time.m_Hour = str_to_i(hour);
	GetSubStr(strline, len, minute, 3, 17, 2);
	time.m_Min = str_to_i(minute);
	GetSubStr(strline, len, sec, 12, 20, 11);
	time.m_Sec = str_to_f(sec);

	time._date2jd();
	time._date2gpst();
}

/// Get line in precise ephemeris
unsigned int Pre_GetLine(FILE* fp, char *strline, unsigned int len)
{
	if (is_eof(fp)) // end of file
	{
		return 0;
	}

	GetFileLine(fp, strline, len);

	switch (strline[0])
	{
	case '*':return 1;       // time string line
	case 'P':return 2;       // position string line
	case 'V':return 3;       // velocity string line
	default:return 4;        // ohter line
	}
}

/// Find Gps Ephemeris by prn
Nav_GpsEph* Nav_GPS_FindPrn(unsigned int prn, NavData* navdata)
{
	if (navdata == NULL) return NULL;

	Nav_GpsEph* ptr = navdata->m_Nav_GpsEph;
	while (ptr != NULL)
	{
		if (ptr->m_Prn == prn)
			break;
		ptr = ptr->m_Next;
	}
	return ptr;
}

/// Find Glo Ephemeris by prn
Nav_GloEph* Nav_GLO_FindPrn(unsigned int prn, NavData* navdata)
{
	Nav_GloEph* ptr = navdata->m_Nav_GloEph;
	while (ptr != NULL)
	{
		if (ptr->m_Prn == prn)
			break;
		ptr = ptr->m_Next;
	}
	return ptr;
}

/// Find Gal Ephemeris by prn
Nav_GalEph* Nav_GAL_FindPrn(unsigned int prn, NavData* navdata)
{
	Nav_GalEph* ptr = navdata->m_Nav_GalEph;
	while (ptr != NULL)
	{
		if (ptr->m_Prn == prn)
			break;
		ptr = ptr->m_Next;
	}
	return ptr;
}

/// Find Bds Ephemeris by prn
Nav_BdsEph* Nav_BDS_FindPrn(unsigned int prn, NavData* navdata)
{
	Nav_BdsEph* ptr = navdata->m_Nav_BdsEph;
	while (ptr != NULL)
	{
		if (ptr->m_Prn == prn)
			break;
		ptr = ptr->m_Next;
	}
	return ptr;
}

void Pre_EpochData_ptr_Advance(PreEpochData **ptr, int count)
{
	int n = abs(count);
	int i = 0;

	if (count > 0)
	{
		for (i = 0; i < n; ++i)
		{
			*ptr = (*ptr)->m_Next;
		}
	}
	else
	{
		for (i = 0; i < n; ++i)
		{
			*ptr = (*ptr)->m_Previous;
		}
	}
}

int Precise_Interp_Pos(gnsstime& time, PreEpochData* ptr, int len, double pos[3], int n)
{
	if (ptr == NULL || len == 0 || len != (n + 1))
	{
		pos[0] = 0;
		pos[1] = 0;
		pos[2] = 0;
		return 0;
	}

	if (time.m_Jd<ptr[0].m_GnssTime.m_Jd || time.m_Jd>ptr[len - 1].m_GnssTime.m_Jd)
	{
		pos[0] = 0;
		pos[1] = 0;
		pos[2] = 0;
		return 0;
	}

	int i = 0, j = 0;
	for (i = 0; i<3; ++i)
	{
		pos[i] = 0;
	}

	for (i = 0; i != len; ++i)
	{
		double temp = 1.0;

		for (j = 0; j != len; ++j)
		{
			if (i == j)
			{
				continue;
			}
			
			temp *= (time - (ptr[j].m_GnssTime)) / ((ptr[i].m_GnssTime) - (ptr[j].m_GnssTime));
		}
		pos[0] += (temp*ptr[i].m_X);
		pos[1] += (temp*ptr[i].m_Y);
		pos[2] += (temp*ptr[i].m_Z);
	}

	return 1;
}

int Precise_Interp_Vel(gnsstime& time, PreEpochData* ptr, int len, double vel[3], int n)
{
	if (ptr == NULL || len == 0 || len != (n + 1))
	{
		vel[0] = 0;
		vel[1] = 0;
		vel[2] = 0;
		return 0;
	}

	if (time.m_Jd<ptr[0].m_GnssTime.m_Jd || time.m_Jd>ptr[len - 1].m_GnssTime.m_Jd)
	{
		vel[0] = 0;
		vel[1] = 0;
		vel[2] = 0;
		return 0;
	}

	int i = 0, j = 0, k = 0;

	for (i = 0; i<3; ++i)
	{
		vel[i] = 0.0;
	}

	// interpolate position.
	for (i = 0; i != len; ++i)
	{
		double temp_t = 0.0;

		for (j = 0; j != len; ++j)
		{
			if (i == j)
			{
				continue;
			}
			double t_t = ((ptr[i].m_GnssTime) - (ptr[j].m_GnssTime));

			double p_t = 1.0;

			for (k = 0; k != len; ++k)
			{
				if (k == i || k == j)
				{
					continue;
				}
				p_t *= (time - (ptr[k].m_GnssTime)) / ((ptr[i].m_GnssTime) - (ptr[k].m_GnssTime));
			}

			temp_t += (p_t / t_t);
		}
		vel[0] += (temp_t*ptr[i].m_X);
		vel[1] += (temp_t*ptr[i].m_Y);
		vel[2] += (temp_t*ptr[i].m_Z);
	}

	return 1;
}

int GpsEph::Nav_ana_GpsEph(FILE* fp, char* strline, unsigned int len, int flag)
{
	Nav_ana_Time(strline, len, flag, m_Toc_Gps);
	m_Toc_Gps._date2jd();
	m_Toc_Gps._date2gpst();

	char str_tmp[20] = "";
	unsigned int i = 4, j = 23, k = 0;

	if (flag == 2)
	{
		i = 3;
		j = 22;
	}

	// first line.
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); ++k;
	m_clock_bias = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); ++k;
	m_clock_drift = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); k = 0;
	m_clock_driftrate = atof(str_tmp);

	// 1-orbits
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_IODE = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Crs = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Delta_n = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_M0 = atof(str_tmp);

	// 2-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Cuc = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_e = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Cus = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_sqrt_A = atof(str_tmp);

	// 3-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_TOE = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Cic = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_OMEGA0 = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_Cis = atof(str_tmp);

	// 4-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_i0 = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Crc = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_omega = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_OMEGADOT = atof(str_tmp);

	// 5-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_IDOT = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_L2code = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_GPSweek = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_L2P = atof(str_tmp);

	// 6-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_precision = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_SVhealth = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_TGD = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_IODC = atof(str_tmp);

	// 7-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_TOW = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_h = atof(str_tmp);

	return 1;
}

/// Display gps ephemeris
void GpsEph::Display_GpsEph()
{
	printf("%03d %04d %02d %02d %02d %02d %5.2f",
		m_prn, m_Toc_Gps.m_Year, m_Toc_Gps.m_Month, m_Toc_Gps.m_Day,
		m_Toc_Gps.m_Hour, m_Toc_Gps.m_Min, m_Toc_Gps.m_Sec);
	printf("%22.12e%22.12e%22.12e\n", m_clock_bias, m_clock_drift, m_clock_driftrate);

	printf("    %22.12e%22.12e%22.12e%22.12e\n", m_IODE, m_Crs, m_Delta_n, m_M0);
	printf("    %22.12e%22.12e%22.12e%22.12e\n", m_Cuc, m_e, m_Cus, m_sqrt_A);
	printf("    %22.12e%22.12e%22.12e%22.12e\n", m_TOE, m_Cic, m_OMEGA0, m_Cis);
	printf("    %22.12e%22.12e%22.12e%22.12e\n", m_i0, m_Crc, m_omega, m_OMEGADOT);
	printf("    %22.12e%22.12e%22.12e%22.12e\n", m_IDOT, m_L2code, m_GPSweek, m_L2P);
	printf("    %22.12e%22.12e%22.12e%22.12e\n", m_precision, m_SVhealth, m_TGD, m_IODC);
	printf("    %22.12e%22.12e\n", m_TOW, m_h);
}

/// Write gps ephemeris in file
void GpsEph::Write_GpsEph(FILE *fp)
{
	if (fp)
	{
		fprintf(fp, "%c%02d %04d %02d %02d %02d %02d %5.2f",
			GetSystem(m_prn), m_prn, m_Toc_Gps.m_Year, m_Toc_Gps.m_Month, m_Toc_Gps.m_Day,
			m_Toc_Gps.m_Hour, m_Toc_Gps.m_Min, m_Toc_Gps.m_Sec);
		fprintf(fp, "%22.12e%22.12e%22.12e\n", m_clock_bias, m_clock_drift, m_clock_driftrate);

		fprintf(fp, "    %22.12e%22.12e%22.12e%22.12e\n", m_IODE, m_Crs, m_Delta_n, m_M0);
		fprintf(fp, "    %22.12e%22.12e%22.12e%22.12e\n", m_Cuc, m_e, m_Cus, m_sqrt_A);
		fprintf(fp, "    %22.12e%22.12e%22.12e%22.12e\n", m_TOE, m_Cic, m_OMEGA0, m_Cis);
		fprintf(fp, "    %22.12e%22.12e%22.12e%22.12e\n", m_i0, m_Crc, m_omega, m_OMEGADOT);
		fprintf(fp, "    %22.12e%22.12e%22.12e%22.12e\n", m_IDOT, m_L2code, m_GPSweek, m_L2P);
		fprintf(fp, "    %22.12e%22.12e%22.12e%22.12e\n", m_precision, m_SVhealth, m_TGD, m_IODC);
		fprintf(fp, "    %22.12e%22.12e\n", m_TOW, m_h);
	}
}

int GloEph::Nav_ana_GloEph(FILE* fp, char* strline, unsigned int len, int flag)
{
	Nav_ana_Time(strline, len, flag, m_Toc_Glo);
	m_Toc_Glo._date2jd();
	m_Toc_Glo._date2gpst();

	gnsstime tmp = m_Toc_Glo;
	tmp.Utc2Gpst();
	m_Toc_Gps = tmp;

	char str_tmp[20] = "";
	unsigned int i = 4, j = 23, k = 0;

	if (flag == 2)
	{
		i = 3;
		j = 22;
	}

	// first line.
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); ++k;
	m_tau = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); ++k;
	m_gamma = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); k = 0;
	m_tk = atof(str_tmp);

	// 1-orbits
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_X = atof(str_tmp)*1000.0;
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Vx = atof(str_tmp)*1000.0;
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Ax = atof(str_tmp)*1000.0;
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_health = atof(str_tmp);

	// 2-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Y = atof(str_tmp)*1000.0;
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Vy = atof(str_tmp)*1000.0;
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Ay = atof(str_tmp)*1000.0;
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_frequency_number = atof(str_tmp);

	// 3-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Z = atof(str_tmp)*1000.0;
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Vz = atof(str_tmp)*1000.0;
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Az = atof(str_tmp)*1000.0;
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_E = atof(str_tmp);

	return 1;
}

/// Display glonass ephemeris
void GloEph::Display_GloEph()
{
	printf("%04d %02d %02d %02d %02d %5.2f",
		m_Toc_Glo.m_Year, m_Toc_Glo.m_Month, m_Toc_Glo.m_Day,
		m_Toc_Glo.m_Hour, m_Toc_Glo.m_Min, m_Toc_Glo.m_Sec);
	printf("%22.12e%22.12e%22.12e\n", m_tau, m_gamma, m_tk);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_X, m_Vx, m_Ax, m_health);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_Y, m_Vy, m_Ay, m_frequency_number);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_Z, m_Vz, m_Az, m_E);

}

int BdsEph::Nav_ana_BdsEph(FILE* fp, char* strline, unsigned int len, int flag)
{
	Nav_ana_Time(strline, len, flag, m_Toc_Bds);
	m_Toc_Bds._date2jd();
	m_Toc_Bds._date2bdst();

	gnsstime timetmp = m_Toc_Bds;
	timetmp.Bdst2Gpst();
	m_Toc_Gps = timetmp;

	char str_tmp[20] = "";
	unsigned int i = 4, j = 23, k = 0;

	if (flag == 2)
	{
		i = 3;
		j = 22;
	}

	// first line.
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); ++k;
	m_clock_bias = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); ++k;
	m_clock_drift = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); k = 0;
	m_clock_driftrate = atof(str_tmp);

	// 1-orbits
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_IODE = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Crs = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Delta_n = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_M0 = atof(str_tmp);

	// 2-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Cuc = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_e = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Cus = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_sqrt_A = atof(str_tmp);

	// 3-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_TOE = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Cic = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_OMEGA0 = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_Cis = atof(str_tmp);

	// 4-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_i0 = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Crc = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_omega = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_OMEGADOT = atof(str_tmp);

	// 5-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_IDOT = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_reserved1 = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_bdweek = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_reserved2 = atof(str_tmp);

	// 6-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_precision = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_svhealth = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_TGD1 = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_TGD2 = atof(str_tmp);

	// 7-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_IODC = atof(str_tmp);

	return 1;
}

/// Display bds ephemeris
void BdsEph::Display_BdsEph()
{
	printf("%04d %02d %02d %02d %02d %5.2f",
		m_Toc_Bds.m_Year, m_Toc_Bds.m_Month, m_Toc_Bds.m_Day,
		m_Toc_Bds.m_Hour, m_Toc_Bds.m_Min, m_Toc_Bds.m_Sec);
	printf("%22.12e%22.12e%22.12e\n", m_clock_bias, m_clock_drift, m_clock_driftrate);

	printf("%22.12e%22.12e%22.12e%22.12e\n", m_IODE, m_Crs, m_Delta_n, m_M0);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_Cuc, m_e, m_Cus, m_sqrt_A);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_TOE, m_Cic, m_OMEGA0, m_Cis);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_i0, m_Crc, m_omega, m_OMEGADOT);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_IDOT, m_reserved1, m_bdweek, m_reserved2);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_precision, m_svhealth, m_TGD1, m_TGD2);
	printf("%22.12e%22.12e\n", m_IODC, 0.0);

}

int GalEph::Nav_ana_GalEph(FILE* fp, char* strline, unsigned int len, int flag)
{
	Nav_ana_Time(strline, len, flag, m_Toc_Gal);
	m_Toc_Gal._date2jd();
	m_Toc_Gal._date2gpst();

	int week = m_Toc_Gal.m_Week;
	double tow = m_Toc_Gal.m_SecondofWeek;

	m_Toc_Gps.m_Week = week;
	m_Toc_Gps.m_SecondofWeek = tow;
	m_Toc_Gps._gpst2date();
	m_Toc_Gps._date2jd();

	char str_tmp[20] = "";
	unsigned int i = 4, j = 23, k = 0;

	if (flag == 2)
	{
		i = 3;
		j = 22;
	}

	// first line.
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); ++k;
	m_clock_bias = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); ++k;
	m_clock_drift = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, j + k * 19, 19); k = 0;
	m_clock_driftrate = atof(str_tmp);

	// 1-orbits
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_IODnav = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Crs = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Delta_n = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_M0 = atof(str_tmp);

	// 2-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Cuc = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_e = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Cus = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_sqrt_A = atof(str_tmp);

	// 3-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_TOE = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Cic = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_OMEGA0 = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_Cis = atof(str_tmp);

	// 4-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_i0 = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_Crc = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_omega = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_OMEGADOT = atof(str_tmp);

	// 5-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_IDOT = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_data_source = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_GALweek = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_reserved1 = atof(str_tmp);

	// 6-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_SISA = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_SVhealth = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_BGD_E1_E5a = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); k = 0;
	m_BGD_E1_E5b = atof(str_tmp);

	// 7-orbit
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_TOW = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_reserved2 = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_reserved3 = atof(str_tmp);
	GetSubStr(strline, len, str_tmp, 20, i + k * 19, 19); ++k;
	m_reserved4 = atof(str_tmp);

	return 1;
}

/// Display gal ephemeris
void GalEph::Display_GalEph()
{
	printf("%04d %02d %02d %02d %02d %5.2f",
		m_Toc_Gal.m_Year, m_Toc_Gal.m_Month, m_Toc_Gal.m_Day,
		m_Toc_Gal.m_Hour, m_Toc_Gal.m_Min, m_Toc_Gal.m_Sec);
	printf("%22.12e%22.12e%22.12e\n", m_clock_bias, m_clock_drift, m_clock_driftrate);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_IODnav, m_Crs, m_Delta_n, m_M0);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_Cuc, m_e, m_Cus, m_sqrt_A);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_TOE, m_Cic, m_OMEGA0, m_Cis);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_i0, m_Crc, m_omega, m_OMEGADOT);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_IDOT, m_data_source, m_GALweek, m_reserved1);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_SISA, m_SVhealth, m_BGD_E1_E5a, m_BGD_E1_E5b);
	printf("%22.12e%22.12e%22.12e%22.12e\n", m_TOW, m_reserved2, m_reserved3, m_reserved4);
}

/// Navigation Header Constructor
NavHeader::NavHeader()
{
	m_Sys = 0;
	m_Type = 0;
	m_Gps_Ion_flag = 0;
	m_Gal_Ion_flag = 0;
	m_Version = 0.0;
	memset((void*)m_Gps_ionA, 0, sizeof(double)* 4);
	memset((void*)m_Gps_ionB, 0, sizeof(double)* 4);
	memset((void*)m_Gal_ionA, 0, sizeof(double)* 4);
}

/// Read navigation file header
int NavHeader::Read_NavigationHeader(char* filename)
{
	if (filename == NULL) return 0;

	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("Error: Navigation file open failed.\n");
		return 0;
	}

	unsigned int len = NAV_FILE_LINE + 1;
	char strline[NAV_FILE_LINE + 1] = "";
	char sel1[21] = "RINEX VERSION / TYPE";

	if (SetPointerPos(fp, sel1, strline, len))
	{
		char version[10] = "";
		GetSubStr(strline, len, version, 10, 0, 9);
		m_Version = str_to_f(version);
		m_Type = strline[20];
		m_Sys = strline[40];
		if (m_Sys == ' ')
		{
			m_Sys = (m_Type == 'N') ? 'G' : 'R';
		}
		if ((m_Type != 'N') && (m_Type != 'G'))
		{
			printf("Warning: This file is not navigation file.\n");
			fclose(fp);
			return 0;
		}
	}
	else
	{
		printf("This navigation file is invalid.\n");
		fclose(fp);
		return 0;
	}

	int i = 0;
	char sel2[21] = "IONOSPHERIC CORR    ";

	if (SetPointerPos(fp, sel2, strline, len))
	{
		char iontype[5] = "";
		GetSubStr(strline, len, iontype, 5, 0, 4);
		if (strcmp(iontype, "GAL ") == 0)
		{
			m_Gal_Ion_flag = 1;
			for (i = 0; i<4; ++i)
			{
				char ionecf[13] = "";
				GetSubStr(strline, len, ionecf, 13, 5 + i * 12, 12);
				m_Gal_ionA[i] = atof(ionecf);
			}
		}
		else
		{
			if (strcmp(iontype, "GPSA") == 0)
			{
				m_Gps_Ion_flag = 1;
				for (i = 0; i<4; ++i)
				{
					char ionecf[13] = "";
					GetSubStr(strline, len, ionecf, 13, 5 + i * 12, 12);
					m_Gps_ionA[i] = atof(ionecf);
				}
				GetFileLine(fp, strline, len);
				for (i = 0; i<4; ++i)
				{
					char ionecf[13] = "";
					GetSubStr(strline, len, ionecf, 13, 5 + i * 12, 12);
					m_Gps_ionB[i] = atof(ionecf);
				}
			}
			if (strcmp(iontype, "GPSB") == 0)
			{
				m_Gps_Ion_flag = 1;
				for (i = 0; i<4; ++i)
				{
					char ionecf[13] = "";
					GetSubStr(strline, len, ionecf, 13, 5 + i * 12, 12);
					m_Gps_ionB[i] = atof(ionecf);
				}
				GetFileLine(fp, strline, len);
				for (i = 0; i<4; ++i)
				{
					char ionecf[13] = "";
					GetSubStr(strline, len, ionecf, 13, 5 + i * 12, 12);
					m_Gps_ionA[i] = atof(ionecf);
				}
			}
		}
	}

	fclose(fp);
	return 1;
}

/// Display navigation header information
void NavHeader::Display_NavHeadr()
{
	printf("version: %10.4f %c  %c\n", m_Version, m_Type, m_Sys);
	printf("   GPSA: %15.4e%15.4e%15.4e%15.4e\n", m_Gps_ionA[0], m_Gps_ionA[1], m_Gps_ionA[2], m_Gps_ionA[3]);
	printf("   GPSB: %15.4e%15.4e%15.4e%15.4e\n", m_Gps_ionB[0], m_Gps_ionB[1], m_Gps_ionB[2], m_Gps_ionB[3]);
	printf("   GAL : %15.4e%15.4e%15.4e%15.4e\n", m_Gal_ionA[0], m_Gal_ionA[1], m_Gal_ionA[2], m_Gal_ionA[3]);
}

NavData::NavData()
{
	m_Nav_Buff_Size = 0;

	m_GPS_PRN_Count = 0;
	m_GLO_PRN_Count = 0;
	m_BDS_PRN_Count = 0;
	m_GAL_PRN_Count = 0;

	/*memset((void*)(&gpsorb), 0, sizeof(rtcm_orbitv));
	memset((void*)(&gpsclk), 0, sizeof(rtcm_clockv));
	memset((void*)(&gloorb), 0, sizeof(rtcm_orbitv));
	memset((void*)(&gloclk), 0, sizeof(rtcm_clockv));
	memset((void*)(&galorb), 0, sizeof(rtcm_orbitv));
	memset((void*)(&galclk), 0, sizeof(rtcm_clockv));
	memset((void*)(&bdsorb), 0, sizeof(rtcm_orbitv));
	memset((void*)(&bdsclk), 0, sizeof(rtcm_clockv));*/
}

NavData::~NavData()
{
	this->Clear();
}

/// Add gps ephemeris to navigation
int NavData::Nav_Add_GpsEph(GpsEph* eph)
{
	Nav_GpsEph* ptr = Nav_GPS_FindPrn(eph->m_prn, this);

	if (m_Nav_GpsEph == NULL || m_GPS_PRN_Count == 0)
	{
		gpstail = NULL;
	}

	if (ptr == NULL)
	{
		Nav_GpsEph* p = (Nav_GpsEph*)(malloc(sizeof(Nav_GpsEph)));
		memset((void*)(p), 0, sizeof(Nav_GpsEph));

		p->m_EpochCount++;
		p->m_GpsEph = eph;
		p->m_Prn = eph->m_prn;
		p->m_Tail = eph;

		m_GPS_PRN_Count++;
		if (m_Nav_GpsEph == NULL)
		{
			gpstail = p;
			m_Nav_GpsEph = p;
		}
		else
		{
			gpstail->m_Next = p;
			gpstail = p;
		}
	}
	else
	{
		GpsEph* pt = ptr->m_GpsEph;

		while (pt != NULL)
		{
			if (pt->m_Toc_Gps.m_Week == eph->m_Toc_Gps.m_Week &&
				(int)(pt->m_Toc_Gps.m_SecondofWeek) == (int)(eph->m_Toc_Gps.m_SecondofWeek) &&
				(int)(pt->m_IODE) == (int)(eph->m_IODE))
			{
				return 0;
			}
			pt = pt->m_Next;
		}

		ptr->m_EpochCount++;
		ptr->m_Tail->m_Next = eph;
		ptr->m_Tail = eph;

		if (m_Nav_Buff_Size > 0)
		{
			if (ptr->m_EpochCount > m_Nav_Buff_Size)
			{
				GpsEph* m_GpsEph_tmp = ptr->m_GpsEph->m_Next;

				if (ptr->m_GpsEph)
					free(ptr->m_GpsEph);

				ptr->m_GpsEph = m_GpsEph_tmp;
				ptr->m_EpochCount--;
			}
		}
	}

	return 1;
}

/// Add glonass ephmeris to navigation
int NavData::Nav_Add_GloEph(GloEph* eph)
{
	Nav_GloEph* ptr = Nav_GLO_FindPrn(eph->m_prn, this);

	if (m_Nav_GloEph == NULL || m_GLO_PRN_Count == 0)
	{
		glotail = NULL;
	}

	if (ptr == NULL)
	{
		Nav_GloEph* p = (Nav_GloEph*)(malloc(sizeof(Nav_GloEph)));
		memset((void*)(p), 0, sizeof(Nav_GloEph));

		p->m_EpochCount++;
		p->m_GloEph = eph;
		p->m_Prn = eph->m_prn;
		p->m_Tail = eph;

		m_GLO_PRN_Count++;
		if (m_Nav_GloEph == NULL)
		{
			glotail = p;
			m_Nav_GloEph = p;
		}
		else
		{
			glotail->m_Next = p;
			glotail = p;
		}
	}
	else
	{
		GloEph* pt = ptr->m_GloEph;

		while (pt != NULL)
		{
			if (pt->m_Toc_Glo.m_Week == eph->m_Toc_Glo.m_Week &&
				(int)(pt->m_Toc_Glo.m_SecondofWeek) == (int)(eph->m_Toc_Glo.m_SecondofWeek))
			{
				return 0;
			}
			pt = pt->m_Next;
		}

		ptr->m_EpochCount++;
		ptr->m_Tail->m_Next = eph;
		ptr->m_Tail = eph;

		if (m_Nav_Buff_Size > 0)
		{
			if (ptr->m_EpochCount > m_Nav_Buff_Size)
			{
				GloEph* m_GloEph_tmp = ptr->m_GloEph->m_Next;

				if (ptr->m_GloEph)
					free(ptr->m_GloEph);

				ptr->m_GloEph = m_GloEph_tmp;
				ptr->m_EpochCount--;
			}
		}
	}

	return 1;
}

/// Add galileo ephemeris to navigation
int NavData::Nav_Add_GalEph(GalEph* eph)
{
	Nav_GalEph* ptr = Nav_GAL_FindPrn(eph->m_prn, this);

	if (m_Nav_GalEph == NULL || m_GAL_PRN_Count == 0)
	{
		galtail = NULL;
	}

	if (ptr == NULL)
	{
		Nav_GalEph* p = (Nav_GalEph*)(malloc(sizeof(Nav_GalEph)));
		memset((void*)(p), 0, sizeof(Nav_GalEph));

		p->m_EpochCount++;
		p->m_GalEph = eph;
		p->m_Prn = eph->m_prn;
		p->m_Tail = eph;

		m_GAL_PRN_Count++;
		if (m_Nav_GalEph == NULL)
		{
			galtail = p;
			m_Nav_GalEph = p;
		}
		else
		{
			galtail->m_Next = p;
			galtail = p;
		}
	}
	else
	{
		GalEph* pt = ptr->m_GalEph;

		while (pt != NULL)
		{
			if (pt->m_Toc_Gps.m_Week == eph->m_Toc_Gps.m_Week &&
				(int)(pt->m_Toc_Gps.m_SecondofWeek) == (int)(eph->m_Toc_Gps.m_SecondofWeek) &&
				(int)(pt->m_IODnav) == (int)(eph->m_IODnav))
			{
				return 0;
			}
			pt = pt->m_Next;
		}

		ptr->m_EpochCount++;
		ptr->m_Tail->m_Next = eph;
		ptr->m_Tail = eph;

		if (m_Nav_Buff_Size > 0)
		{
			if (ptr->m_EpochCount > m_Nav_Buff_Size)
			{
				GalEph* m_GalEph_tmp = ptr->m_GalEph->m_Next;

				if (ptr->m_GalEph)
					free(ptr->m_GalEph);

				ptr->m_GalEph = m_GalEph_tmp;
				ptr->m_EpochCount--;
			}
		}
	}

	return 1;
}

/// Add bds ephemeris to navigation
int NavData::Nav_Add_BdsEph(BdsEph* eph)
{
	Nav_BdsEph* ptr = Nav_BDS_FindPrn(eph->m_prn, this);

	if (m_Nav_BdsEph == NULL || m_BDS_PRN_Count == 0)
	{
		bdstail = NULL;
	}

	if (ptr == NULL)
	{
		Nav_BdsEph* p = (Nav_BdsEph*)(malloc(sizeof(Nav_BdsEph)));
		memset((void*)(p), 0, sizeof(Nav_BdsEph));

		p->m_EpochCount++;
		p->m_BdsEph = eph;
		p->m_Prn = eph->m_prn;
		p->m_Tail = eph;

		m_BDS_PRN_Count++;
		if (m_Nav_BdsEph == NULL)
		{
			bdstail = p;
			m_Nav_BdsEph = p;
		}
		else
		{
			bdstail->m_Next = p;
			bdstail = p;
		}
	}
	else
	{
		BdsEph* pt = ptr->m_BdsEph;

		while (pt != NULL)
		{
			if (pt->m_Toc_Gps.m_Week == eph->m_Toc_Gps.m_Week &&
				(int)(pt->m_Toc_Gps.m_SecondofWeek) == (int)(eph->m_Toc_Gps.m_SecondofWeek) &&
				(int)(pt->m_IODE) == (int)(eph->m_IODE))
			{
				return 0;
			}
			pt = pt->m_Next;
		}

		ptr->m_EpochCount++;
		ptr->m_Tail->m_Next = eph;
		ptr->m_Tail = eph;

		if (m_Nav_Buff_Size > 0)
		{
			if (ptr->m_EpochCount > m_Nav_Buff_Size)
			{
				BdsEph* m_BdsEph_tmp = ptr->m_BdsEph->m_Next;

				if (ptr->m_BdsEph)
					free(ptr->m_BdsEph);

				ptr->m_BdsEph = m_BdsEph_tmp;
				ptr->m_EpochCount--;
			}
		}
	}

	return 1;
}

/// Read navigaton file body
int NavData::Read_NavigationFile_Body(char *filename)
{
	if (filename == NULL) return 0;

	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("Error: Navigation file open failed.\n");
		return 0;
	}

	unsigned int len = NAV_FILE_LINE + 1;
	char strline[NAV_FILE_LINE + 1] = "";

	char sel[21] = "END OF HEADER       ";
	if (SetPointerPos(fp, sel, strline, len) == 0)
	{
		printf("Warning: Navigation file open failed. Can not find 'END OF HEADER'.\n");
	}

	while (!is_eof(fp))
	{
		GetFileLine(fp, strline, len);

		char str_prn[4] = "";
		if (m_Header.m_Version + 0.1<3)
		{
			GetSubStr(strline, len, str_prn + 1, 3, 0, 2);
			str_prn[0] = m_Header.m_Sys;
		}
		else
		{
			GetSubStr(strline, len, str_prn, 4, 0, 3);
			if (str_prn[0] == ' ')
			{
				str_prn[0] = m_Header.m_Sys;
			}
		}

		unsigned int prn = GetSatNo(str_prn, 4, 0);

		int flag = m_Header.m_Version<3 ? 2 : 3;

		if (prn >= MIN_GPS_SATNO && prn<MAX_GPS_SATNO)
		{
			GpsEph *ptr = (GpsEph*)(malloc(sizeof(GpsEph)));
			if (ptr != NULL)
			{
				memset((void*)(ptr), 0, sizeof(GpsEph));
				ptr->m_prn = prn;
				ptr->Nav_ana_GpsEph(fp, strline, len, flag);
				Nav_Add_GpsEph(ptr);
			}
		}
		else if (prn >= MIN_GLO_SATNO && prn<MAX_GLO_SATNO)
		{
			GloEph *ptr = (GloEph*)(malloc(sizeof(GloEph)));
			if (ptr != NULL)
			{
				memset((void*)(ptr), 0, sizeof(GloEph));
				ptr->m_prn = prn;
				ptr->Nav_ana_GloEph(fp, strline, len, flag);
				Nav_Add_GloEph(ptr);
			}
		}
		else if (prn >= MIN_GAL_SATNO && prn<MAX_GAL_SATNO)
		{
			GalEph *ptr = (GalEph*)(malloc(sizeof(GalEph)));
			if (ptr != NULL)
			{
				memset((void*)(ptr), 0, sizeof(GalEph));
				ptr->m_prn = prn;
				ptr->Nav_ana_GalEph(fp, strline, len, flag);
				Nav_Add_GalEph(ptr);
			}
		}
		else if (prn >= MIN_BDS_SATNO && prn<MAX_BDS_SATNO)
		{
			BdsEph *ptr = (BdsEph*)(malloc(sizeof(BdsEph)));
			if (ptr != NULL)
			{
				memset((void*)(ptr), 0, sizeof(BdsEph));
				ptr->m_prn = prn;
				ptr->Nav_ana_BdsEph(fp, strline, len, flag);
				Nav_Add_BdsEph(ptr);
			}
		}
		else
		{
			// Other GNSS System Code
		}
	}

	fclose(fp);
	return 1;
}

/// Read navigation file
int NavData::Read_NavigationFile(char* filename)
{
	if (m_Header.Read_NavigationHeader(filename))
	{
		Read_NavigationFile_Body(filename);
	}
	else
		return 0;

	return 1;
}

/// Clear navigation data
int NavData::Clear()
{
	// clear gps data.
	Nav_GpsEph* ptr_gps = m_Nav_GpsEph;
	while (ptr_gps != NULL)
	{
		GpsEph* ptr = ptr_gps->m_GpsEph;
		while (ptr != NULL)
		{
			GpsEph* p = ptr->m_Next;
			free(ptr);
			ptr = p;
		}
		Nav_GpsEph* p = ptr_gps->m_Next;
		free(ptr_gps);
		ptr_gps = p;
	}

	// clear glo data.
	Nav_GloEph *ptr_glo = m_Nav_GloEph;
	while (ptr_glo != NULL)
	{
		GloEph* ptr = ptr_glo->m_GloEph;
		while (ptr != NULL)
		{
			GloEph* p = ptr->m_Next;
			free(ptr);
			ptr = p;
		}
		Nav_GloEph* p = ptr_glo->m_Next;
		free(ptr_glo);
		ptr_glo = p;
	}

	// clear bds data.
	Nav_BdsEph *ptr_bds = m_Nav_BdsEph;
	while (ptr_bds != NULL)
	{
		BdsEph* ptr = ptr_bds->m_BdsEph;
		while (ptr != NULL)
		{
			BdsEph* p = ptr->m_Next;
			free(ptr);
			ptr = p;
		}
		Nav_BdsEph* p = ptr_bds->m_Next;
		free(ptr_bds);
		ptr_bds = p;
	}

	// clear gal data.
	Nav_GalEph* ptr_gal = m_Nav_GalEph;
	while (ptr_gal != NULL)
	{
		GalEph* ptr = ptr_gal->m_GalEph;
		while (ptr != NULL)
		{
			GalEph* p = ptr->m_Next;
			free(ptr);
			ptr = p;
		}
		Nav_GalEph* p = ptr_gal->m_Next;
		free(ptr_gal);
		ptr_gal = p;
	}
	memset((void*)(this), 0, sizeof(NavData));
	return 1;
}


int PreFileHeader::Read_PreciseFileHeader(char* filename)
{
	FILE *fp = fopen(filename, "r");
	if (fp == NULL)	return 0;

	int len = PRECISE_LINE_CHARACTER_NUM + 1;//len = 101;
	char strline[PRECISE_LINE_CHARACTER_NUM + 1] = "";//strline[101] = "";

	char version[2] = "", type[2] = "";
	GetFileLine(fp, strline, len);
	GetSubStr(strline, len, version, 2, 1, 1);
	m_Version = version[0];
	GetSubStr(strline, len, type, 2, 2, 1);
	m_FileType = type[0];

	char year[5] = "", month[3] = "", day[3] = "";
	char hour[3] = "", minute[3] = "", second[12] = "";

	GetSubStr(strline, len, year, 5, 3, 4);
	m_FirstEpoch.m_Year = str_to_i(year);
	GetSubStr(strline, len, month, 3, 8, 2);
	m_FirstEpoch.m_Month = str_to_i(month);
	GetSubStr(strline, len, day, 3, 11, 2);
	m_FirstEpoch.m_Day = str_to_i(day);

	GetSubStr(strline, len, hour, 3, 14, 2);
	m_FirstEpoch.m_Hour = str_to_i(hour);
	GetSubStr(strline, len, minute, 3, 17, 2);
	m_FirstEpoch.m_Min = str_to_i(minute);
	GetSubStr(strline, len, second, 12, 20, 11);
	m_FirstEpoch.m_Sec = str_to_i(second);

	m_FirstEpoch._date2jd();

	char count[8] = "";
	GetSubStr(strline, len, count, 8, 32, 7);
	m_TotalEpochNum = str_to_i(count);

	GetSubStr(strline, len, m_DataType, 6, 40, 5);
	GetSubStr(strline, len, m_CoordSystem, 6, 46, 5);

	GetSubStr(strline, len, m_OrbitType, 4, 52, 3);
	GetSubStr(strline, len, m_Organization, 5, 56, 4);

	// read second line.
	GetFileLine(fp, strline, len);
	char week[5] = "", sec[16] = "";

	GetSubStr(strline, len, week, 5, 3, 4);
	m_FirstEpoch_Week = (int)str_to_i(week);
	GetSubStr(strline, len, sec, 16, 8, 15);
	m_FirstEpoch_Sec = str_to_f(sec);

	char interval[15] = "", mjd_u[6] = "", mjd_m[16] = "";

	GetSubStr(strline, len, interval, 15, 24, 14);
	m_Interval = str_to_f(interval);
	GetSubStr(strline, len, mjd_u, 6, 39, 5);
	m_Mjd_u = str_to_f(mjd_u);
	GetSubStr(strline, len, mjd_m, 16, 45, 15);
	m_Mjd_m = str_to_f(mjd_m);

	// read third line.
	GetFileLine(fp, strline, len);
	char satcount[4] = "";
	GetSubStr(strline, len, satcount, 4, 3, 3);
	m_StaCount = (int)str_to_i(satcount);

	// read line 3-7.
	unsigned int i = 0, j = 0;
	unsigned int line_num = 3;

	for (i = 0; i< m_StaCount; ++i)
	{
		if (line_num == 11)//8
		{
			break;
		}
		char prn[4] = "";
		GetSubStr(strline, len, prn, 4, 9 + j * 3, 3);
		m_Satllite[i].m_Prn = GetSatNo(prn, 4, 0);
		if ((++j) == 17)
		{
			j = 0;
		}
		if ((((i + 1) % 17) == 0) && ((i + 1)< m_StaCount))
		{
			GetFileLine(fp, strline, len);
			++line_num;
		}
	}
	for (i = 0; i<11 - line_num; ++i)//8
	{
		GetFileLine(fp, strline, len);
	}

	//   read line 8-22
	line_num = 11;
	j = 0;
	for (i = 0; i< m_StaCount; ++i)
	{
		if (line_num == 16)//13
		{
			break;
		}
		char pre[4] = "";
		GetSubStr(strline, len, pre, 4, 9 + j * 3, 3);
		m_Satllite[i].m_Precise = str_to_i(pre);
		if ((++j) == 17)
		{
			j = 0;
		}
		if ((((i + 1) % 17) == 0) && ((i + 1)< m_StaCount))
		{
			GetFileLine(fp, strline, len);
			++line_num;
		}
	}
	for (i = 0; i<22 - line_num; ++i)
	{
		GetFileLine(fp, strline, len);
	}

	fclose(fp);
	return 1;
}

unsigned int PreEpochData::Pre_ana_Stringline_p(char *strline, unsigned int len)
{
	char prn_str[4] = "";
	GetSubStr(strline, len, prn_str, 4, 1, 3);
	unsigned int prn = GetSatNo(prn_str, 4, 0);

	char x[15] = "", y[15] = "", z[15] = "", cor[15] = "";
	GetSubStr(strline, len, x, 15, 4, 14);
	m_X = str_to_f(x)*1000.0;
	GetSubStr(strline, len, y, 15, 18, 14);
	m_Y = str_to_f(y)*1000.0;
	GetSubStr(strline, len, z, 15, 32, 14);
	m_Z = str_to_f(z)*1000.0;
	GetSubStr(strline, len, cor, 15, 46, 14);

	double ds = sqrt(m_X* m_X + m_Y* m_Y + m_Z* m_Z);

	if (ds< 10000) return 0;

	double tmp_cor = str_to_f(cor);
	if (tmp_cor>500)
	{
		m_ClockCorr = 999999;
	}
	else
	{
		m_ClockCorr = tmp_cor*1.0e-6;
	}

	return prn;
}

PreciseData::PreciseData()
{
	m_PrnCount = 0;
	m_Predata = NULL;
	memset((void*)(&m_Header), 0, sizeof(PreFileHeader));
	m_NavData = (NavData*)(calloc(1, sizeof(NavData)));
	//m_NavData = NULL;
}
PreciseData::~PreciseData()
{
	this->Clear();
}

PreDataSat* PreciseData::Precise_Search_Prn(unsigned int prn)
{
	PreDataSat *p = m_Predata;
	while (p != NULL)
	{
		if (p->m_Prn == prn)
		{
			break;
		}
		p = p->m_Next;
	}
	return p;
}

PreEpochData* PreciseData::Precise_Search_Nearest_Epoch(unsigned int prn, gnsstime& time)
{
	if (time.m_Jd == 0)
		time._date2jd();

	PreDataSat * sat = Precise_Search_Prn(prn);
	if (sat != NULL)
	{
		if (sat->m_FirstEpoch != NULL)
		{
			PreEpochData *rigt_t = NULL;
			PreEpochData *left_t = NULL;
			PreEpochData *end_t = NULL;
			PreEpochData *ref_t = sat->m_MidEpoch;

			if (time.m_Jd>ref_t->m_GnssTime.m_Jd)
			{
				rigt_t = ref_t;
				end_t = NULL;
			}
			else
			{
				rigt_t = sat->m_FirstEpoch;
				end_t = ref_t->m_Next;
			}

			while (rigt_t->m_Next != end_t)
			{
				if (time.m_Jd<rigt_t->m_GnssTime.m_Jd)
				{
					break;
				}
				rigt_t = rigt_t->m_Next;
			}

			left_t = rigt_t->m_Previous;

			if (left_t == NULL)
			{
				return rigt_t;
			}
			else
			{
				double f1 = fabs(left_t->m_GnssTime.m_Jd - time.m_Jd);
				double f2 = fabs(rigt_t->m_GnssTime.m_Jd - time.m_Jd);

				return f1<f2 ? left_t : rigt_t;
			}
		}
		else
		{
			return NULL;
		}
	}
	else
	{
		return NULL;
	}
}

PreEpochData* PreciseData::Precise_Search_Nearest_Epoch_w(unsigned int prn, gnsstime& time)
{
	if (time.m_Jd == 0)
		time._date2jd();

	PreDataSat * sat = Precise_Search_Prn(prn);
	if (sat != NULL)
	{
		if (sat->m_FirstEpoch != NULL)
		{
			if (sat->m_IndexTable == NULL)
			{
				return Precise_Search_Nearest_Epoch(prn, time);
			}
			int i = 0, j = sat->m_EpochCount - 1;
			int d = abs(j - i);

			while (d>1)
			{
				unsigned int k = d / 2 + i;
				if (time.m_Jd<sat->m_IndexTable[k]->m_GnssTime.m_Jd)
				{
					j = k;
				}
				else
				{
					i = k;
				}

				d = abs(j - i);
			}

			if (d == 0)
			{
				return sat->m_IndexTable[i];
			}
			else
			{
				double f1 = fabs(sat->m_IndexTable[i]->m_GnssTime.m_Jd - time.m_Jd);
				double f2 = fabs(sat->m_IndexTable[j]->m_GnssTime.m_Jd - time.m_Jd);

				return f1<f2 ? sat->m_IndexTable[i] : sat->m_IndexTable[j];
			}
		}
		else
		{
			return NULL;
		}
	}
	else
	{
		return NULL;
	}
}


PreEpochData* PreciseData::Get_Precise_InterpData(unsigned int prn, gnsstime& time, int& len, int n)
{
	len = 0;
	int m = n + 1;
	PreEpochData * ptr_t = Precise_Search_Nearest_Epoch_w(prn, time);
	if (ptr_t == NULL)
	{
		return NULL;
	}

	PreEpochData *ptr_d = (PreEpochData *)(malloc(sizeof(PreEpochData)*m));
	if (ptr_d == NULL)
	{
		return NULL;
	}
	memset((void*)(ptr_d), 0, sizeof(PreEpochData)*m);

	PreEpochData *ptr_left = ptr_t;
	PreEpochData *ptr_rigt = ptr_t;

	len += 1;
	while (n)
	{
		if (ptr_left->m_Previous != NULL)
		{
			ptr_left = ptr_left->m_Previous;
			len += 1;
			--n;
		}

		if (n == 0)
		{
			break;
		}

		if (ptr_rigt->m_Next != NULL)
		{
			ptr_rigt = ptr_rigt->m_Next;
			len += 1;
			--n;
		}
	}

	double tolre = m_Header.m_Interval*(m - 1)*1.5;
	if (ptr_rigt->m_GnssTime.m_Jd - ptr_left->m_GnssTime.m_Jd>tolre)
	{
		len = 0;
		free(ptr_d);
		ptr_d = NULL;
		return ptr_d;
	}

	int i = 0;
	for (i = 0; i<len; ++i)
	{
		ptr_d[i] = *ptr_left;
		ptr_left = ptr_left->m_Next;
	}

	return ptr_d;
}

/// Search precise epoch ephemeris by time
PreEpochData*  PreciseData::Precise_Search_by_Epoch(unsigned int prn, gnsstime &time)
{
	PreEpochData *rsat = NULL;
	PreDataSat *sat = Precise_Search_Prn(prn);
	if (sat != NULL)
	{
		rsat = sat->m_FirstEpoch;
		while (rsat)
		{
			if (rsat->m_GnssTime == time)
				return rsat;

			rsat = rsat->m_Next;
		}
	}

	return NULL;
}

int PreciseData::Precise_Interp_Pos_w(unsigned int prn, gnsstime& time, double pos[3], int n)
{
	int len = 0;
	PreEpochData *ptr = Get_Precise_InterpData(prn, time, len, n);

	int re = 1;
	re = Precise_Interp_Pos(time, ptr, len, pos, n);
	free(ptr);
	ptr = NULL;

	return re;
}

int PreciseData::Precise_Interp_Vel_w(unsigned int prn, gnsstime& time, double vel[3], int n)
{
	int len = 0;
	PreEpochData *ptr = Get_Precise_InterpData(prn, time, len, n);

	int re = 1;
	re = Precise_Interp_Vel(time, ptr, len, vel, n);
	free(ptr);
	ptr = NULL;

	return re;
}

int PreciseData::Precise_Interp_PosVel_w(unsigned int prn, gnsstime& time, double pos[3], double vel[3], int n)
{
	int len = 0;
	PreEpochData *ptr = Get_Precise_InterpData(prn, time, len, n);

	int re = 1;
	re = Precise_Interp_Pos(time, ptr, len, vel, n);
	re = Precise_Interp_Vel(time, ptr, len, pos, n);
	free(ptr);
	ptr = NULL;

	return re;
}

int PreciseData::Precise_Add_Data_p(unsigned int prn, PreEpochData *epochdata, PreDataSat **tail_sat)
{
	if (prn == 0)
		return 0;

	PreDataSat *p = this->Precise_Search_Prn(prn);
	if (m_Predata == NULL || p == NULL)
	{
		PreDataSat *sat_data = (PreDataSat *)(malloc(sizeof(PreDataSat)));
		if (sat_data == NULL)
		{
			return 0;
		}
		memset(sat_data, 0, sizeof(PreDataSat));
		sat_data->m_Prn = prn;
		sat_data->m_EpochCount = 0;
		if (m_Predata == NULL)
		{
			m_Predata = sat_data;
			*tail_sat = sat_data;
		}
		else
		{
			(*tail_sat)->m_Next = sat_data;
			(*tail_sat) = sat_data;
		}
		m_PrnCount += 1;
		p = sat_data;
	}

	if (p->m_FirstEpoch == NULL)
	{
		p->m_FirstEpoch = epochdata;
		p->LastEpoch = epochdata;
	}
	else
	{
		epochdata->m_Previous = p->LastEpoch;
		p->LastEpoch->m_Next = epochdata;
		p->LastEpoch = epochdata;
	}

	p->m_EpochCount += 1;

	return 1;
}

/// Read precise ephemeris body
int PreciseData::Read_PreciseFileBody(char* filename)
{
	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
		return 0;

	unsigned int len = PRECISE_LINE_CHARACTER_NUM + 1;
	char strline[PRECISE_LINE_CHARACTER_NUM + 1] = "";

	if (strlen(filename)<3)
	{
		fclose(fp);
		return 0;
	}

	int is_igg = 0;
	unsigned int skip_line = 22;
	if (!strcmp(m_Header.m_Organization, (char*)(" IGG")) || !strcmp(m_Header.m_Organization, (char*)(" igg")))
	{
		skip_line = 32;
		is_igg = 1;
	}

	unsigned int i = 0;
	for (; i<skip_line; ++i)
	{
		GetFileLine(fp, strline, len);
	}

	unsigned int flags = 0;
	gnsstime epochtime;

	PreDataSat *tail_sat = NULL;
	if (m_PrnCount == 0)
	{
		tail_sat = NULL;
		m_Predata = NULL;
	}
	else
	{
		PreDataSat *ptr = m_Predata;
		while (ptr != NULL)
		{
			tail_sat = ptr;
			ptr = ptr->m_Next;
		}
	}

	while (flags = Pre_GetLine(fp, strline, len))
	{
		if (flags == 1)       // Time information
		{
			Pre_ana_Epochtime(epochtime, strline, len);
		}
		else if (flags == 4)  // unknown information. ignore this line
		{
			continue;
		}
		else                  // position and velocity information
		{
 			PreEpochData* epochdata = (PreEpochData*)(malloc(sizeof(PreEpochData)));
			if (epochdata == NULL)
			{
				continue;
			}
			memset((void*)(epochdata), 0, sizeof(PreEpochData));

			if (flags == 2)   // position information
			{
				unsigned int prn = epochdata->Pre_ana_Stringline_p(strline, len);
				epochdata->m_GnssTime = epochtime;
				if (this->Precise_Add_Data_p(prn, epochdata, &tail_sat) == 0) free(epochdata);
			}
			else
			{
				// flags=3. velocity information.
			}
		}
	}

	PreDataSat *p = m_Predata;
	while (p != NULL)
	{
		if (p->m_IndexTable != NULL)
		{
			free(p->m_IndexTable);
			p->m_IndexTable = NULL;
		}
		p->m_IndexTable = (PreEpochdata_ptr *)(malloc(p->m_EpochCount*sizeof(PreEpochdata_ptr)));
		memset((void *)(p->m_IndexTable), 0, p->m_EpochCount*sizeof(PreEpochdata_ptr));
		if (p->m_IndexTable == NULL)
		{
			p->m_MidEpoch = p->m_FirstEpoch;
			Pre_EpochData_ptr_Advance(&(p->m_MidEpoch), p->m_EpochCount / 2);
		}
		else
		{
			p->m_MidEpoch = p->m_FirstEpoch;
			PreEpochData *ptr = p->m_FirstEpoch;
			unsigned int i = 0;
			while (ptr != NULL)
			{
				p->m_IndexTable[i] = ptr;
				++i;
				ptr = ptr->m_Next;
			}
			p->m_MidEpoch = p->m_IndexTable[p->m_EpochCount / 2];
		}
		p = p->m_Next;
	}

	fclose(fp);
	return 1;
}

/// Read precise data --------------------------------------------------------------
/* Parameter:
* char *filename            I            file name
* int   flag                I            read file flag:1 only pre data  0:nav+pre
* Return:                                0: read faild
*                                        1: successfully
* --------------------------------------------------------------------------------- */
int PreciseData::Read_PreciseFile(char* filename, int flag)
{
	// Read navigation data
	if (flag == 0)
	{
		NavData *nav = NULL;
		if (m_NavData == NULL)
		{
			nav = (NavData *)(malloc(sizeof(NavData)));
			memset((void*)(nav), 0, sizeof(NavData));

			if (nav->Read_NavigationFile(filename))
			{
				m_NavData = nav;
			}
			else
			{
				free(nav); return 0;
			}
		}
		else
		{
			if (m_NavData->Read_NavigationFile(filename))
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
		return 1;
	}

	// read precise ephemeris
	if (m_Header.Read_PreciseFileHeader(filename) == 0)
	{
		return 0;
	}
	if (Read_PreciseFileBody(filename) == 0)
	{
		memset((void*)(this), 0, sizeof(PreciseData));
		return 0;
	}

	return 1;
}

int PreciseData::Extract_PreciseClock(char* filename)
{
	int n = m_PrnCount;
	int i = 0, j = 0;

	FILE *fp = fopen(filename, "r");
	if (fp == NULL)	return 0;

	fprintf(fp, "%         s\n", "     3.00           C                                       RINEX VERSION / TYPE");
	fprintf(fp, "%s\n", "                                                            END OF HEADER       ");

	PreDataSat *sat_dat = m_Predata;

	while (sat_dat != NULL)
	{
		unsigned int prn = sat_dat->m_Prn;
		PreEpochData *p = sat_dat->m_FirstEpoch;
		while (p != NULL)
		{
			double cor = p->m_ClockCorr;
			if (cor>1000)
			{
				p = p->m_Next;
				continue;
			}
			char strline[100] = "";

			sprintf(strline, "AS %03d  %04d %02d %02d %02d %02d %9.6f  1",
				prn, p->m_GnssTime.m_Year, p->m_GnssTime.m_Month,
				p->m_GnssTime.m_Day, p->m_GnssTime.m_Hour, p->m_GnssTime.m_Min, p->m_GnssTime.m_Sec);

			char flag = GetSystem(prn);
			strline[3] = flag;
			fprintf(fp, "%s", strline);

			char str_cor[30] = "";

			sprintf(str_cor, "%23.12e", cor);

			int i = 0, len = strlen(str_cor);
			for (; i<len; ++i)
			{
				if (str_cor[i] == 'e')
				{
					break;
				}
			}
			i += 2;
			int j = i + 1;
			int k = i + 2;
			str_cor[i] = str_cor[j];
			str_cor[j] = str_cor[k];
			str_cor[k] = '\0';

			fprintf(fp, "%s\n", str_cor);

			p = p->m_Next;
		}
		sat_dat = sat_dat->m_Next;
	}

	fclose(fp);
	return 1;
}

int PreciseData::Clear()
{
	memset((void*)(&m_Header), 0, sizeof(PreFileHeader));

	m_PrnCount = 0;
	PreDataSat *p = m_Predata;
	while (p != NULL)
	{
		PreEpochData *eptr = p->m_FirstEpoch;
		while (eptr != NULL)
		{
			PreEpochData *eptr_t = eptr->m_Next;
			free(eptr);
			eptr = eptr_t;
		}

		PreDataSat *ptr = p->m_Next;
		free(p->m_IndexTable);
		free(p);
		p = ptr;
	}

	if (m_NavData)
	{
		m_NavData->Clear();
		free(m_NavData);
	}

	return 1;
}

/// Get Glonass Slot No.
int Get_Nav_Orbitn(unsigned int prn, NavData* navdata)
{
	if (navdata == NULL) return -1;

	if (GetSystem(prn) != 'R') return 0;

	Nav_GloEph* PtrGlo = navdata->m_Nav_GloEph;

	while (PtrGlo != NULL)
	{
		if (PtrGlo->m_Prn == prn)
			break;
		PtrGlo = PtrGlo->m_Next;
	}

	if (PtrGlo == NULL) return -1;

	GloEph* p = PtrGlo->m_GloEph;

	if (p == NULL) return -1;

	return (int)(p->m_frequency_number);
}
