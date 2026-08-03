#include "SUPREME_ClkRNX.h"
#include "SUPREME_CommonFunction.h"

#define CLOCK_LINE_CHARACTER_NUM    100

int Clock_Search_As(FILE* fp, char* strline, unsigned int len)
{
	char strflag[] = "AS";
	char stras[3] = "";
	do
	{
		GetFileLine(fp, strline, len);
		GetSubStr(strline, len, stras, 3, 0, 2);
		if (is_eof(fp))                             // comes to end of file.
		{
			return 0;                               // do not find the flag. ma be file is not integrity.
		}
	} while (strcmp(stras, strflag));

	return 1;                                       // find the flag, and "string line" save a line data.
}

int Clock_EpochData_ptr_Advance(Clk_EpochData **ptr, int count)
{
	int n = abs(count);
	int i = 0;
	if (count>0)
	{
		for (i = 0; i<n; ++i)
		{
			if (*ptr == NULL)
			{
				break;
			}
			*ptr = (*ptr)->m_next;
		}
	}
	else
	{
		for (i = 0; i<n; ++i)
		{
			if (*ptr == NULL)
			{
				break;
			}
			*ptr = (*ptr)->m_previous;
		}
	}

	return 1;
}

int Clock_EpochData_ptr_Disdance(Clk_EpochData *ptr, Clk_EpochData *des)
{
	int d_t = 0;
	while (ptr != des)
	{
		ptr = ptr->m_next;
		++d_t;
	}
	return d_t;
}

double Clock_Interpolate(gnsstime& time, Clk_EpochData *ptr, int len, int n)//卫星钟差内插
{
	if (ptr == NULL || len == 0 || len != (n + 1))
	{
		return 999999;
	}
	if (time.m_Jd<ptr[0].m_time.m_Jd || time.m_Jd>ptr[len - 1].m_time.m_Jd)
	{
		return 999999;
	}

	double result = 0.0;

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
		result += (temp*ptr[i].m_clockcor);
	}

	return result;
}

double Clock_Interpolate_w(unsigned int prn, gnsstime& time, ClockData *clkdata, int n)
{
	int len = 0;
	Clk_EpochData *ptr = clkdata->Get_Clock_Interp_Data(prn, time, len, n);

	double result = 0.0;
	result = Clock_Interpolate(time, ptr, len, n);
	free(ptr);
	ptr = NULL;

	return result;
}

int Clk_File_Header::Read_Clk_FileHeader(char* filename,double &verson_)
{
	FILE *fp = fopen(filename, "r");
	if (fp == NULL)	return 0;

	unsigned int len = CLOCK_LINE_CHARACTER_NUM + 1;
	char strline[CLOCK_LINE_CHARACTER_NUM + 1] = "";
	char sel1[21] = "RINEX VERSION / TYPE";
	GetFileLine(fp, strline, len);
	char version[21] = "";
	GetSubStr(strline, len, version, 21, 0, 20);
	m_Version = str_to_f(version);
	if (m_Version == 3.04) 
	{
		char tempstr[21] = "";
		GetSubStr(strline, len, tempstr, 21, 65, 20);//原来60；
		if (strcmp(tempstr, sel1) == 0)
		{
			char filetype = strline[21];//yunali  20
			if ((filetype != 'c') && (filetype != 'C'))
			{
				fclose(fp);
				return 0;
			}
			m_FileType = filetype;
			
		}
	}
	else if (m_Version == 3.00)
	{
		char tempstr[21] = "";
		GetSubStr(strline, len, tempstr, 21, 60, 20);//原来60；
		if (strcmp(tempstr, sel1) == 0)
		{
			char filetype = strline[20];//yunali  20
			if ((filetype != 'c') && (filetype != 'C'))
			{
				fclose(fp);
				return 0;
			}
			m_FileType = filetype;
		}	
	}
	else
	{
		fclose(fp);
		return 0;
	}
	char sel2[21] = "LEAP SECONDS        ";
	char se13[21] = "LEAP SECONDS GNSS   ";
	char se14[21] = " ";
	if (m_Version == 3.04) 
	{
		for (int i = 0; i < 21; i++)
		{
			se14[i] = se13[i];
		}
	}
	else 
	{
		for (int i = 0; i < 21; i++)
		{
			se14[i] = sel2[i];
		}
	}
	if (SetPointerPos_(fp, se14, strline, len,m_Version))      // find the string.
	{
		char leap[21] = "";
		GetSubStr(strline, len, leap, 21, 0, 20);
		m_LeapSec = str_to_i(leap);                 // get leap seconds.
	}
	verson_ = m_Version;
	fclose(fp);
	return 1;
}

unsigned int Clk_EpochData::Clk_ana_StringLine(char* strline, unsigned int len,double verson)
{
	char str_prn[4] = "";
	unsigned int prn = 0;

	GetSubStr(strline, len, str_prn, 4, 3, 3);
	prn = GetSatNo(str_prn, 4, 0);// satellite index
	int lenth = 0;
	char year[5] = "", month[3] = "", day[3] = "";
	char hour[3] = "", minute[3] = "", sec[10] = "";
	if (verson == 3.04)lenth = 5;
	GetSubStr(strline, len, year, 5, 8+lenth, 4);
	m_time.m_Year = str_to_i(year);
	GetSubStr(strline, len, month, 3, 13 + lenth, 2);
	m_time.m_Month = str_to_i(month);
	GetSubStr(strline, len, day, 3, 16 + lenth, 2);
	m_time.m_Day = str_to_i(day);
	GetSubStr(strline, len, hour, 3, 19 + lenth, 2);
	m_time.m_Hour = str_to_i(hour);
	GetSubStr(strline, len, minute, 3, 22 + lenth, 2);
	m_time.m_Min = str_to_i(minute);
	GetSubStr(strline, len, sec, 10, 25 + lenth, 9);
	m_time.m_Sec = str_to_f(sec);

	m_time._date2jd();
	m_time._date2gpst();

	char str_num[4] = "";
	GetSubStr(strline, len, str_num, 4, 34 + lenth, 3);
	int num = str_to_i(str_num);

	char str_clockcor[19] = "", str_e[4] = "";
	GetSubStr(strline, len, str_clockcor, 19, 37 + lenth, 18);
	GetSubStr(strline, len, str_e, 4, 56 + lenth, 3);
	m_clockcor = str_to_f(str_clockcor);
	m_clockcor *= pow((double)(10.0), str_to_f(str_e));

	if (num>1)
	{
		char str_clockvel[17] = "", str_e[4] = "";
		GetSubStr(strline, len, str_clockvel, 17, 59 + lenth, 16);
		GetSubStr(strline, len, str_e, 4, 76 + lenth, 3);
		m_clockvel = str_to_f(str_clockvel);
		m_clockvel *= pow((double)(10.0), str_to_f(str_e));
	}
	return prn;
}

ClockData::ClockData()
{
	memset((void*)this, 0, sizeof(ClockData));
}

ClockData::~ClockData()
{
	Clear();
}

Clk_Data_Satellite* ClockData::Search_clock_Prn(unsigned int prn)
{
	Clk_Data_Satellite *p = m_clkdata;
	while (p != NULL)
	{
		if (p->m_prn == prn)
		{
			break;
		}
		p = p->m_next;
	}
	return p;
}

Clk_EpochData* ClockData::Search_Nearest_clock_Epoch(unsigned int prn, gnsstime& time)
{
	if (time.m_Jd == 0)
	{
		time._date2jd();
	}

	Clk_Data_Satellite * sat = Search_clock_Prn(prn);
	if (sat != NULL)
	{
		if (sat->m_FirstEpoch != NULL)
		{
			Clk_EpochData *rigt_t = NULL;
			Clk_EpochData *left_t = NULL;
			Clk_EpochData *end_t = NULL;
			Clk_EpochData *ref_t = sat->m_MidEpoch;

			if (time.m_Jd>ref_t->m_time.m_Jd)
			{
				rigt_t = ref_t;
				end_t = NULL;
			}
			else
			{
				rigt_t = sat->m_FirstEpoch;
				end_t = ref_t->m_next;
			}

			while (rigt_t->m_next != end_t)
			{
				if (time.m_Jd<rigt_t->m_time.m_Jd)
				{
					break;
				}
				rigt_t = rigt_t->m_next;
			}

			left_t = rigt_t->m_previous;

			if (left_t == NULL)
			{
				return rigt_t;
			}
			else
			{
				double f1 = fabs(left_t->m_time.m_Jd - time.m_Jd);
				double f2 = fabs(rigt_t->m_time.m_Jd - time.m_Jd);

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

Clk_EpochData* ClockData::Search_Nearest_clock_Epoch_w(unsigned int prn, gnsstime& time)
{
	if (time.m_Jd == 0)
	{
		time._date2jd();
	}
	Clk_Data_Satellite * sat = Search_clock_Prn(prn);
	if (sat != NULL)
	{
		if (sat->m_FirstEpoch != NULL)
		{
			if (sat->m_Index_Table == NULL)
			{
				return Search_Nearest_clock_Epoch(prn, time);
			}
			int i = 0, j = sat->m_epochcount - 1;
			int d = abs(j - i);

			while (d > 1)
			{
				unsigned int k = d / 2 + i;
				if (time.m_Jd < sat->m_Index_Table[k]->m_time.m_Jd)
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
				return sat->m_Index_Table[i];
			}
			else
			{
				double f1 = fabs(sat->m_Index_Table[i]->m_time.m_Jd - time.m_Jd);
				double f2 = fabs(sat->m_Index_Table[j]->m_time.m_Jd - time.m_Jd);

				return f1 < f2 ? sat->m_Index_Table[i] : sat->m_Index_Table[j];
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

Clk_EpochData* ClockData::Get_Clock_Interp_Data(unsigned int prn, gnsstime& time, int& len, int n)
{
	len = 0;
	int m = n + 1;

	Clk_EpochData * ptr_t = Search_Nearest_clock_Epoch_w(prn, time);
	if (ptr_t == NULL)
	{
		return NULL;
	}

	Clk_EpochData *ptr_d = (Clk_EpochData *)(malloc(sizeof(Clk_EpochData)*m));
	if (ptr_d == NULL)
	{
		return NULL;
	}
	memset((void*)(ptr_d), 0, sizeof(Clk_EpochData)*m);

	Clk_EpochData *ptr_left = ptr_t;
	Clk_EpochData *ptr_rigt = ptr_t;

	len += 1;
	while (n)
	{
		if (ptr_left->m_previous != NULL)
		{
			ptr_left = ptr_left->m_previous;
			len += 1;
			--n;
		}

		if (n == 0)
		{
			break;
		}

		if (ptr_rigt->m_next != NULL)
		{
			ptr_rigt = ptr_rigt->m_next;
			len += 1;
			--n;
		}
	}

	if (ptr_rigt->m_time.m_Jd - ptr_left->m_time.m_Jd>3600)
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
		ptr_left = ptr_left->m_next;
	}

	return ptr_d;
}

Clk_EpochData* ClockData::Search_clock_by_Epoch(unsigned int prn, gnsstime& time)
{
	Clk_EpochData *rsat = NULL;
	Clk_Data_Satellite *sat = Search_clock_Prn(prn);

	if (sat != NULL)
	{
		rsat = sat->m_FirstEpoch;

		while (rsat)
		{
			if (rsat->m_time == time)
				return rsat;

			rsat = rsat->m_next;
		}
	}

	return NULL;
}

int ClockData::Clock_Add_Data(unsigned int prn, Clk_EpochData* epochdata, Clk_Data_Satellite **tail_sat)
{
	Clk_Data_Satellite *p = this->Search_clock_Prn(prn);
	if (m_clkdata == NULL || p == NULL)
	{
		Clk_Data_Satellite *sat_data = (Clk_Data_Satellite *)(malloc(sizeof(Clk_Data_Satellite)));
		if (sat_data == NULL)
		{
			return 0;
		}
		memset(sat_data, 0, sizeof(Clk_Data_Satellite));

		sat_data->m_prn = prn;
		sat_data->m_epochcount = 0;

		if (m_clkdata == NULL)
		{
			m_clkdata = sat_data;
			*tail_sat = sat_data;
		}
		else
		{
			(*tail_sat)->m_next = sat_data;
			(*tail_sat) = sat_data;
		}
		m_PrnCount += 1;
		p = sat_data;
	}
	if (p->m_FirstEpoch == NULL)
	{
		p->m_FirstEpoch = epochdata;
		p->m_LastEpoch = epochdata;
	}
	else
	{
		epochdata->m_previous = p->m_LastEpoch;
		p->m_LastEpoch->m_next = epochdata;
		p->m_LastEpoch = epochdata;
	}
	p->m_epochcount += 1;

	return 1;
}

/// Read clock file body
int ClockData::Read_ClkFileBody(char* filename,double verson)
{
	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
		return 0;

	unsigned int len = CLOCK_LINE_CHARACTER_NUM + 1;
	char strline[CLOCK_LINE_CHARACTER_NUM + 1] = "";

	// skip the header.  if can not find "sel",close the file.
	char sel[21] = "END OF HEADER       ";
	if (SetPointerPos_(fp, sel, strline, len,verson) == 0) // find the string.
	{
		fclose(fp);
	}

	Clk_Data_Satellite * tail_sat = NULL;
	if (m_PrnCount == 0)
	{
		tail_sat = NULL;
		m_clkdata = NULL;
	}
	else
	{
		Clk_Data_Satellite *ptr = m_clkdata;
		while (ptr != NULL)
		{
			tail_sat = ptr;
			ptr = ptr->m_next;
		}
	}
	while (Clock_Search_As(fp, strline, len))
	{
		Clk_EpochData *epochdata = (Clk_EpochData *)(malloc(sizeof(Clk_EpochData)));
		if (epochdata == NULL)
		{
			continue;
		}
		memset((void *)(epochdata), 0, sizeof(Clk_EpochData));

		unsigned int prn = epochdata->Clk_ana_StringLine(strline, len,verson);//from string line to struct.
		if (prn == 0)
		{
			free(epochdata);
			continue;
		}
		this->Clock_Add_Data(prn, epochdata, &tail_sat);
	}
	Clk_Data_Satellite *p = m_clkdata;

	while (p != NULL)
	{
		if (p->m_Index_Table != NULL)
		{
			free(p->m_Index_Table);
			p->m_Index_Table = NULL;
		}
		p->m_Index_Table = (Clk_EpochData_ptr *)(malloc(p->m_epochcount*sizeof(Clk_EpochData_ptr)));
		memset((void *)(p->m_Index_Table), 0, p->m_epochcount*sizeof(Clk_EpochData_ptr));
		if (p->m_Index_Table == NULL)
		{
			p->m_MidEpoch = p->m_FirstEpoch;
			Clock_EpochData_ptr_Advance(&(p->m_MidEpoch), p->m_epochcount / 2);
		}
		else
		{
			p->m_MidEpoch = p->m_FirstEpoch;
			Clk_EpochData *ptr = p->m_FirstEpoch;
			unsigned int i = 0;
			while (ptr != NULL)
			{
				p->m_Index_Table[i] = ptr;
				++i;
				ptr = ptr->m_next;
			}
			p->m_MidEpoch = p->m_Index_Table[p->m_epochcount / 2];
		}
		p = p->m_next;
	}

	fclose(fp);
	return 1;
}

/// Read precise clock data
int ClockData::Read_ClockFile(char* filename)
{
	double verson = 0;

	if (m_Header.Read_Clk_FileHeader(filename,verson) == 0)
		return 0;

	if (Read_ClkFileBody(filename,verson) == 0)
	{
		memset((void*)(this), 0, sizeof(ClockData));
		return 0;
	}

	return 1;
}

int ClockData::Clear()
{
	memset((void*)(&m_Header), 0, sizeof(Clk_File_Header));
	m_PrnCount = 0;
	Clk_Data_Satellite *p = m_clkdata;
	while (p != NULL)
	{
		Clk_EpochData *eptr = p->m_FirstEpoch;
		while (eptr != NULL)
		{
			Clk_EpochData *eptr_t = eptr->m_next;
			free(eptr);
			eptr = eptr_t;
		}
		Clk_Data_Satellite *ptr = p->m_next;
		free(p->m_Index_Table);
		free(p);
		p = ptr;
	}
	m_clkdata = NULL;

	return 1;
}

/// Clock data compare
int Clock_Compare(char* base_file, char* compare_file, char* result_file, char svn[4], int n)
{
	if (base_file == NULL || compare_file == NULL || result_file == NULL)
	{
		return 0;
	}
	ClockData base_data;
	if (base_data.Read_ClockFile(base_file) == 0)
	{
		return 0;
	}
	ClockData compare_data;
	if (compare_data.Read_ClockFile(compare_file) == 0)
	{
		return 0;
	}
	unsigned int prn = GetSatNo(svn, 4, 0);

	Clk_Data_Satellite * sat = base_data.Search_clock_Prn(prn);
	if (sat == NULL)
	{
		printf("cant not find %s in base file.\n", svn);
		return 0;
	}

	FILE *fp = fopen(result_file, "w");
	if (fp == NULL)
	{
		printf("can not open file: %s\n", result_file);
		return 0;
	}

	Clk_EpochData *p = sat->m_FirstEpoch;

	int clk06 = 0, clk12 = 0, clk18 = 0;

	while (p != NULL)
	{
		if (p->m_time.m_Hour >= 6 && clk06 == 0)
		{
			clk06 = 1;
			int len = strlen(compare_file);
			compare_file[len - 6] = '0';
			compare_file[len - 5] = '6';
			compare_data.Clear();
			compare_data.Read_ClockFile(compare_file);
		}
		if (p->m_time.m_Hour >= 12 && clk12 == 0)
		{
			clk12 = 1;
			int len = strlen(compare_file);
			compare_file[len - 6] = '1';
			compare_file[len - 5] = '2';
			compare_data.Clear();
			compare_data.Read_ClockFile(compare_file);
		}
		if (p->m_time.m_Hour >= 18 && clk18 == 0)
		{
			clk18 = 1;
			int len = strlen(compare_file);
			compare_file[len - 6] = '1';
			compare_file[len - 5] = '8';
			compare_data.Clear();
			compare_data.Read_ClockFile(compare_file);
		}
		double re_clk = Clock_Interpolate_w(prn, p->m_time, &compare_data, n);
		double base_clk = p->m_clockcor;
		double clk_diff = re_clk - base_clk;

		if (base_clk>1000 || re_clk>1000)
		{
			p = p->m_next; continue;
			clk_diff = 20;
		}
		else
		{
			clk_diff *= 1.0e9;
		}

		if (base_clk<1000)
		{
			base_clk *= 1.0e9;
		}
		if (re_clk<1000)
		{
			re_clk *= 1.0e9;
		}

		fprintf(fp, "%4d%15.3f%15.3f%15.3f\n", prn, base_clk, re_clk, clk_diff);
		p = p->m_next;
	}

	fclose(fp);

	return 1;
}