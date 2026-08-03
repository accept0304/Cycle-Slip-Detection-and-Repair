#include <string.h>

#include "SUPREME_Anntenna.h"
#include "SUPREME_CommonFunction.h"
#include "SUPREME_Matrix.h"
#include "SUPREME_Coordinate.h"
#include "SUPREME_Options.h"

///---------------------------------
//#ifndef _NO_NAMESPACE
//using namespace math;
//#define STD std
//#else
//#define STD
//#endif
////
//#ifndef _NO_TEMPLATE
//typedef MatrixT<double> MatrixT;
//#else
//typedef MatrixT MatrixT;
//#endif
//---------------------------------

#define SAT_REC_ANTENNA_LINE_NUM         400

PCV_Data::PCV_Data()
{
	m_ZenithNum = 0;
	m_AzimuthNum = 0;
	m_ZenithDiff = 0;
	m_AzimuthDiff = 0;
}

PCV_Data::~PCV_Data()
{
	//
}

AntennaInfo_Sat::AntennaInfo_Sat()
{
	memset(m_Ant_Type, 0, sizeof(char)* 21);
	memset(m_Sat_Type, 0, sizeof(char)* 4);
	m_Prn = 0;
	m_PCO[0] = 0; m_PCO[1] = 0; m_PCO[2] = 0;
}
AntennaInfo_Sat::~AntennaInfo_Sat()
{
	this->Clear();
}

int AntennaInfo_Sat::Clear()
{
	/*if (this != NULL)
	{
	free(&(m_PCV));
	memset((void*)(this), 0, sizeof(AntennaInfo_Sat));
	}*/

	return 1;
}

SatelliteAntenna::SatelliteAntenna()
{
	m_PrnCount = 0;

}
SatelliteAntenna::~SatelliteAntenna()
{
	this->Clear();
}

int SatelliteAntenna::Clear()
{
	if (m_AntInfo != NULL)
	{
		AntennaInfo_Sat *ptr = m_AntInfo;

		while (ptr != NULL)
		{
			AntennaInfo_Sat * p = ptr->m_Next;
			//ptr->Clear();
			free(ptr);
			ptr = p;
		}

		if (m_IndexTable != NULL)
		{
			free(m_IndexTable);
		}
	}

	memset((void*)(this), 0, sizeof(SatelliteAntenna));

	return 1;
}

ReceiverAntenna::ReceiverAntenna()
{
	memset(m_AntennaType, 0, sizeof(char)* 21);
	m_PCO1[0] = 0; m_PCO1[1] = 0; m_PCO1[2] = 0;
	m_PCO2[0] = 0; m_PCO2[1] = 0; m_PCO2[2] = 0;
}
ReceiverAntenna::~ReceiverAntenna()
{
	//
}

unsigned int sat_ana_prn(char *strline, unsigned int len)
{
	if ((strline[0]<'0' || strline[0]>'9') && (strline[0] != ' '))
	{
		return 0;
	}
	char prn[4] = "";
	GetSubStr(strline, len, prn, 4, 0, 3);

	return str_to_i(prn);
}

int sat_ana_time(char *strline, unsigned int len, double ref_jd)
{
	char flags[4] = "";
	char str_t[4] = "SLR";
	GetSubStr(strline, len, flags, 4, 6, 3);
	if (strcmp(flags, str_t) == 0)
	{
		return 0;
	}

	char year1[5] = "", month1[3] = "", day1[3] = "", hour1[3] = "", min1[3] = "", sec1[3] = "";
	char year2[5] = "", month2[3] = "", day2[3] = "", hour2[3] = "", min2[3] = "", sec2[3] = "";

	gnsstime t1, t2;

	GetSubStr(strline, len, year1, 5, 31, 4);
	GetSubStr(strline, len, month1, 3, 36, 2);
	GetSubStr(strline, len, day1, 3, 39, 2);
	GetSubStr(strline, len, hour1, 3, 42, 2);
	GetSubStr(strline, len, min1, 3, 45, 2);
	GetSubStr(strline, len, sec1, 3, 48, 2);

	t1.m_Year = str_to_i(year1);
	t1.m_Month = str_to_i(month1);
	t1.m_Day = str_to_i(day1);
	t1.m_Hour = str_to_i(hour1);
	t1.m_Min = str_to_i(min1);
	t1.m_Sec = str_to_f(sec1);
	t1._date2jd();

	GetSubStr(strline, len, year2, 5, 52, 4);
	GetSubStr(strline, len, month2, 3, 57, 2);
	GetSubStr(strline, len, day2, 3, 60, 2);
	GetSubStr(strline, len, hour2, 3, 63, 2);
	GetSubStr(strline, len, min2, 3, 66, 2);
	GetSubStr(strline, len, sec2, 3, 69, 2);

	t2.m_Year = str_to_i(year2);
	t2.m_Month = str_to_i(month2);
	t2.m_Day = str_to_i(day2);
	t2.m_Hour = str_to_i(hour2);
	t2.m_Min = str_to_i(min2);
	t2.m_Sec = str_to_f(sec2);

	if (t2.m_Year == 0)
	{
		t2.m_Year = t1.m_Year + 100;
	}

	t2._date2jd();

	if (ref_jd<t1.m_Jd || ref_jd>t2.m_Jd)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

void sat_ana_pco(char *strline, unsigned int len, AntennaInfo_Sat *ptr)
{
	char x[11] = "", y[11] = "", z[11] = "";

	GetSubStr(strline, len, ptr->m_Sat_Type, 4, 16, 3);
	GetSubStr(strline, len, ptr->m_Ant_Type, 21, 6, 20);

	GetSubStr(strline, len, x, 11, 74, 10);
	GetSubStr(strline, len, y, 11, 84, 10);
	GetSubStr(strline, len, z, 11, 94, 10);

	ptr->m_PCO[0] = str_to_f(x);
	ptr->m_PCO[1] = str_to_f(y);
	ptr->m_PCO[2] = str_to_f(z);
}

int sat_add_pco(AntennaInfo_Sat *ptr, AntennaInfo_Sat **tail, SatelliteAntenna *sat_ant_info)
{
	if (ptr == NULL)
	{
		return 0;
	}

	if (sat_ant_info->m_AntInfo == NULL)
	{
		sat_ant_info->m_AntInfo = ptr;
		*tail = ptr;
		++sat_ant_info->m_PrnCount;
	}
	else
	{
		if (ptr->m_Prn > (*tail)->m_Prn)  // sort .
		{
			(*tail)->m_Next = ptr;
			(*tail) = ptr;
			++sat_ant_info->m_PrnCount;
		}
		else // sort.
		{
			AntennaInfo_Sat *p = sat_ant_info->m_AntInfo;
			AntennaInfo_Sat *q = sat_ant_info->m_AntInfo;

			while (p != NULL)
			{
				if (ptr->m_Prn < p->m_Prn)
				{
					break;
				}
				q = p;
				p = p->m_Next;
			}
			ptr->m_Next = q->m_Next;
			q->m_Next = ptr;
			++sat_ant_info->m_PrnCount;
		}
	}

	return 1;
}

/// Read satellite antenna PCO
int Read_Satellite_Antenna_PCO(char *SatAntFile, double jd, SatelliteAntenna& sat_ant_info)
{
	FILE *fp_s = fopen(SatAntFile, "r");

	if (fp_s == NULL)
	{
		printf("Error:SATELLITE.Ixx can not be opened.\n");
		return 0;
	}

	char sat_file_flags[] = "PART 2";
	char sat_strline[SAT_REC_ANTENNA_LINE_NUM + 1] = "";
	unsigned int sat_len = SAT_REC_ANTENNA_LINE_NUM + 1;

	while (!is_eof(fp_s))
	{
		char flags[7] = "";

		GetFileLine(fp_s, sat_strline, sat_len);
		GetSubStr(sat_strline, sat_len, flags, 7, 0, 6);

		if (strcmp(flags, sat_file_flags) == 0)
		{
			break;
		}
	}

	// skip four lines;
	GetFileLine(fp_s, sat_strline, sat_len);
	GetFileLine(fp_s, sat_strline, sat_len);
	GetFileLine(fp_s, sat_strline, sat_len);
	GetFileLine(fp_s, sat_strline, sat_len);

	// read pco;
	unsigned int old_prn = 0, prn = 0;
	GetFileLine(fp_s, sat_strline, sat_len);

	prn = sat_ana_prn(sat_strline, sat_len);
	old_prn = prn;

	AntennaInfo_Sat *sat_tail = NULL;

	do
	{
		int state = 0;
		while (old_prn == prn)
		{
			if (is_eof(fp_s))
			{
				old_prn = 1000;
				break;
			}
			if (state == 0 && sat_ana_time(sat_strline, sat_len, jd) == 1)
			{
				AntennaInfo_Sat *ptr = (AntennaInfo_Sat *)(malloc(sizeof(AntennaInfo_Sat)));
				if (ptr != NULL)
				{
					memset((void*)(ptr), 0, sizeof(AntennaInfo_Sat));
					sat_ana_pco(sat_strline, sat_len, ptr);
					ptr->m_Prn = prn;
					sat_add_pco(ptr, &sat_tail, &sat_ant_info);
				}
				state = 1;
			}
			GetFileLine(fp_s, sat_strline, sat_len);
			prn = sat_ana_prn(sat_strline, sat_len);

		}
		old_prn = prn;

	} while (old_prn<300 && old_prn>0); // gps (1-49)  glonass (101-149) gallileo (201-249)

	// create index table.
	AntennaInfo_Sat *p = sat_ant_info.m_AntInfo;
	if (p != NULL)
	{
		sat_ant_info.m_IndexTable = (AntennaInfo_Sat **)(malloc(sat_ant_info.m_PrnCount*sizeof(AntennaInfo_Sat *)));
		memset((void *)(sat_ant_info.m_IndexTable), 0, sat_ant_info.m_PrnCount*sizeof(struct sat_antenna_info *));
	}

	unsigned int i = 0;
	while (p != NULL)
	{
		sat_ant_info.m_IndexTable[i] = p;
		++i;
		p = p->m_Next;
	}

	fclose(fp_s);

	return 1;
}

int sat_find_pcv_record(FILE *fp, char *strline, unsigned int len, char *ant_type)
{
	int state = 0;

	while (!is_eof(fp))
	{
		char flags[21] = "";

		GetFileLine(fp, strline, len);
		GetSubStr(strline, len, flags, 21, 0, 20);

		if (strcmp(flags, ant_type) == 0)
		{
			state = 1;
			break;
		}
	}

	return state;
}

void ana_pcv_record(char *strline, unsigned int len, PCV_Data *p)
{
	char z_diff[4] = "", a_diff[4] = "", m_z[4] = "";
	unsigned int max_z = 0;

	GetSubStr(strline, len, z_diff, 4, 64, 3);
	GetSubStr(strline, len, a_diff, 4, 69, 3);
	GetSubStr(strline, len, m_z, 4, 74, 3);

	p->m_ZenithDiff = str_to_i(z_diff);
	p->m_AzimuthDiff = str_to_i(a_diff);
	max_z = str_to_i(m_z);

	p->m_ZenithNum = max_z / p->m_ZenithDiff + 1;
	p->m_AzimuthNum = 360 / p->m_AzimuthDiff + 1;

	if (p->m_AzimuthDiff == 360)
	{
		p->m_AzimuthNum = 1;
	}
}

int free_pcv(PCV_Data *p)
{
	unsigned int i = 0;

	if (p->m_IndexCol != NULL)
	{
		free(p->m_IndexCol);
	}

	if (p->m_IndexRow != NULL)
	{
		free(p->m_IndexRow);
	}

	if (p->m_PCV1 != NULL)
	{
		for (i = 0; i<p->m_AzimuthNum; ++i)
		{
			if (p->m_PCV1[i] != NULL)
			{
				free(p->m_PCV1[i]);
			}
		}
		free(p->m_PCV1);
		p->m_PCV1 = NULL;
	}

	if (p->m_PCV2 != NULL)
	{
		for (i = 0; i<p->m_AzimuthNum; ++i)
		{
			if (p->m_PCV2[i] != NULL)
			{
				free(p->m_PCV2[i]);
			}
		}
		free(p->m_PCV2);
		p->m_PCV2 = NULL;
	}
	p->m_AzimuthNum = 0;
	p->m_ZenithNum = 0;

	return 1;
}

int allocate_pcv(PCV_Data *p)
{
	if (p->m_AzimuthNum == 0 || p->m_ZenithNum == 0)
	{
		return 0;
	}

	p->m_IndexCol = (double *)malloc(p->m_ZenithNum*sizeof(double));
	if (p->m_IndexCol == NULL)
	{
		free_pcv(p);
		return 0;
	}

	p->m_IndexRow = (double *)malloc(p->m_AzimuthNum*sizeof(double));
	if (p->m_IndexRow == NULL)
	{
		free_pcv(p);
		return 0;
	}

	p->m_PCV1 = (double **)malloc(p->m_AzimuthNum*sizeof(double *));

	if (p->m_PCV1 == NULL)
	{
		free_pcv(p);
		return 0;
	}
	memset((void*)(p->m_PCV1), 0, sizeof(p->m_AzimuthNum*sizeof(double *)));

	p->m_PCV2 = (double **)malloc(p->m_AzimuthNum*sizeof(double *));

	if (p->m_PCV2 == NULL)
	{
		free_pcv(p);
		return 0;
	}
	memset((void*)(p->m_PCV2), 0, sizeof(p->m_AzimuthNum*sizeof(double *)));

	unsigned int i = 0, j = 0;
	int state = 1;

	for (i = 0; i<p->m_AzimuthNum; ++i)
	{
		p->m_PCV1[i] = (double *)malloc(p->m_ZenithNum*sizeof(double));
		p->m_PCV2[i] = (double *)malloc(p->m_ZenithNum*sizeof(double));
		memset((void*)(p->m_PCV1[i]), 0, p->m_ZenithNum*sizeof(double));
		memset((void*)(p->m_PCV2[i]), 0, p->m_ZenithNum*sizeof(double));

		if (p->m_PCV1[i] == NULL || p->m_PCV2[i] == NULL)
		{
			state = 0;
			break;
		}
	}

	if (state == 0)
	{
		free_pcv(p);
		return 0;
	}

	return 1;
}

int sat_find_az(FILE *fp, char *strline, unsigned int len, char *flags_az)
{
	int state = 0;
	char flags[2] = "";
	char flags_c[] = "N";

	do
	{
		GetFileLine(fp, strline, len);
		GetSubStr(strline, len, flags, 2, 4, 1);

		if (strcmp(flags, flags_az) == 0)
		{
			state = 1;
			break;
		}

	} while (strcmp(flags, flags_c) || is_eof(fp));

	return state;
}

int read_pcv_data(FILE *fp, char *strline, unsigned int len, PCV_Data *p)
{
	char flags_c[] = "AN";
	char flags_l1[] = "L1";
	char flags_l2[] = "L2";
	char flags[3] = "";

	unsigned int i = 0, j = 0, k = 0;
	char index[5] = "";

	for (k = 0; k<p->m_ZenithNum; ++k)
	{
		GetSubStr(strline, len, index, 5, 9 + k * 7, 4);
		p->m_IndexCol[k] = str_to_f(index);
	}

	do
	{
		GetFileLine(fp, strline, len);
		GetSubStr(strline, len, flags, 3, 0, 2);

		if (strcmp(flags, flags_l1) == 0)
		{
			char index[4] = "";
			GetSubStr(strline, len, index, 4, 3, 3);
			p->m_IndexRow[i] = str_to_f(index);

			char data[8] = "";
			for (k = 0; k<p->m_ZenithNum; ++k)
			{
				GetSubStr(strline, len, data, 8, 6 + k * 7, 7);
				p->m_PCV1[i][k] = str_to_f(data);
			}
			++i;
		}
		if (strcmp(flags, flags_l2) == 0)
		{
			char data[8] = "";
			for (k = 0; k<p->m_ZenithNum; ++k)
			{
				GetSubStr(strline, len, data, 8, 6 + k * 7, 7);
				p->m_PCV2[j][k] = str_to_f(data);
			}
			++j;
		}

	} while (strcmp(flags, flags_c) || is_eof(fp));

	return 1;
}

/// Read satellite antenna PCV
int Read_Satellite_Antenna_PCV(char *PhaseFile, SatelliteAntenna& sat_ant_info)
{
	FILE *fp_p = fopen(PhaseFile, "r");

	if (fp_p == NULL)
	{
		printf("Error:PHAS_COD.Ixx can not be opened.\n");
		return 0;
	}

	char   str_flag[] = "ANTENNA TYPE         DUMMY"; // 26 characters.
	char   sat_strline[SAT_REC_ANTENNA_LINE_NUM + 1] = "";
	unsigned int sat_len = SAT_REC_ANTENNA_LINE_NUM + 1;

	while (!is_eof(fp_p))
	{
		char flags[27] = "";

		GetFileLine(fp_p, sat_strline, sat_len);
		GetSubStr(sat_strline, sat_len, flags, 27, 0, 26);

		if (strcmp(flags, str_flag) == 0)
		{
			break;
		}
	}// end of    "    while(!is_eof(fp_s))   "

	unsigned int nsize = ftell(fp_p);

	AntennaInfo_Sat *ptr = sat_ant_info.m_AntInfo;

	while (ptr != NULL)
	{
		if (sat_find_pcv_record(fp_p, sat_strline, sat_len, ptr->m_Ant_Type) == 1)
		{
			ana_pcv_record(sat_strline, sat_len, &(ptr->m_PCV));
			if (allocate_pcv(&(ptr->m_PCV)))
			{
				char flags[] = "A";
				if (sat_find_az(fp_p, sat_strline, sat_len, flags))
				{
					read_pcv_data(fp_p, sat_strline, sat_len, &(ptr->m_PCV));
				}
			}
		}
		fseek(fp_p, nsize, SEEK_SET);
		ptr = ptr->m_Next;
	}


	fclose(fp_p);

	return 1;
}

void rec_clear(ReceiverAntenna *rec_ant_info)
{
	free_pcv(&(rec_ant_info->m_PCV));

	memset((void*)(rec_ant_info), 0, sizeof(ReceiverAntenna));
}

void rec_get_ant_pco(FILE *fp, char *strline, unsigned int len, ReceiverAntenna *rec_ant_info)
{
	char  flags[2] = "";
	GetSubStr(strline, len, rec_ant_info->m_AntennaType, 21, 0, 20);
	do
	{
		char freq[2] = "";
		int f = 0;
		GetSubStr(strline, len, freq, 2, 38, 1);
		f = str_to_i(freq);

		if (f == 1)
		{
			char data[9] = "";
			unsigned int i = 0;
			for (i = 0; i<3; ++i)
			{
				GetSubStr(strline, len, data, 9, 41 + 8 * i, 8);
				rec_ant_info->m_PCO1[i] = str_to_f(data);
			}
		}
		if (f == 2)
		{
			char data[9] = "";
			unsigned int i = 0;
			for (i = 0; i<3; ++i)
			{
				GetSubStr(strline, len, data, 9, 41 + 8 * i, 8);
				rec_ant_info->m_PCO2[i] = str_to_f(data);
			}
		}
		GetFileLine(fp, strline, len);
		GetSubStr(strline, len, flags, 2, 0, 1);

	} while (strcmp(flags, " ") == 0);
}

/// Read receiver antenna information
int Read_Receiver_Antenna(char* PhaseFile, char* ant_type, ReceiverAntenna& rec_ant_info)
{
	FILE *fp_p = fopen(PhaseFile, "r");

	if (fp_p == NULL)
	{
		printf("Error:PHAS_COD.Ixx can not be opened.\n");
		return 0;
	}

	char   str_flag[] = "ANTENNA TYPE         DUMMY"; // 26 characters.
	char   strline[SAT_REC_ANTENNA_LINE_NUM + 1] = "";
	unsigned int len = SAT_REC_ANTENNA_LINE_NUM + 1;

	while (!is_eof(fp_p))
	{
		char flags[27] = "";

		GetFileLine(fp_p, strline, len);
		GetSubStr(strline, len, flags, 27, 0, 26);

		if (strcmp(flags, str_flag) == 0)
		{
			break;
		}
	}

	unsigned int nsize = ftell(fp_p);

	fseek(fp_p, 0, SEEK_SET); // go to beginning of the file.

	unsigned int psize = ftell(fp_p);
	int ant_find = 0;

	while (psize<nsize)
	{
		char ant[21] = "";
		GetFileLine(fp_p, strline, len);
		GetSubStr(strline, len, ant, 21, 0, 20);
		if (strcmp(ant, ant_type) == 0)
		{
			ant_find = 1;
			break;
		}
		psize = ftell(fp_p);
	}
	if (ant_find == 0)
	{
		fclose(fp_p);
		rec_clear(&rec_ant_info);
		return 0;
	}

	// get receiver antenna data.
	rec_get_ant_pco(fp_p, strline, len, &rec_ant_info);

	fseek(fp_p, nsize, SEEK_SET);

	if (sat_find_pcv_record(fp_p, strline, len, ant_type) == 1)
	{
		ana_pcv_record(strline, len, &rec_ant_info.m_PCV);
		if (allocate_pcv(&rec_ant_info.m_PCV))
		{
			char flags[] = "A";
			if (sat_find_az(fp_p, strline, len, flags))
			{
				read_pcv_data(fp_p, strline, len, &rec_ant_info.m_PCV);
			}
		}
	}

	fclose(fp_p);

	return 1;
}

/// Find Satellite Antenna Infomation by prn
AntennaInfo_Sat* Find_SatAntennaInfo_by_Prn(SatelliteAntenna& satant, unsigned int prn)
{
	if (satant.m_AntInfo == NULL)
		return NULL;

	int i = 0, j = satant.m_PrnCount - 1;
	int d = abs(j - i);

	while (d > 1)
	{
		unsigned int k = d / 2 + i;
		if (prn < satant.m_IndexTable[k]->m_Prn)
			j = k;
		else
			i = k;

		d = abs(j - i);
	}

	if (prn = satant.m_IndexTable[i]->m_Prn)
		return satant.m_IndexTable[i];
	else if (prn == satant.m_IndexTable[j]->m_Prn)
		return satant.m_IndexTable[j];
	else
		return NULL;
}

/* decode antenna parameter field --------------------------------------------*/
static int decodef(char *p, int n, double *v)
{
	int i;

	for (i = 0; i<n; i++) v[i] = 0.0;
	for (i = 0, p = strtok(p, " "); p&&i<n; p = strtok(NULL, " ")) {
		v[i++] = atof(p)*1E-3;
	}
	return i;
}

/* add antenna parameter -----------------------------------------------------*/
static void addpcv(const pcv_t *pcv, pcvs_t *pcvs)
{
	pcv_t *pcvs_pcv;

	if (pcvs->nmax <= pcvs->n) {
		pcvs->nmax += 256;
		if (!(pcvs_pcv = (pcv_t *)realloc(pcvs->pcv, sizeof(pcv_t)*pcvs->nmax))) {
			printf("*** ERROR: addpcv: memory allocation error\n");
			free(pcvs->pcv); pcvs->pcv = NULL; pcvs->n = pcvs->nmax = 0;
			return;
		}
		pcvs->pcv = pcvs_pcv;
	}
	pcvs->pcv[pcvs->n++] = *pcv;

	//pcvs->pcv[pcv->sat] = *pcv;
}

int ReadAntex(const char *file, pcvs_t *pcvs)
{
	FILE *fp;
	static const pcv_t pcv0 = { 0 };
	pcv_t pcv;
	double neu[3], dd;
	int i, f, prn, freq = 0, j, id;
	//char buff[256], csys;
	char buff[350], csys; // modified by zcb

	if (!(fp = fopen(file, "r"))) {
		return 0;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		if (strlen(buff)<60 || strstr(buff + 60, "COMMENT")) continue;

		if (strstr(buff + 60, "START OF ANTENNA")) {
			pcv = pcv0;
		}
		if (strstr(buff + 60, "END OF ANTENNA")) {
			addpcv(&pcv, pcvs);
			continue;
		}

		if (strstr(buff + 60, "TYPE / SERIAL NO")) {
			strncpy(pcv.type, buff, 20); pcv.type[20] = '\0';
			strncpy(pcv.code, buff + 20, 20); pcv.code[20] = '\0';
			if (!(prn = (int)str2num(pcv.code, 1, 2))) continue;
			if (!strncmp(pcv.code + 3, "        ", 8)) {
				pcv.sat = GetSatNo(pcv.code, MAXANT, 0);
			}
		}
		else if (strstr(buff + 60, "VALID FROM")) {
			if (!str2time(buff, 0, 43, &pcv.ts)) continue;
		}
		else if (strstr(buff + 60, "VALID UNTIL")) {
			if (!str2time(buff, 0, 43, &pcv.te)) continue;
		}
		else if (strstr(buff + 60, "DAZI")) {
			pcv.dazi = str2num(buff, 2, 6); continue;
		}
		else if (strstr(buff + 60, "ZEN1 / ZEN2 / DZEN")) {
			pcv.zen1 = str2num(buff, 2, 6);
			pcv.zen2 = str2num(buff, 8, 6);
			pcv.dzen = str2num(buff, 14, 6);
			continue;
		}
		else if (strstr(buff + 60, "START OF FREQUENCY")) {
			if (sscanf(buff + 4, "%d", &f)<1) continue;
			//for (i=0;i<NFREQ;i++) if (freqs[i]==f) break;
			//if (i<NFREQ) freq=i+1;
			if (sscanf(buff + 3, "%c", &csys)<1) continue;
			if (csys == 'G')
				freq = f;
			else if (csys == 'R')
				freq = f + NFREQ;
			else if (csys == 'C')
				freq = f + 2 * NFREQ;
			else if (csys == 'E') {
				if (f == 1) freq = f + 3 * NFREQ;
				else if (f == 5) freq = 2 + 3 * NFREQ;
				else if (f == 6) freq = 3 + 3 * NFREQ;
				else freq = 0;
			}
			else if (csys == 'J') {
				if (f<5) freq = f + 4 * NFREQ;
				else if (f == 5) freq = 3 + 4 * NFREQ;
				else freq = 0;

			}
			else freq = 0;
		}
		else if (strstr(buff + 60, "END OF FREQUENCY")) {
			freq = 0;
		}
		else if (strstr(buff + 60, "NORTH / EAST / UP")) {
			//if (freq<1||NFREQ<freq) continue;
			if (decodef(buff, 3, neu)<3) continue;
			if (freq<1) continue;
			pcv.off[freq - 1][0] = neu[pcv.sat ? 0 : 1]; /* x or e */
			pcv.off[freq - 1][1] = neu[pcv.sat ? 1 : 0]; /* y or n */
			pcv.off[freq - 1][2] = neu[2];           /* z or u */
		}
		//接收机PCV考虑随方位角的变化
		else if (strstr(buff, "NOAZI")) {
			//if (SYS_CMP==PPP_Glo.sFlag[pcv.sat-1].sys) continue;

			//if (freq<1||NFREQ<freq) continue;
			if (freq<1) continue;

			dd = (pcv.zen2 - pcv.zen1) / pcv.dzen + 1;

			if (dd != round(dd) || dd <= 1) {
				printf("*** WARNING: zen in atx file error (d!=round(d)||d<1)!\n");
				continue;
			}

			if (pcv.dazi == 0.0) {
				//pcv.var[freq-1]=new double[int(dd)];
				i = decodef(buff + 8, (int)dd, pcv.var[freq - 1]);
				if (i <= 0) {
					printf("*** ERROR: error in reading atx (i<=0)!\n");
					continue;
				}
				else if (i != (int)dd) {
					printf("*** ERROR: error in reading atx (i!=(int)dd)!\n");
					continue;
				}
			}
			else {
				id = (int)((360 - 0) / pcv.dazi) + 1;
				//pcv.var[freq-1]=new double[int(dd)*id];

				for (i = 0; i<id; i++) {
					fgets(buff, sizeof(buff), fp);
					j = decodef(buff + 8, (int)dd, &pcv.var[freq - 1][i*(int)dd]);
					if (j <= 0) {
						printf("*** ERROR: error in reading atx (j<=0)!\n");
						continue;
					}
					else if (j != (int)dd) {
						printf("*** ERROR: error in reading atx (j!=(int)dd)!\n");
						continue;
					}
				}
			}
		}
	}
	fclose(fp);

	return 1;
}

/* search antenna parameter ----------------------------------------------------
* read satellite antenna phase center position
* args   : int    sat         I   satellite number (0: receiver antenna)
*          char   *type       I   antenna type for receiver antenna
*          gtime_t time       I   time to search parameters
*          pcvs_t *pcvs       IO  antenna parameters
* return : antenna parameter (NULL: no antenna)
*-----------------------------------------------------------------------------*/
pcv_t *searchpcv(int sat, const char *type, gtime_t time, const pcvs_t *pcvs)
{
	pcv_t *pcv;
	char buff[MAXANT], *types[2], *p;
	int i, j, n = 0;

	if (sat) { /* search satellite antenna */
		for (i = 0; i<pcvs->n; i++) {
			pcv = pcvs->pcv + i;
			if (pcv->sat != sat) continue;
			if (pcv->ts.time != 0 && timediff(pcv->ts, time)>0.0) continue;
			if (pcv->te.time != 0 && timediff(pcv->te, time)<0.0) continue;
			return pcv;
		}
	}
	else {
		strcpy(buff, type);
		for (p = strtok(buff, " "); p&&n<2; p = strtok(NULL, " ")) types[n++] = p;
		if (n <= 0) return NULL;

		/* search receiver antenna with radome at first */
		for (i = 0; i<pcvs->n; i++) {
			pcv = pcvs->pcv + i;
			for (j = 0; j<n; j++) if (!strstr(pcv->type, types[j])) break;
			if (j >= n) return pcv;
		}
		/* search receiver antenna without radome */
		for (i = 0; i<pcvs->n; i++) {
			pcv = pcvs->pcv + i;
			if (strstr(pcv->type, types[0]) != pcv->type) continue;

			printf("*** WARNING: pcv without radome is used type=%s\n", type);

			return pcv;
		}
	}
	return NULL;
}


void sat_vector(double satpos[3], double sunpos[3], char* sat_type, double iz[3], double jz[3], double kz[3])
{
	double ue[3] = { 0.0 };

	double r1 = 0.0, r2 = 0.0;
	unsigned int i = 0;

	for (i = 0; i<3; ++i)
	{
		ue[i] = sunpos[i] - satpos[i];
		r1 += (ue[i] * ue[i]);
	}

	r1 = sqrt(r1);
	r2 = sqrt(satpos[0] * satpos[0] + satpos[1] * satpos[1] + satpos[2] * satpos[2]);

	for (i = 0; i<3; ++i)
	{
		ue[i] /= r1;
		kz[i] = -1 * satpos[i] / r2;
	}

	Cross3(kz, ue, jz);
	Norm3(jz);
	Cross3(jz, kz, iz);
	Norm3(iz);

	if (strcmp("IIR", sat_type) == 0)
	{
		for (i = 0; i<3; ++i)
		{
			iz[i] = -1 * iz[i];
			jz[i] = -1 * jz[i];
		}
	}
}

/* interpolate antenna phase center variation --------------------------------*/
double interpvar0(int sat, double ang, const double *var, int bsat)
{
	int i, sys, limit = 18;
	double a;

	if (bsat)
	{
		//sys = PPP_Glo.sFlag[sat - 1].sys;

		ang = ang / 5.0;

		//if (sys==SYS_GPS) limit=14;
		//else if (sys==SYS_GLO) limit=15;

		if (ang >= limit)
		{
			if (ang>limit + 0.25)
			{
				//
			}
			return var[limit];
		}
		if (ang<0)
		{
			/*sprintf(PPP_Glo.chMsg, "*** ERROR: compute satellite antenna offset: nadir < 0\n");
			outDebug(OUTWIN, OUTFIL, 0);*/
			return var[0];
		}

		i = (int)ang;

		return var[i] * (1.0 + i - ang) + var[i + 1] * (ang - i);
	}
	else
	{
		a = ang / 5.0; /* ang=0-90 */
		i = (int)a;
		if (i<0)
		{
			/*sprintf(PPP_Glo.chMsg, "*** ERROR: compute receiver antenna offset: i<0\n");
			outDebug(OUTWIN, OUTFIL, 0);*/
			return var[0];
		}
		else if (i>18)
		{
			/*sprintf(PPP_Glo.chMsg, "*** ERROR: compute receiver antenna offset: i>18\n");
			outDebug(OUTWIN, OUTFIL, 0);*/
			return var[18];
		}
		return var[i] * (1.0 - a + i) + var[i + 1] * (a - i);
	}
}

/* interpolate antenna phase center variation --------------------------------*/
double interpvar1(double azim, double zeni, const pcv_t *pcv, int f)
{
	double p, q, pcvr = 0.0;
	int izeni, iazim;
	int i = (int)((pcv->zen2 - pcv->zen1) / pcv->dzen) + 1;

	if (i != 19) {
		printf("i!=19\n");
		//getchar();
	}

	izeni = (int)((zeni - pcv->zen1) / pcv->dzen);
	iazim = (int)(azim / pcv->dazi);

	p = zeni / pcv->dzen - izeni;
	q = azim / pcv->dazi - iazim;

	if (p >= 1 || p<0 || q >= 1 || q<0)
	{
		printf("interpvar %f\t%f\n", p, q);
		//getchar();
	}

	pcvr = (1.0 - p) * (1.0 - q) * pcv->var[f][(iazim + 0)*i + (izeni + 0)]
		+ p      * (1.0 - q) * pcv->var[f][(iazim + 0)*i + (izeni + 1)]
		+ q      * (1.0 - p) * pcv->var[f][(iazim + 1)*i + (izeni + 0)]
		+ p      * q       * pcv->var[f][(iazim + 1)*i + (izeni + 1)];

	return pcvr;
}

/* satellite antenna model ------------------------------------------------------
* compute satellite antenna phase center parameters
* args   : pcv_t *pcv       I   antenna phase center parameters
*          double nadir     I   nadir angle for satellite (rad)
*          double *dant     O   range offsets for each frequency (m)
* return : none
*-----------------------------------------------------------------------------*/
void antmodel_s(char sys, int sat, const pcv_t *pcv, double nadir, double *dant)
{
	int i=0;

	//sys = PPP_Glo.sFlag[sat - 1].sys;
	for (i = 0; i<NFREQ; i++)
	{
		//在interpvar函数里对nadir也除以了5.0，这是正确的吗？
		//输出nadir*R2D的值，都在14°以内；
		//分别使用两种方案计算alrt站2010年的数据，发现除以5.0的结果高程方向偏差系统为正；
		//不除以5.0的结果高程方向无系统偏差
		//dant[i]=interpvar(sat, nadir*R2D*5.0, pcv->var[i], true); v
		if (sys == 'G')
		{
			dant[i] = interpvar0(sat, nadir*R2D*5.0, pcv->var[i], 1);
			if (i == 2)
			{
				dant[i] = interpvar0(sat, nadir*R2D*5.0, pcv->var[1], 1);
			}
		}
		else if (sys == 'R')
		{
			dant[i] = interpvar0(sat, nadir*R2D*5.0, pcv->var[i + NFREQ], 1);
			if (i == 2)
			{
				dant[i] = interpvar0(sat, nadir*R2D*5.0, pcv->var[1 + NFREQ], 1);
			}
		}
		else if (sys == 'C')
		{
			dant[i] = interpvar0(sat, nadir*R2D*5.0, pcv->var[i + 2 * NFREQ], 1);
			if (i == 2)
			{
				dant[i] = interpvar0(sat, nadir*R2D*5.0, pcv->var[1 + 2 * NFREQ], 1);
			}
		}
		else if (sys == 'E')
		{
			dant[i] = interpvar0(sat, nadir*R2D*5.0, pcv->var[i + 3 * NFREQ], 1);
			if (i == 2)
			{
				dant[i] = interpvar0(sat, nadir*R2D*5.0, pcv->var[1 + 3 * NFREQ], 1);
			}
		}
		else if (sys == 'J')
		{
			dant[i] = interpvar0(sat, nadir*R2D*5.0, pcv->var[i + 4 * NFREQ], 1);
			if (i == 2)
			{
				dant[i] = interpvar0(sat, nadir*R2D*5.0, pcv->var[1 + 4 * NFREQ], 1);
			}
		}
	}
}

void SatAntPCV(char sys, int sat, const double *rs, const double *rr, const pcv_t *pcv, double *dant)
{
	double ru[3], rz[3], eu[3], ez[3], nadir, cosa;
	int i;

	for (i = 0; i<3; i++)
	{
		ru[i] = rr[i] - rs[i];
		rz[i] = -rs[i];
	}
	if (!Norm3(ru, eu) || !Norm3(rz, ez)) return;

	cosa = Dot(eu, ez, 3);
	cosa = cosa<-1.0 ? -1.0 : (cosa>1.0 ? 1.0 : cosa);
	nadir = acos(cosa);

	antmodel_s(sys, sat, pcv+sat, nadir, dant);
}

/* satellite antenna phase center offset --------------------------------------
* compute satellite antenna phase center offset in ecef
* args   :
*          double *rsun       I   sun position
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          int    sat         I   satellite number
*          double *dant       O   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m) (iono-free LC value)
* return : none
* ---------------------------------------------------------------------------- */
void SatAntPCO(pcv_t *pcvs, const double *rsun, const double *rs, char sys, int sat, int orb_n, double *dant)
{
	pcv_t *pcv = NULL;
	double ex[3], ey[3], ez[3], es[3], r[3];
	double gamma, C1, C2, dant1, dant2;
	int i, j = 0, k = 1;

	pcv = pcvs + sat; // modifided by zcb 2020-06-11

	/* unit vectors of satellite fixed coordinates */
	for (i = 0; i < 3; i++) r[i] = -rs[i];
	if (!Norm3(r, ez)) return;
	for (i = 0; i < 3; i++) r[i] = rsun[i] - rs[i];
	if (!Norm3(r, es)) return;
	Cross3(ez, es, r);
	if (!Norm3(r, ey)) return;
	Cross3(ey, ez, ex);

	double lam1 = 0.0, lam2 = 0.0;
	lam1 = Get_WaveLength(sat, 1, orb_n);
	lam2 = Get_WaveLength(sat, 2, orb_n);

	gamma = SQR(lam2) / SQR(lam1);
	C1 = gamma / (gamma - 1.0);
	C2 = -1.0 / (gamma - 1.0);

	int index = 0;
	switch (sys)
	{
		case 'G':{index = 0; break; }
		case 'R':{index = NFREQ; break; }
		case 'C':{index = 2 * NFREQ; break; }
		case 'E':{index = 3 * NFREQ; break; }
		case 'J':{index = 4 * NFREQ; break; }
		default:
			break;
	}

	for (i = 0; i < 3; ++i)
	{
		dant1 = pcv->off[0 + index][0] * ex[i] + pcv->off[0 + index][1] * ey[i] + pcv->off[0 + index][2] * ez[i];
		dant2 = pcv->off[1 + index][0] * ex[i] + pcv->off[1 + index][1] * ey[i] + pcv->off[1 + index][2] * ez[i];
		dant[i] = C1*dant1 + C2*dant2;
	}

}

/// Satellite Antenna PCO Correction
int SatAntPCO(char sys, int sat, double satpos[3], double sunpos[3], pcv_t *pcvs, double SatPCO[3])
{
	pcv_t *pcv = NULL;
	pcv = pcvs + sat;

	double kz[3] = { 0.0 }, jz[3] = { 0.0 }, iz[3] = { 0.0 };
	// kz , jz, iz    z -axis , y -axis, x -axis of satellite coordinate system
	sat_vector(satpos, sunpos, pcv->type, iz, jz, kz);

	MatrixT mr(3, 3), pco(3, 1), re(3, 1), vmr(3, 3);

	mr(1, 1) = iz[0];
	mr(1, 2) = iz[1];
	mr(1, 3) = iz[2];

	mr(2, 1) = jz[0];
	mr(2, 2) = jz[1];
	mr(2, 3) = jz[2];

	mr(3, 1) = kz[0];
	mr(3, 2) = kz[1];
	mr(3, 3) = kz[2];

	int index = 0;
	switch (sys)
	{
		case 'G':{index = 0; break; }
		case 'R':{index = NFREQ; break; }
		case 'C':{index = 2 * NFREQ; break; }
		case 'E':{index = 3 * NFREQ; break; }
		case 'J':{index = 4 * NFREQ; break; }
		default:
			break;
	}

	pco(1, 1) = pcv->off[index][0];
	pco(2, 1) = pcv->off[index][1];
	pco(3, 1) = pcv->off[index][2];

	vmr = mr.Inv();
	re = vmr*pco;

	SatPCO[0] = re(1, 1);
	SatPCO[1] = re(2, 1);
	SatPCO[2] = re(3, 1);

	//satpos[0] += SatPCO[0];
	//satpos[1] += SatPCO[1];
	//satpos[2] += SatPCO[2];

	return 1;
}

/// Satellite Antenna PCO Correction
int SatAntPCO(double satpos[3], double sunpos[3], AntennaInfo_Sat& sat_info, double SatPCO[3])
{
	if (strlen(sat_info.m_Sat_Type) == 0)
		return 0;

	double kz[3] = { 0.0 }, jz[3] = { 0.0 }, iz[3] = { 0.0 };
	// kz , jz, iz    z -axis , y -axis, x -axis of satellite coordinate system
	sat_vector(satpos, sunpos, sat_info.m_Sat_Type, iz, jz, kz);

	MatrixT mr(3, 3), pco(3, 1), re(3, 1), vmr(3, 3);

	mr(1, 1) = iz[0];
	mr(1, 2) = iz[1];
	mr(1, 3) = iz[2];

	mr(2, 1) = jz[0];
	mr(2, 2) = jz[1];
	mr(2, 3) = jz[2];

	mr(3, 1) = kz[0];
	mr(3, 2) = kz[1];
	mr(3, 3) = kz[2];

	pco(1, 1) = sat_info.m_PCO[0];
	pco(2, 1) = sat_info.m_PCO[1];
	pco(3, 1) = sat_info.m_PCO[2];

	vmr = mr.Inv();
	re = vmr*pco;

	SatPCO[0] = re(1, 1);
	SatPCO[1] = re(2, 1);
	SatPCO[2] = re(3, 1);

	satpos[0] += SatPCO[0];
	satpos[1] += SatPCO[1];
	satpos[2] += SatPCO[2];

	return 1;
}

/// BDS satellite antenna pco
int SatAntPCO_BDS(unsigned int prn, double satpos[3], double sunpos[3], double SatPCO[3])
{
	double kz[3] = { 0.0 }, jz[3] = { 0.0 }, iz[3] = { 0.0 };
	// kz , jz, iz    z -axis , y -axis, x -axis of satellite coordinate system
	double ue[3] = { 0.0 };

	double r1 = 0.0, r2 = 0.0;
	unsigned int i = 0;

	for (i = 0; i<3; ++i)
	{
		ue[i] = sunpos[i] - satpos[i];
		r1 += (ue[i] * ue[i]);
	}

	r1 = sqrt(r1);
	r2 = sqrt(satpos[0] * satpos[0] + satpos[1] * satpos[1] + satpos[2] * satpos[2]);

	for (i = 0; i<3; ++i)
	{
		ue[i] /= r1;
		kz[i] = -1 * satpos[i] / r2;
	}

	Cross3(kz, ue, jz);
	Norm3(jz);
	Cross3(jz, kz, iz);
	Norm3(iz);

	MatrixT mr(3, 3), pco(3, 1), re(3, 1), vmr(3, 3);

	mr(1, 1) = iz[0];
	mr(1, 2) = iz[1];
	mr(1, 3) = iz[2];

	mr(2, 1) = jz[0];
	mr(2, 2) = jz[1];
	mr(2, 3) = jz[2];

	mr(3, 1) = kz[0];
	mr(3, 2) = kz[1];
	mr(3, 3) = kz[2];

	pco(1, 1) = BDS_PCO_WUM[prn - 1][0];
	pco(2, 1) = BDS_PCO_WUM[prn - 1][1];
	pco(3, 1) = BDS_PCO_WUM[prn - 1][2];

	vmr = mr.Inv();
	re = vmr*pco;

	SatPCO[0] = re(1, 1);
	SatPCO[1] = re(2, 1);
	SatPCO[2] = re(3, 1);

	satpos[0] += SatPCO[0];
	satpos[1] += SatPCO[1];
	satpos[2] += SatPCO[2];

	return 1;
}

int sat_nadir_azimu(double station[3], double satpos[3], double sunpos[3], AntennaInfo_Sat& sat_info, double *nadir, double *azimu)
{
	double kz[3] = { 0.0 }, jz[3] = { 0.0 }, iz[3] = { 0.0 };
	// kz , jz, iz    z -axis , y -axis, x -axis of satellite coordinate system
	sat_vector(satpos, sunpos, sat_info.m_Sat_Type, iz, jz, kz);

	MatrixT mr(3, 3), dp(3, 1), re(3, 1);

	mr(1, 1) = iz[0];
	mr(1, 2) = iz[1];
	mr(1, 3) = iz[2];

	mr(2, 1) = jz[0];
	mr(2, 2) = jz[1];
	mr(2, 3) = jz[2];

	mr(3, 1) = kz[0];
	mr(3, 2) = kz[1];
	mr(3, 3) = kz[2];

	double rr = sqrt((station[0] - satpos[0])*(station[0] - satpos[0])
		+ (station[1] - satpos[1])*(station[1] - satpos[1])
		+ (station[2] - satpos[2])*(station[2] - satpos[2]));

	dp(1, 1) = (station[0] - satpos[0]) / rr;
	dp(2, 1) = (station[1] - satpos[1]) / rr;
	dp(3, 1) = (station[2] - satpos[2]) / rr;

	re = mr*dp;

	double pi = 4.0*atan(1.0);

	double ptry = re(3, 1);
	if (fabs(ptry)>1)
	{
		ptry = 1.0;
	}

	*nadir = acos(ptry)*180.0 / pi;

	double ptru = re(1, 1);
	if (fabs(ptru)<1.0e-10)
	{
		ptru += 1.0e-5;
	}
	*azimu = atan2(re(2, 1), ptru)*180.0 / pi;

	if (*azimu<0)
	{
		*azimu += 360.0;
	}

	if (fabs(*azimu)>360)
	{
		*azimu = 360;
	}
	if (fabs(*nadir)>90)
	{
		*nadir = 90;
	}

	return 1;
}

int Interp_PCV_Data(PCV_Data& pcvdata, double zenith, double azimuth, unsigned int freq, double& pcv)
{
	pcv = 0.0;

	if (zenith<0)
	{
		zenith = -1 * zenith;
	}
	if (fabs(zenith)>90)
	{
		zenith = 90;
	}

	if (azimuth<0)
	{
		azimuth += 360;
	}

	if (fabs(azimuth)>360)
	{
		azimuth = 360;
	}

	if (pcvdata.m_AzimuthNum == 0 || pcvdata.m_ZenithNum == 0 || pcvdata.m_PCV1 == NULL || pcvdata.m_PCV2 == NULL)
	{
		return 0;
	}

	int col_l = ((int)(zenith)) / ((int)(pcvdata.m_ZenithDiff));  // 0 5 10 15 20     0 1 2 3 4
	int col_r = col_l + 1;

	int row_u = ((int)(azimuth)) / ((int)(pcvdata.m_AzimuthDiff));
	int row_d = row_u + 1;

	if (col_l<0)
	{
		col_l = 0; col_r = 1;
	}

	if (row_u<0)
	{
		row_u = 0; row_d = 1;
	}

	if (col_l >= (int)(pcvdata.m_ZenithNum - 1))
	{
		col_l = pcvdata.m_ZenithNum - 1;
		col_r = col_l;
		zenith = pcvdata.m_IndexCol[col_l];
	}
	if (row_u >= (int)(pcvdata.m_AzimuthNum - 1))
	{
		row_u = pcvdata.m_AzimuthNum - 1;
		row_d = row_u;
		azimuth = pcvdata.m_IndexRow[row_u];
	}

	double temp1_1 = 0.0, temp1_2 = 0.0;
	double temp2_1 = 0.0, temp2_2 = 0.0;

	double bili = 0.0;

	if (row_d == row_u)
	{
		bili = 0.0;
	}
	else
	{
		bili = (azimuth - pcvdata.m_IndexRow[row_u]) / (pcvdata.m_IndexRow[row_d] - pcvdata.m_IndexRow[row_u]);
	}

	temp1_1 = pcvdata.m_PCV1[row_u][col_l] + bili*(pcvdata.m_PCV1[row_d][col_l] - pcvdata.m_PCV1[row_u][col_l]);
	temp1_2 = pcvdata.m_PCV1[row_u][col_r] + bili*(pcvdata.m_PCV1[row_d][col_r] - pcvdata.m_PCV1[row_u][col_r]);

	temp2_1 = pcvdata.m_PCV2[row_u][col_l] + bili*(pcvdata.m_PCV2[row_d][col_l] - pcvdata.m_PCV2[row_u][col_l]);
	temp2_2 = pcvdata.m_PCV2[row_u][col_r] + bili*(pcvdata.m_PCV2[row_d][col_r] - pcvdata.m_PCV2[row_u][col_r]);

	if (col_r == col_l)
	{
		bili = 0.0;
	}
	else
	{
		bili = (zenith - pcvdata.m_IndexCol[col_l]) / (pcvdata.m_IndexCol[col_r] - pcvdata.m_IndexCol[col_l]);
	}
	if (freq == 1)
	{
		pcv = temp1_1 + bili*(temp1_2 - temp1_1);
	}
	else
	{
		pcv = temp2_1 + bili*(temp2_2 - temp2_1);
	}

	pcv = (pcv)*0.001;        // convert (mm) to (m)

	return 1;
}

/// Satellite Antenna PCV Correction
int SatAntPCV(double sta[3], double sat[3], double sun[3], AntennaInfo_Sat& sat_info, unsigned int freq, double& pcv)
{
	pcv = 0.0;

	double nadir = 0.0, azimuth = 360.0;;

	double **pcv_p = (freq == 1) ? sat_info.m_PCV.m_PCV1 : sat_info.m_PCV.m_PCV2;

	if (&sat_info == NULL || pcv_p == NULL)
	{
		printf("Error:Satellite Antenna Information is Invalid!\n");
		return 0;
	}

	//nadir=sat_nadir(sta,sat)*180.0/4*atan(1.0); // covert to degree.
	sat_nadir_azimu(sta, sat, sun, sat_info, &nadir, &azimuth);

	Interp_PCV_Data((sat_info.m_PCV), nadir, azimuth, freq, pcv);

	return 1;
}

/// Receiver Antenna PCO Correction
int RecAntPCO(double a, double e, double sta[3], double ant[3], ReceiverAntenna& rec_info, double dv1[3], double dv2[3])
{
	Coordinate StaCoor;
	StaCoor.XYZ.setXYZ(sta[0], sta[1], sta[2]);
	StaCoor._xyz2blh();

	//double dlx[3] = { -1 * sin(StaCoor.BLH.getB())*cos(StaCoor.BLH.getL()), -1 * sin(StaCoor.BLH.getL()), cos(StaCoor.BLH.getB())*cos(StaCoor.BLH.getL()) };   // North
	//double dly[3] = { -1 * sin(StaCoor.BLH.getB())*sin(StaCoor.BLH.getL()), cos(StaCoor.BLH.getL()), cos(StaCoor.BLH.getB())*sin(StaCoor.BLH.getL()) };        // East
	//double dlz[3] = { cos(StaCoor.BLH.getB()), 0, sin(StaCoor.BLH.getB()) };                                                                                   // Up

	double dlx[3] = { sin(StaCoor.BLH._B)*cos(StaCoor.BLH._L), sin(StaCoor.BLH._L), cos(StaCoor.BLH._B)*cos(StaCoor.BLH._L) };              // North
	double dly[3] = { sin(StaCoor.BLH._B)*sin(StaCoor.BLH._L), cos(StaCoor.BLH._L), cos(StaCoor.BLH._B)*sin(StaCoor.BLH._L) };              // East
	double dlz[3] = { cos(StaCoor.BLH._B), 0, sin(StaCoor.BLH._B) };

	double pco1[3] = { rec_info.m_PCO1[0] + ant[0], rec_info.m_PCO1[1] + ant[1], rec_info.m_PCO1[2] + ant[2] };
	double pco2[3] = { rec_info.m_PCO2[0] + ant[0], rec_info.m_PCO2[1] + ant[1], rec_info.m_PCO2[2] + ant[2] };

	dv1[0] = dlx[0] * pco1[0] + dlx[1] * pco1[1] + dlx[2] * pco1[2];
	dv1[1] = dly[0] * pco1[0] + dly[1] * pco1[1] + dly[2] * pco1[2];
	dv1[2] = dlz[0] * pco1[0] + dlz[1] * pco1[1] + dlz[2] * pco1[2];

	dv2[0] = dlx[0] * pco2[0] + dlx[1] * pco2[1] + dlx[2] * pco2[2];
	dv2[1] = dly[0] * pco2[0] + dly[1] * pco2[1] + dly[2] * pco2[2];
	dv2[2] = dlz[0] * pco2[0] + dlz[1] * pco2[1] + dlz[2] * pco2[2];

	return 1;
}

/* receiver antenna model ------------------------------------------------------
* compute antenna offset by antenna phase center parameters
* args   : pcv_t *pcv       I   antenna phase center parameters
*          double *azel     I   azimuth/elevation for receiver {az,el} (rad)
*          int     opt      I   option (0:only offset,1:offset+pcv)
*          double *dant     O   range offsets for each frequency (m)
* return : none
* notes  : current version does not support azimuth dependent terms
*-----------------------------------------------------------------------------*/
void RecAntModel(char sys, const pcv_t *pcv, const double *del, const double *azel, int opt, double *dant)
{
#if 1
	double e[3], off[3], cosel = cos(azel[1]);
	int i, j, ii = 0;

	e[0] = sin(azel[0])*cosel;
	e[1] = cos(azel[0])*cosel;
	e[2] = sin(azel[1]);

	if (strlen(pcv->type) == 0)
	{
		for (i = 0; i<NFREQ; i++)
		{
			if (sys == 'G' || sys == 'C' || sys == 'E' || sys == 'J')
			{
				ii = i;
				if (i == 2) ii = 1;
			}
			else if (sys == 'R')
			{
				ii = i + NFREQ;
				if (i == 2) ii = 1 + NFREQ;
			}

			for (j = 0; j<3; j++) off[j] = pcv->off[ii][j] + del[j];

			if (Norm(pcv->off[ii], 3)>0.0)
			{
				/*sprintf(PPP_Glo.chMsg, "norm(pcv->off[ii],3)>0.0\n");
				outDebug(OUTWIN, OUTFIL, 0);*/
			}

			dant[i] = -Dot(off, e, 3) + (opt ? interpvar0(0, 90.0 - azel[1] * R2D, pcv->var[ii], 0) : 0.0);
		}

		return;
	}

	/* 考虑方位角因素 */
	for (i = 0; i<NFREQ; i++)
	{
		if (sys == 'G' || sys == 'C' || sys == 'E' || sys == 'J')
		{
			ii = i;
			if (i == 2) ii = 1;
		}
		else if (sys == 'R')
		{
			ii = i + NFREQ;
			if (i == 2) ii = 1 + NFREQ;
		}

		for (j = 0; j<3; j++) off[j] = pcv->off[ii][j] + del[j];

		if (pcv->dazi != 0.0)
			dant[i] = -Dot(off, e, 3) + interpvar1(azel[0] * R2D, 90 - azel[1] * R2D, pcv, ii);
		else
			dant[i] = -Dot(off, e, 3) + interpvar0(0, 90.0 - azel[1] * R2D, pcv->var[ii], 0);
	}
#endif
}

/// Receiver Antenna PCV Correction
int RecAntPCV(double a, double e, double sta[3], double sat[3], ReceiverAntenna& rec_info, unsigned int freq, double& pcv)
{
	pcv = 0.0;

	double **pcv_p = (freq == 1) ? rec_info.m_PCV.m_PCV1 : rec_info.m_PCV.m_PCV2;

	if (pcv_p == NULL)
	{
		return 0;
	}

	Coordinate local_pos;
	local_pos.XYZ.setXYZ(sat[0], sat[1], sat[2]);
	Cart_Crd station(sta[0], sta[1], sta[2]);

	local_pos._geo2local(a, e, station);
	local_pos.XYZ.setXYZ(local_pos.NEU._N, local_pos.NEU._E, local_pos.NEU._U);
	local_pos._xyz2pol();

	double pi = 4 * atan(1.0);
	double zenith = 90 - local_pos.Pol._E*180.0 / pi;
	double azimuth = local_pos.Pol._A*180.0 / pi;

	Interp_PCV_Data(rec_info.m_PCV, zenith, azimuth, freq, pcv);

	return 1;
}


int RecAntOffset(double a, double e, double sta[3], double ant[3], double dv[3])
{
	Coordinate StaCoor;
	StaCoor.XYZ.setXYZ(sta[0], sta[1], sta[2]);
	StaCoor._xyz2blh();

	double dlx[3] = { sin(StaCoor.BLH._B)*cos(StaCoor.BLH._L), sin(StaCoor.BLH._L), cos(StaCoor.BLH._B)*cos(StaCoor.BLH._L) };              // North
	double dly[3] = { sin(StaCoor.BLH._B)*sin(StaCoor.BLH._L), cos(StaCoor.BLH._L), cos(StaCoor.BLH._B)*sin(StaCoor.BLH._L) };              // East
	double dlz[3] = { cos(StaCoor.BLH._B), 0, sin(StaCoor.BLH._B) };

	double antdiff[3] = { ant[0], ant[1], ant[2] };

	dv[0] = dlx[0] * antdiff[0] + dlx[1] * antdiff[1] + dlx[2] * antdiff[2];
	dv[1] = dly[0] * antdiff[0] + dly[1] * antdiff[1] + dly[2] * antdiff[2];
	dv[2] = dlz[0] * antdiff[0] + dlz[1] * antdiff[1] + dlz[2] * antdiff[2];

	return 1;
}

/// Phase windup correction
int WindupCorr(double rsun[3], double rs[3], double rr[3], double& phw)
{
	double ek[3], exs[3], eys[3], ezs[3], ess[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
	double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph;
	int i;

	/* unit vector satellite to receiver */
	for (i = 0; i<3; i++) r[i] = rr[i] - rs[i];
	if (!Norm3(r, ek)) return 0;

	/* unit vectors of satellite antenna */
	for (i = 0; i<3; i++) r[i] = -rs[i];
	if (!Norm3(r, ezs)) return 0;
	for (i = 0; i<3; i++) r[i] = rsun[i] - rs[i];
	if (!Norm3(r, ess)) return 0;
	Cross3(ezs, ess, r);
	if (!Norm3(r, eys)) return 0;
	Cross3(eys, ezs, exs);

	/* unit vectors of receiver antenna */
	ecef2pos(rr, pos);
	xyz2enu(pos, E);
	exr[0] = E[1]; exr[1] = E[4]; exr[2] = E[7];    /* x = north */
	eyr[0] = -E[0]; eyr[1] = -E[3]; eyr[2] = -E[6]; /* y = west  */

	/* phase windup effect */
	Cross3(ek, eys, eks);
	Cross3(ek, eyr, ekr);
	for (i = 0; i<3; i++)
	{
		ds[i] = exs[i] - ek[i] * Dot(ek, exs, 3) - eks[i];
		dr[i] = exr[i] - ek[i] * Dot(ek, exr, 3) + ekr[i];
	}
	cosp = Dot(ds, dr, 3) / Norm(ds, 3) / Norm(dr, 3);
	if (cosp<-1.0) cosp = -1.0;
	else if (cosp> 1.0) cosp = 1.0;
	ph = acos(cosp) / 2.0 / PI;
	Cross3(ds, dr, drs);
	if (Dot(ek, drs, 3)<0.0) ph = -ph;

	phw = ph + floor(phw - ph + 0.5);      /* in cycle */

	return 1;
}
