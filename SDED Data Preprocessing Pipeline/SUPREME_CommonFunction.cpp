#include "SUPREME_CommonFunction.h"
#include "SUPREME_Options.h"
#include "SUPREME_Constant.h"

#define MAX_PATH_LEN     400        // Max line length

double Dot(const double *a, const double *b, int n)
{
	double c = 0.0;

	while (--n >= 0) c += a[n] * b[n];
	return c;
}

double Norm(const double *a, int n)
{
	return sqrt(Dot(a, a, n));
}

void Cross3(double a[3], double b[3], double des[3])
{
	des[0] = a[1] * b[2] - b[1] * a[2];
	des[1] = a[2] * b[0] - a[0] * b[2];
	des[2] = a[0] * b[1] - a[1] * b[0];
}

void Norm3(double v[3])
{
	double tmp = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

	v[0] /= tmp;
	v[1] /= tmp;
	v[2] /= tmp;
}

int Norm3(const double *a, double *b)
{
	double r;
	if ((r = Norm(a, 3)) <= 0.0) return 0;
	b[0] = a[0] / r;
	b[1] = a[1] / r;
	b[2] = a[2] / r;
	return 1;
}

/// File Exist  -----------------------------
/* Parameter
*       char  *path    I    input file path
*       return 1->exist  0->not exit
* ------------------------------------------ */
int FileExist(char *path)
{
#ifdef _WIN32
	if (_access(path, 0) == 0) return 1;
#else
	if (access(path, 0) == 0)  return 1;
#endif

	return 0;
}

/// Create file path
int CreatePath(char *path)
{
	char comd[1000] = "";

	unsigned int i = 0, len = 0;

#ifdef _WIN32
	if (_access(path, 0) == 0)  return 0;
#else
	if (access(path, 0) == 0)  return 0;
#endif

#if defined(PLATF_WVS)
	sprintf_s(comd, 999, "mkdir %s", path);
	len = strlen(comd);
	for (i = 0; i < len; ++i)
	{
		if (comd[i] == '/') comd[i] = '\\';
	}
	system(comd);
#else
	sprintf(comd, "mkdir -p %s", path);
	system(comd);
#endif

	return 1;
}

/// Is comes to end of file
int is_eof(FILE * fp)
{
	if (feof(fp))
		return 1;   // end of file
	else
		return 0;   // not end of file
}

/// string to intenger
int str_to_i(char* p)
{
	int re = 0, neg_flag = 0;  // '+' flag.

	if (p == NULL) return 0;

	while (*p != '\0')
	{
		if (*p == '+' || *p == '-')
		{
			neg_flag = (*p != '+');
		}
		else
		{
			if (*p > 47 && *p < 58)
			{
				re = re * 10 + (*p - '0');
			}
		}
		p++;
	}
	return neg_flag ? -1 * re : re;
}

/// string to double/float
double str_to_f(char* p)
{
	double re = 0, e = 10.0;

	int neg_flag = 0;  // '+' flag.

	if (p == NULL) return 0;

	while (*p != '\0')
	{
		if (*p == '+' || *p == '-')
		{
			neg_flag = (*p != '+');
		}
		else
		{
			if (*p == '.') e = 0.1;
			else
			{
				if (*p > 47 && *p<58)
				{
					if (e>1) re = re*e + (*p - '0');
					else
					{
						re += ((*p - '0')*e);
						e *= 0.1;
					}
				}
			}
		}
		p++;
	}
	return neg_flag ? -1 * re : re;
}

/* string to number ------------------------------------------------------------
* convert substring in string to number
* args   : char   *s        I   string ("... nnn.nnn ...")
*          int    i,n       I   substring position and width
* return : converted number (0.0:error)
* ----------------------------------------------------------------------------- */
double str2num(const char *s, int i, int n)
{
	double value;
	char str[256], *p = str;

	if (i < 0 || (int)strlen(s) < i || (int)sizeof(str)-1 < n) return 0.0;
	for (s += i; *s&&--n >= 0; s++) *p++ = *s == 'd' || *s == 'D' ? 'E' : *s; *p = '\0';
	return sscanf(str, "%lf", &value) == 1 ? value : 0.0;
}

int Lower2Cap(char low[5], char cap[5])
{
	for (int i = 0; i <=4; ++i)
	{
		if ((low[i] >= 'a') && (low[i] <= 'z'))
		{
			cap[i] = low[i] - 32;
		}
		else
			cap[i] = low[i];
	}

	return 1;
}

int Cap2Lower(char cap[5], char low[5])
{
	for (int i = 0; i < 4; ++i)
	{
		if ((cap[i] >= 'A') && (cap[i] <= 'Z'))
		{
			low[i] = cap[i] + 32;
		}
		else
			low[i] = cap[i];
	}

	return 1;
}

/// String copy
void strcpy_v(char* dst, int dstlen, char* src)
{
#ifdef PLATF_WVS
	strcpy_s(dst, dstlen, src);
#else
	strcpy(dst, src);
#endif
}

/// Memory copy
void memcpy_v(void* dst, int dstlen, void* src, int maxcount)
{
#ifdef PLATF_WVS
	memcpy_s(dst, dstlen, src, maxcount);
#else
	memcpy(dst, src, maxcount);
#endif
}

/// Get file string line
int GetFileLine(FILE *fp, char *strline, unsigned int len)
{
	if (feof(fp))
	{
		//comes to end of file.
		return 0;
	}
	else
	{
		// get line.
		int k = (int)(len - 1);
		int i = 0, j = 0;
		for (i = 0; i < k; ++i)
		{
			strline[i] = ' ';       // fill with space.
		}
		strline[k] = '\0';
		fgets(strline, len, fp);
		for (j = k; j >= 0; --j)
		{
			if (strline[j] == '\n')
			{
				strline[j] = ' ';   // fill with space.
				strline[j + 1] = ' ';
				break;
			}
		}
		return 1;
	}
}

/* get a specific substring in the string ------------------
* const char   *s      I        string
* int           i      I        start position of the substring
* int           n      I        length of the substring
* return : substring
* --------------------------------------------------------- */
char* GetSubStr(const char *s, int i, int n)
{
	//double value;
	char str[256], *p = str;

	for (s += i; *s&&--n >= 0; s++)
		*p++ = *s;

	*p = '\0';

	return str;
}

/* get a specific substring in the string ---------------------------------
* char           *res            I         resource string
* unsigned int   len1            I         res string length
* char           *des            O         destination string
* char           len2            I         des string length
* unsigned int   start_index     I         start position of substring
* unsigned int   count           I         substring length
* ------------------------------------------------------------------------ */
int GetSubStr(char *res, unsigned int len1, char *des, unsigned int len2,
	unsigned int start_index, unsigned int count)
{
	unsigned int i = 0, k = 0;
	for (k = 0; k<len2; ++k)
	{
		des[k] = ' ';
	}
	des[len2 - 1] = '\0';
	if (start_index>(len1 - 1) || (len1 - 1) < count || (len2 - 1) < count)
	{
		return 0;
	}
	for (i = 0; i < count; ++i)
	{
		des[i] = res[start_index + i];
	}
	return 1;
}

/* multiply MatrixT -----------------------------------------------------------*/
void matmul(const char *tr, int n, int k, int m, double alpha,
	const double *A, const double *B, double beta, double *C)
{
	double d;
	int i, j, x, f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);

	for (i = 0; i<n; i++) for (j = 0; j<k; j++) {
		d = 0.0;
		switch (f) {
		case 1: for (x = 0; x<m; x++) d += A[i + x*n] * B[x + j*m]; break;
		case 2: for (x = 0; x<m; x++) d += A[i + x*n] * B[j + x*k]; break;
		case 3: for (x = 0; x<m; x++) d += A[x + i*m] * B[x + j*m]; break;
		case 4: for (x = 0; x<m; x++) d += A[x + i*m] * B[j + x*k]; break;
		}
		if (beta == 0.0) C[i + j*n] = alpha*d; else C[i + j*n] = alpha*d + beta*C[i + j*n];
	}
}

/// Set file pointer to selection position
int SetPointerPos(FILE *fp, char *selection, char *strline, unsigned int len)
{
	char end[21] = "END OF HEADER       ";
	char strp[21] = "";

	fseek(fp, 0, SEEK_SET);

	do
	{
		if (GetFileLine(fp, strline, len) == 0)
		{
			return 0;
		}
		char tempstr[21] = "";
		GetSubStr(strline, len, tempstr, 21, 60, 20);//原来60；
		if (strcmp(tempstr, selection) == 0)
		{
			return 1;
		}
		unsigned int i = 0;
		for (; i < 20; ++i)
		{
			strp[i] = tempstr[i];
		}
	} while (strcmp(strp, end) != 0);

	return 0;
}
int SetPointerPos_(FILE* fp, char* selection, char* strline, unsigned int len,double verson_)
{
	char end[21] = "END OF HEADER       ";
	char strp[21] = "";

	fseek(fp, 0, SEEK_SET);
	int size = 0;
	if (verson_ == 3.00)size = 60;
	else size = 65;
	do
	{
		if (GetFileLine(fp, strline, len) == 0)
		{
			return 0;
		}
		char tempstr[21] = "";
		GetSubStr(strline, len, tempstr, 21, size, 20);//原来60；
		if (strcmp(tempstr, selection) == 0)
		{
			return 1;
		}
		unsigned int i = 0;
		for (; i < 20; ++i)
		{
			strp[i] = tempstr[i];
		}
	} while (strcmp(strp, end) != 0);

	return 0;
}

void sleep_v(int milisecond)
{
#ifdef SYS_NTRIP_WIND
	Sleep(milisecond);
#else
	usleep(milisecond * 1000);
#endif
}

/* string to time --------------------------------------------------------------
* convert substring in string to gtime_t struct
* args   : char   *s        I   string ("... yyyy mm dd hh mm ss ...")
*          int    i,n       I   substring position and width
*          gtime_t *t       O   gtime_t struct
* return : status (0:ok,0>:error)
* ----------------------------------------------------------------------------- */
int str2time(const char *s, int i, int n, gtime_t *t)
{
	double ep[6];
	char str[256], *p = str;

	if (i<0 || (int)strlen(s)<i || (int)sizeof(str)-1<i) return -1;
	for (s += i; *s&&--n >= 0;) *p++ = *s++; *p = '\0';
	if (sscanf(str, "%lf %lf %lf %lf %lf %lf", ep, ep + 1, ep + 2, ep + 3, ep + 4, ep + 5)<6)
		return -1;
	if (ep[0]<100.0) ep[0] += ep[0]<80.0 ? 2000.0 : 1900.0;
	*t = epoch2time(ep);
	return 0;
}

//%%% Get Median number start
int QuickSortOnce(double a[], int low, int high)
{
	// 将首元素作为枢轴。
	double pivot = a[low];
	int i = low, j = high;

	while (i < j) {
		// 从右到左，寻找首个小于pivot的元素。
		while (a[j] >= pivot && i < j) {
			j--;
		}

		// 执行到此，j已指向从右端起首个小于或等于pivot的元素。
		// 执行替换。
		a[i] = a[j];

		// 从左到右，寻找首个大于pivot的元素。
		while (a[i] <= pivot && i < j) {
			i++;
		}

		// 执行到此，i已指向从左端起首个大于或等于pivot的元素。
		// 执行替换。
		a[j] = a[i];
	}

	// 退出while循环，执行至此，必定是i=j的情况。
	// i（或j）指向的即是枢轴的位置，定位该趟排序的枢轴并将该位置返回。
	a[i] = pivot;

	return i;
}

void QuickSort(double a[], int low, int high) 
{
	if (low >= high) 
	{
		return;
	}

	int pivot = QuickSortOnce(a, low, high);

	// 对枢轴的左端进行排序。
	QuickSort(a, low, pivot - 1);

	// 对枢轴的右端进行排序。
	QuickSort(a, pivot + 1, high);
}

int GetMedian(double a[], int n) 
{
	QuickSort(a, 0, n - 1);

	if (n % 2 != 0) 
	{
		return a[n / 2];
	}
	else 
	{
		return (a[n / 2] + a[n / 2 - 1]) / 2;
	}
}
//%%% Get Median number end

unsigned int GetCountry(char *path, char sta[4], char *country, char &flag)
{
	FILE *fp = fopen(path, "r");

	if (fp)
	{
		int lenl = MAX_PATH_LEN;                 // config string length
		char strl[MAX_PATH_LEN] = { 0 };         // config string line

		while (!feof(fp))
		{
			if (GetFileLine(fp, strl, lenl) == 0) continue;

			char StaName[5] = { 0 };
			memcpy(StaName, strl, sizeof(char)* 4);

			if (strcmp(sta, StaName)==0)
			{
				memcpy(country, strl + 6, sizeof(char)* 3);
				flag = strl[10];
				
				break;
			}
		}
	}
	else
	{
		printf("Error: Open Station_Country_List.dat file failed.\n ");
		return 0;
	}

	if (fp) fclose(fp);

	return 1;
}

/// Get station from station list
unsigned int GetDownloadStationList(char *path, char sta[400][10], char &flag)
{
	unsigned int count = 0;
	int lenl = 100 + 1;                 // config string length
	char strl[100 + 1] = { 0 };         // config string line

	FILE *fp = fopen(path, "r");
	if (fp)
	{
		while (!feof(fp))
		{
			if (count >= 400) return count;

			if (GetFileLine(fp, strl, lenl) == 0) continue;
			if (strl[0] == '#' || strl[0] == ' ' || strl[0] == '!' || strl[0] == '%') continue;

			char sta_name_country[10] = { 0 };

			memcpy(sta_name_country, strl, sizeof(char)* 9);
			strcpy(sta[count], sta_name_country);
			count++;
		}
	}
	if (fp) fclose(fp);

	return count;
}

unsigned int GetStationList(char *path, char sta[400][5])
{
	unsigned int count = 0;
	int lenl = 100 + 1;                 // config string length
	char strl[100 + 1] = { 0 };         // config string line

	FILE *fp = fopen(path, "r");
	if (fp)
	{
		while (!feof(fp))
		{
			if (count >= 400) return count;

			if (GetFileLine(fp, strl, lenl) == 0) continue;
			if (strl[0] == '#' || strl[0] == ' ' || strl[0] == '!' || strl[0] == '%') continue;

			char sta_name[5] = { 0 };

			memcpy(sta_name, strl, sizeof(char)* 4);
			strcpy(sta[count], sta_name);
			count++;
		}
	}
	if (fp) fclose(fp);

	return count;
}

/// Calculate Satellite Number
unsigned int GetSatNo(char *prnstr, unsigned int len, char sys_flag)
{
	char satsys = ' ';           // record the satellite constellation mark.
	unsigned int satno = 0;      // record satno number,if glonass, PRN+100,if compass,PRN+300.
	satsys = prnstr[0];
	satno = str_to_i(prnstr + 1);
	if (satno == 0)
	{
		return 0;
	}
	switch (satsys)
	{
		case ' ':
		{
			satno = ((sys_flag == 'G' || sys_flag == 0) ? satno : 0);
			break;
		}
		case 'G':
		{
			satno = ((sys_flag == 'G' || sys_flag == 0) ? (satno + MIN_GPS_SATNO - 1) : 0);
			break;
		}
		case 'R':
		{
			satno = ((sys_flag == 'R' || sys_flag == 0) ? (satno + MIN_GLO_SATNO - 1) : 0);
			break;
		}
		case 'C':
		{
			satno = ((sys_flag == 'C' || sys_flag == 0) ? (satno + MIN_BDS_SATNO - 1) : 0);
			break;
		}
		case 'E':
		{
			satno = ((sys_flag == 'E' || sys_flag == 0) ? (satno + MIN_GAL_SATNO - 1) : 0);
			break;
		}
		default:
		{
			satno = 0;
			break;
		}
	}
	return satno;
}

int get_fcbindex(string type) 
{
	if (type == "B2")	return 1;
	if (type == "B1C")	 return 2;
	if (type == "B2a")	 return 3;
	if (type == "E5b")	 return 4;
	if (type == "E5(ab)")return 5;
	if (type == "E6")	 return 6;
	if (type == "L5")	return 7;
	return 0;
}
int get_fcbcofe(string type)
{
	if (type == "B2")	return 1;
	if (type == "B1C")	 return 1;
	if (type == "B2a")	 return 1;
	if (type == "E5b")	 return 1;
	if (type == "E5(ab)")return 1;
	if (type == "E6")	 return 1;
	if (type == "L5")	return 1;
	return 0;
}
/// Get system by satno (return system char flag)
char GetSystem(unsigned int satno)
{
	if (satno > (MIN_GPS_SATNO - 1) && satno<(MAX_GPS_SATNO))
	{
		return 'G';
	}
	else if (satno>(MIN_GLO_SATNO - 1) && satno<(MAX_GLO_SATNO))
	{
		return 'R';
	}
	else if (satno>(MIN_GAL_SATNO - 1) && satno<(MAX_GAL_SATNO))
	{
		return 'E';
	}
	else if (satno>(MIN_BDS_SATNO - 1) && satno < (MAX_BDS_SATNO))
	{
		return 'C';
	}
	else
	{
		return 'S';
	}
}
char GetSystem_BDS_2_3(unsigned int satno)
{
	
	if (satno > (MIN_BDS_SATNO - 1) && satno < (MIN_BDS_SATNO - 1+17))
	{
		return 'B';
	}
	else if((satno > (MIN_BDS_SATNO - 1 + 17) && satno < (MAX_BDS_SATNO - 1 )))
	{
		return 'C';
	}
}
char GetSystem_GREC2C3(unsigned int satno)
{
	if (satno > (MIN_GPS_SATNO - 1) && satno < (MAX_GPS_SATNO))
	{
		return 'G';
	}
	else if (satno > (MIN_GLO_SATNO - 1) && satno < (MAX_GLO_SATNO))
	{
		return 'R';
	}
	else if (satno > (MIN_GAL_SATNO - 1) && satno < (MAX_GAL_SATNO))
	{
		return 'E';
	}
	else if (satno > (MIN_BDS_SATNO - 1) && satno < (MAX_BDS_SATNO))
	{
		return GetSystem_BDS_2_3(satno);
	}
	else
	{
		return 'S';
	}
}
char GetSystem_BDS_2_3_prn(const unsigned int& satno, unsigned int& prn)
{

	if (satno > (MIN_BDS_SATNO - 1) && satno < (MIN_BDS_SATNO - 1 + 17))
	{
		prn = satno - (MIN_BDS_SATNO - 1);
		return 'B';
	}
	else if ((satno > (MIN_BDS_SATNO - 1 + 17) && satno < (MAX_BDS_SATNO - 1)))
	{
		prn = satno - (MIN_BDS_SATNO - 1 + 17);
		return 'C';
	}
}
/// Get Satellite Prn and System flag by SatNo
char GetSysPrn(const unsigned int &satno, unsigned int &prn)
{
	if (satno > (MIN_GPS_SATNO - 1) && satno<(MAX_GPS_SATNO))
	{
		prn = satno;
		return 'G';
	}
	else if (satno>(MIN_GLO_SATNO - 1) && satno<(MAX_GLO_SATNO))
	{
		prn = satno - (MIN_GLO_SATNO - 1);
		return 'R';
	}
	else if (satno>(MIN_GAL_SATNO - 1) && satno<(MAX_GAL_SATNO))
	{
		prn = satno - (MIN_GAL_SATNO - 1);
		return 'E';
	}
	else if (satno>(MIN_BDS_SATNO - 1) && satno < (MAX_BDS_SATNO))
	{
		prn = satno - (MIN_BDS_SATNO - 1);
		return 'C';
	}
	else
	{
		return 'S';
	}
}
char GetSysPrnGREC3C2(const unsigned int& satno, unsigned int& prn)
{
	if (satno > (MIN_GPS_SATNO - 1) && satno < (MAX_GPS_SATNO))
	{
		prn = satno;
		return 'G';
	}
	else if (satno > (MIN_GLO_SATNO - 1) && satno < (MAX_GLO_SATNO))
	{
		prn = satno - (MIN_GLO_SATNO - 1);
		return 'R';
	}
	else if (satno > (MIN_GAL_SATNO - 1) && satno < (MAX_GAL_SATNO))
	{
		prn = satno - (MIN_GAL_SATNO - 1);
		return 'E';
	}
	if (satno > (MIN_BDS_SATNO - 1) && satno < (MIN_BDS_SATNO - 1 + 17))
	{
		prn = satno - (MIN_BDS_SATNO - 1);
		return 'B';
	}
	else if ((satno > (MIN_BDS_SATNO - 1 + 17) && satno < (MAX_BDS_SATNO - 1)))
	{
		prn = satno - (MIN_BDS_SATNO - 1 + 17);
		return 'C';
	}
	else
	{
		return 'S';
	}
}
/// Get system number by config setting
unsigned int GetSysNum(unsigned int sys_setting)
{
	int sysn = 0;

	if (sys_setting&PRO_SYS_GPS)
		sysn++;

	if (sys_setting&PRO_SYS_GLO)
		sysn++;

	if (sys_setting&PRO_SYS_BDS)
		sysn++;

	if (sys_setting&PRO_SYS_GAL)
		sysn++;

	return sysn;
}

/// Crx to Rnx (convert d file to o file) -----------------------------------------------
/// Parameter:
///             char  *dfile       file need to be converted
///             int    flag        0: not delete dfile  1: delete dfile after converting
//---------------------------------------------------------------------------------------
int Crx2Rnx(char *dfile, int flag)
{
	char crx2rnx_command[300] = { 0 };

	if (!FileExist(dfile))
	{
		printf("Warrning: %s is not exist...\n", dfile);
		return -1;
	}

	/*if (!FileExist("bin/crx2rnx.exe"))
	{
		printf("Warrning: CRX2RNX executable file is not existing...\n");
		return 0;
	}*/

	sprintf(crx2rnx_command, ".\\bin\\crx2rnx %s", dfile);
	system(crx2rnx_command);

	if (flag) remove(dfile);

	return 1;
}

/// Get day number by start time and end time
int GetDayNumber(gnsstime &start_time, gnsstime &end_time)
{
	int DayN = 0;

	DayN = int(end_time - start_time) / 86400 + 1;

	return DayN;
}

/// Get frequency by satno and slot number
double Get_Frequency(unsigned int satno, unsigned int freq, int orbit_n)
{
	double frequency = -1.0;
	if (satno < MAX_GPS_SATNO)
	{
		switch (freq)
		{
			case 1:{frequency = GPS_FREQU_L1; break; }
			case 2:{frequency = GPS_FREQU_L2; break; }
			case 5:{frequency = GPS_FREQU_L5; break; }
			default:{frequency = -1.0; break; }
		}
	}
	else if (satno < MAX_GLO_SATNO)
	{
		switch (freq)
		{
			case 1:{frequency = GLO_FREQU_L1(orbit_n); break; }
			case 2:{frequency = GLO_FREQU_L2(orbit_n); break; }
			case 3: {frequency = GLO_FREQU_L1a; break; }
			case 4: {frequency = GLO_FREQU_L2a; break; }
			case 5: {frequency = GLO_FREQU_L3; break; }
			default:{frequency = -1.0; break; }
		}
	}
	else if (satno < MAX_GAL_SATNO)
	{
		switch (freq)
		{
			case  1:{frequency = GAL_FREQU_E1; break; }
			case  2:{frequency = GAL_FREQU_E5a; break; } // 51
			case  3:{frequency = GAL_FREQU_E5b; break; } // 52
			case  4: {frequency = GAL_FREQU_E5; break; } // 51
			case  5: {frequency = GAL_FREQU_E6; break; } // 52

			default:{frequency = -1.0; break; }
		}
	}
	else if (satno < MAX_BDS_SATNO)
	{
		switch (freq)
		{
			case 1:{frequency = BDS_FREQU_B1; break; }
			case 2:{frequency = BDS_FREQU_B2; break; }
			case 3:{frequency = BDS_FREQU_B3; break; }
			case 4: {frequency = BDS_FREQU_B1C; break; }
		    case 5: {frequency = BDS_FREQU_B2a; break; }

			default:{frequency = -1.0; break; }
		}
	}
	else
	{
		frequency = -1.0;
	}
	return frequency;
}

double Get_dif_fre_var(unsigned int satno, unsigned int freq, int orbit_n)
{
	double frequency = 1.0;
	if (satno < MAX_GPS_SATNO)
	{
		switch (freq)
		{
		case 1: {frequency = GPS_FREQU_L1/ GPS_FREQU_L1; break; }
		case 2: {frequency = GPS_FREQU_L2/ GPS_FREQU_L1; break; }
		case 5: {frequency = GPS_FREQU_L5/ GPS_FREQU_L1; break; }
		default: {frequency = 1.0; break; }
		}
	}
	else if (satno < MAX_GLO_SATNO)
	{
		switch (freq)
		{
		case 1: {frequency = GLO_FREQU_L1(orbit_n)/ GLO_FREQU_L1(orbit_n); break; }
		case 2: {frequency = GLO_FREQU_L2(orbit_n)/ GLO_FREQU_L1(orbit_n); break; }
		case 3: {frequency = GLO_FREQU_L1a/ GLO_FREQU_L1(orbit_n); break; }
		case 4: {frequency = GLO_FREQU_L2a/ GLO_FREQU_L1(orbit_n); break; }
		case 5: {frequency = GLO_FREQU_L3/ GLO_FREQU_L1(orbit_n); break; }
		default: {frequency = 1.0; break; }
		}
	}
	else if (satno < MAX_GAL_SATNO)
	{
		switch (freq)
		{
		case  1: {frequency = GAL_FREQU_E1/ GAL_FREQU_E1; break; }
		case  2: {frequency = GAL_FREQU_E5a/ GAL_FREQU_E1; break; } // 51
		case  3: {frequency = GAL_FREQU_E5b/ GAL_FREQU_E1; break; } // 52
		case  4: {frequency = GAL_FREQU_E5/ GAL_FREQU_E1; break; } // 51
		case  5: {frequency = GAL_FREQU_E6/ GAL_FREQU_E1; break; } // 52

		default: {frequency = 1.0; break; }
		}
	}
	else if (satno < MAX_BDS_SATNO)
	{
		switch (freq)
		{
		case 1: {frequency = BDS_FREQU_B1/ BDS_FREQU_B1; break; }
		case 2: {frequency = BDS_FREQU_B2/ BDS_FREQU_B1; break; }
		case 3: {frequency = BDS_FREQU_B3/ BDS_FREQU_B1; break; }
		case 4: {frequency = BDS_FREQU_B1C/ BDS_FREQU_B1; break; }
		case 5: {frequency = BDS_FREQU_B2a/ BDS_FREQU_B1; break; }

		default: {frequency = 1.0; break; }
		}
	}
	else
	{
		frequency = 1.0;
	}
	return frequency;


}

/// Get carrier phase wave length
double Get_WaveLength(unsigned int satno, unsigned int freq, int orbit_n)
{
	double frequency = Get_Frequency(satno, freq, orbit_n);
	return (frequency > 0) ? LIGHTSPEED / frequency : -1.0;
}

int MaxVectorIndex(vector<double> vec)
{
	int index = -1;

	if (vec.size() == 0) return index;

	double max_value = fabs(vec[0]);
	for (int i = 0; i < vec.size(); ++i)
	{
		if (fabs(vec[i]) >= max_value)
		{
			max_value = fabs(vec[i]);
			index = i;
		}
	}

	return index;
}

double MaxVectorValue(vector<double> vec)
{
	int index = -1;

	index = MaxVectorIndex(vec);

	if (index >= 0)
		return vec[index];
	else
		return 99999999;
}

int MaxIndex(double *vec, int n)
{
	int index = -1;

	if (!vec || n==0) return index;

	double max_value = vec[0];
	for (int i = 0; i < n; ++i)
	{
		if (vec[i] >= max_value)
		{
			max_value = vec[i];
			index = i;
		}
	}

	return index;
}

double MaxValue(double *vec, int n)
{
	int index = -1;

	index = MaxIndex(vec, n);

	if (index >= 0)
		return vec[index];
	else
		return 99999999;
}

int GetSatNum(vector<unsigned int> prn, unsigned &gpsn, unsigned &glon, unsigned &bdsn, unsigned &galn)
{
	char sys = 0;
	int i = 0, n = 0;
	n = prn.size();

	gpsn = glon = galn = bdsn = 0;

	for (i = 0; i < n; ++i)
	{
		sys = GetSystem(prn[i]);
		if (sys == 'G')
			gpsn++;
		else if (sys == 'R')
			glon++;
		else if (sys == 'C')
			bdsn++;
		else if (sys == 'E')
			galn++;
	}

	return 1;
}

/// Observation Ionosphere Free Combination Model
double Obs_IonFreeCombine(unsigned int prn, int orbit_n, double pr1, double pr2, unsigned int freqn)
{
	double f1 = Get_Frequency(prn, 1, orbit_n);
	double f2 = Get_Frequency(prn, 2, orbit_n);

	f1 = f1*f1;	f2 = f2*f2;
	double f = f1 - f2;

	if (freqn == 2)
		return (f1 / f*pr1 - f2 / f*pr2);
	else if (freqn == 1)
		return (pr1 + pr2) / 2;

	return 0;
}
