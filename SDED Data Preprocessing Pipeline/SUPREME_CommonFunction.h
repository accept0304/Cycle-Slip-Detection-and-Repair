/* -------------------------------------------------------------------------
* SUPREME_CommonFunction.h : Common function used in SUPREME program
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2018.12.09
* ------------------------------------------------------------------------- */
#ifndef SUPREME_COMMONFUNCTION_H_HH
#define SUPREME_COMMONFUNCTION_H_HH
#include <stdio.h>
#include <vector>
#include "SUPREME_Defvar.h"
#include "SUPREME_GnssTime.h"
#include <string>
#pragma warning(disable:4996)
#pragma warning(disable:4172)

using namespace std;

#ifndef PI    
#define PI           (4*atan(1.0))
#endif

#define D2R          (PI/180.0)                    // Deg to Rad
#define R2D          (180.0/PI)                    // Rad to Deg
#define SQR(x)       ((x)*(x))                     // Square value
#define AS2R         (D2R/3600.0)                  // arc sec to radian

#define ELE_THRES(a)  a*D2R                       // Satellite Elevation Threshold

#define lock_t      CRITICAL_SECTION
#define initlock(f) InitializeCriticalSection(f)
#define lock(f)     EnterCriticalSection(f)
#define unlock(f)   LeaveCriticalSection(f)
#define delelock(f) DeleteCriticalSection(f)

static int satCycleSlip[MAX_GNSS_PRN_NUM] = { 0 };

double Dot(const double *a, const double *b, int n);
double Norm(const double *a, int n);
void Cross3(double a[3], double b[3], double des[3]);
void Norm3(double v[3]);
int Norm3(const double *a, double *b);

int FileExist(char *path);
int CreatePath(char *path);
int is_eof(FILE * fp);

int str_to_i(char* p);
double str_to_f(char* p);
double str2num(const char *s, int i, int n);

int Lower2Cap(char low[5], char cap[5]);  // Lower letter to capital letter
int Cap2Lower(char cap[5], char low[5]);  // Capital letter to lower letter

void strcpy_v(char* dst, int dstlen, char* src);
void memcpy_v(void* dst, int dstlen, void* src, int maxcount);

int GetFileLine(FILE *fp, char *strline, unsigned int len);
char* GetSubStr(const char *s, int i, int n);
int GetSubStr(char *res, unsigned int len1, char *des,
	unsigned int len2, unsigned int start_index, unsigned int count);

void matmul(const char *tr, int n, int k, int m, double alpha,
	const double *A, const double *B, double beta, double *C);

int SetPointerPos(FILE *fp, char *selection, char *strline, unsigned int len);
int SetPointerPos_(FILE* fp, char* selection, char* strline, unsigned int len,double version_);
void sleep_v(int milisecond);

int str2time(const char *s, int i, int n, gtime_t *t);

int GetMedian(double a[], int n);

unsigned int GetCountry(char *path, char sta[4], char *country, char &flag);
unsigned int GetStationList(char *path, char sta[200][5]);
unsigned int GetDownloadStationList(char *path, char sta[400][10], char &flag);

unsigned int GetSatNo(char *prnstr, unsigned int len, char sys_flag);

char GetSystem(unsigned int satno);
int get_fcbindex(string type);
int get_fcbcofe(string type);
char GetSystem_BDS_2_3(unsigned int satno);
char GetSystem_GREC2C3(unsigned int satno);
char GetSystem_BDS_2_3_prn(const unsigned int& satno, unsigned int& prn);
char GetSysPrn(const unsigned int &satno, unsigned int &prn);
char GetSysPrnGREC3C2(const unsigned int& satno, unsigned int& prn);
unsigned int GetSysNum(unsigned int sys_setting);

int Crx2Rnx(char *dfile, int flag);

int GetDayNumber(gnsstime &start_time, gnsstime &end_time);

double Get_Frequency(unsigned int satno, unsigned int freq, int orbit_n);
double Get_WaveLength(unsigned int satno, unsigned int freq, int orbit_n);
double Get_dif_fre_var(unsigned int satno, unsigned int freq, int orbit_n);

int MaxVectorIndex(vector<double> vec);
double MaxVectorValue(vector<double> vec);
int MaxIndex(double *vec, int n);
double MaxValue(double *vec, int n);

int GetSatNum(vector<unsigned int> prn, unsigned &gpsn, unsigned &glon, unsigned &bdsn, unsigned &galn);

double Obs_IonFreeCombine(unsigned int prn, int orbit_n, double pr1, double pr2, unsigned int freqn);

#endif
