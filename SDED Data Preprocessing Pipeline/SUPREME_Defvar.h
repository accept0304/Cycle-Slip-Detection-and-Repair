/* -------------------------------------------------------------------------
* SUPREME_Defvar.h : Variance defined in Multi-GNSS Solutions
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2018.12.07
* ------------------------------------------------------------------------- */
#ifndef SUPREME_DEFINE_VARIANCE_H_HH
#define SUPREME_DEFINE_VARIANCE_H_HH

#if defined(_WIN32)
#define SYS_WIND           // windows system
#define SYS_NTRIP_WIND     // ntrip in windows system
#else
#define SYS_LINUX          // linux system
#define SYS_NTRIP_LINX     // ntrip in linux system
#endif

#if defined(_MSC_VER)      // Microsoft C compiler Vertion
#if(_MSC_VER>1310)
#define PLATF_WVS
#else
#define PLATF_NVS
#endif
#endif

#ifdef SYS_WIND
#include <Windows.h>
#include <direct.h>
#include <io.h>
#include <float.h>
#else
#include <pthread.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/types.h>
#include <stdint.h>
#include <netinet/in.h>
#endif

#if defined(SYS_WIND)
typedef __int32 int32_t;
typedef __int64 int64_t;

typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;

typedef LPTHREAD_START_ROUTINE THREAD_ADDR;
typedef SOCKET msocket;
#endif

#ifdef SYS_NTRIP_LINX
typedef unsigned long DWORD;
typedef void*   HANDLE;
typedef void*   THREAD_ADDR;
typedef void*   LPVOID;
typedef int     msocket;
typedef pthread_mutex_t CRITICAL_SECTION;
#endif

/// Satellite Number (SatNo) and PRN defination
#define MIN_GPS_SATNO     1                                              // GPS MIN Satno
#define MAX_GPS_SATNO   100                                              // GPS MAX Satno
#define MIN_GPS_PRN       1                                              // Min GPS sat prn
#define MAX_GPS_RPN      32                                              // Max GPS sat prn
#define GPS_SATNUM     (MAX_GPS_RPN-MIN_GPS_PRN+1)                       // Gps satellite counts

#define MIN_GLO_SATNO   101                                              // GLO MIN Satno
#define MAX_GLO_SATNO   200                                              // GLO MAX Satno
#define MIN_GLO_PRN       1                                              // Min GLO sat prn
#define MAX_GLO_RPN      27                                              // Max GLO sat prn
#define GLO_SATNUM     (MAX_GLO_RPN-MIN_GLO_PRN+1)                       // Glo satellite counts

#define MIN_GAL_SATNO   201                                              // GAL MIN Satno
#define MAX_GAL_SATNO   300                                              // GAL MAX Satno
#define MIN_GAL_PRN       1                                              // Min GAL sat prn
#define MAX_GAL_RPN      36                                              // Max GAL sat prn
#define GAL_SATNUM     (MAX_GAL_RPN-MIN_GAL_PRN+1)                       // Gal satellite counts

#define MIN_BDS_SATNO   301                                              // BDS MIN Satno
#define MAX_BDS_SATNO   400                                              // BDS MAX Satno
#define MIN_BDS_PRN       1                                              // Min BDS sat prn
#define MAX_BDS_RPN      60                                              // Max BDS sat prn
#define BDS_SATNUM     (MAX_BDS_RPN-MIN_BDS_PRN+1)                       // Bds satellite counts

#define GNSS_SATNO_NUM    MAX_BDS_SATNO
#define GNSS_SAT_PRN_NUM  GPS_SATNUM+GLO_SATNUM+GAL_SATNUM+BDS_SATNUM    // GNSS satellite number

#define MAX_GNSS_PRN_NUM        1000                                     // Max GNSS Prn Number

#define ITERNMAX                8                                        // Max itern number  

#define MAXSTATION              300                                      // Max station number in batch process mode

#define SUPREME_VERSION     2.0                                              // SUPREME program version

#define repairstanum	50												//最大卫星数
#define sinulstanum		3000										    //最大模拟周跳数

#endif
