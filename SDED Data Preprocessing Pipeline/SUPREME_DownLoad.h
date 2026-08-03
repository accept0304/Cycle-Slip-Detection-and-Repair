/* -------------------------------------------------------------------------
* SUPREME_DownLoad.h : Download GNSS Data used in SUPREME Program
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2018.12.10
*
* Note:
*      1. Analys Center: cod, cof, emr, esa, gfz, grg, igr, igs, igu, jpl,
*                        mig, mit,ngs, sio, sir, gbm, wum, com, grm, tom
* ------------------------------------------------------------------------- */
#ifndef SUPREME_DOWNLOAD_H_HH
#define SUPREME_DOWNLOAD_H_HH
#include "SUPREME_GnssTime.h"
#pragma warning(disable:4996)

#define  MAX_DOWNLOAD_COMMAND_LENGTH     400

#define _SHORTFNAME  0           // lone file name mode
#define _LONGFNAME   1           // short file name mode

enum RINEX_TYPE{ OBS, MGEX, NAV, OBSNAV, EPH, CLK, EPHCLK, DATAALL };
enum DCB_TYPE{ P1C1DCB, P1P2DCB, P2C2DCB, P1C1_P1P2, DCBALL };

int GetLoneFileName(char *staname, gnsstime &t, char *country, int flag, char *ofile);

int Download_Nav(char* sta, char *des_path, gnsstime &t);
int Download_Obs(char* sta, char *des_path, gnsstime &t, unsigned int fnamemode);

int DownloadEph(char *des_path, char *Center, gnsstime &t);
int DownloadClk(char *des_path, char *Center, gnsstime &t);
int Download_EphClk(char *des_path, char *Center, gnsstime &t, unsigned int typeflag);

int Download_DCB(char *des_path, gnsstime &t, unsigned int typeflag);
int Download_Erp(char *des_path, gnsstime &t);

int Download_SNX(char *des_path, char *Center, gnsstime &t);
int Download_STA(char *des_path);

int Download_GIM(char *des_path, char *Center, gnsstime &t);
int Download_ZPD(char *des_path, char *StaName, gnsstime &t);

#endif
