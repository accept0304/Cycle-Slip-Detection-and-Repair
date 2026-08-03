/* -------------------------------------------------------------------------------
* SUPREME_ClockRNX.cpp : Precise Clock Rinex
*
*                    Copyright (C) by C. Zhao
*                    Contact caszcb@163.com
*
* Create : 2018.03.16
* -------------------------------------------------------------------------------- */

#ifndef SUPREME_CLOCK_RINEX_H_HH
#define SUPREME_CLOCK_RINEX_H_HH

#include "SUPREME_GnssTime.h"

class Clk_File_Header
{
public:
	char m_FileType;            // type of file, "C" or others
	int m_LeapSec;              // leap second
	double m_Version;           // version of file

	int Read_Clk_FileHeader(char* filename,double &verson_);
};

class Clk_EpochData
{
public:
	double m_clockcor;
	double m_clockvel;
	Clk_EpochData *m_previous;
	Clk_EpochData *m_next;
	gnsstime m_time;

	unsigned int Clk_ana_StringLine(char* strline, unsigned int len,double verson);
};

class Clk_Data_Satellite
{
public:
	unsigned int m_prn;
	unsigned int m_epochcount;
	Clk_EpochData *m_FirstEpoch;
	Clk_EpochData *m_MidEpoch;
	Clk_EpochData *m_LastEpoch;
	Clk_EpochData **m_Index_Table;

	Clk_Data_Satellite *m_next;
};

typedef Clk_EpochData* Clk_EpochData_ptr;

class ClockData
{
public:
	ClockData();
	~ClockData();
	unsigned int m_PrnCount;
	Clk_Data_Satellite *m_clkdata;
	Clk_File_Header m_Header;

	Clk_Data_Satellite* Search_clock_Prn(unsigned int prn);
	Clk_EpochData* Search_Nearest_clock_Epoch(unsigned int prn, gnsstime& time);
	Clk_EpochData* Search_Nearest_clock_Epoch_w(unsigned int prn, gnsstime& time);
	Clk_EpochData* Get_Clock_Interp_Data(unsigned int prn, gnsstime& time, int& len, int n);
	Clk_EpochData* Search_clock_by_Epoch(unsigned int prn, gnsstime& time);

	int Clock_Add_Data(unsigned int prn, Clk_EpochData* epochdata, Clk_Data_Satellite **tail_sat);

	int Read_ClkFileBody(char* filename,double verson);
	int Read_ClockFile(char* filename);

	int Clear();
};

int Clock_Search_As(FILE* fp, char* strline, unsigned int len);
int Clock_EpochData_ptr_Advance(Clk_EpochData **ptr, int count);
int Clock_EpochData_ptr_Disdance(Clk_EpochData *ptr, Clk_EpochData *des);
double Clock_Interpolate(gnsstime& time, Clk_EpochData *ptr, int len, int n);
double Clock_Interpolate_w(unsigned int prn, gnsstime& time, ClockData *clkdata, int n);

int Clock_Compare(char* base_file, char* compare_file, char* result_file, char svn[4], int n);

#endif
