/* -------------------------------------------------------------------------
* SUPREME_STACRD.cpp : GNSS Staion and CRD
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact zhaoC.@whigg.ac.cn
*
* Update date: 2017.9.16
* Update date: 2017.12.4
* ------------------------------------------------------------------------- */

#ifndef SUPREME_STATION_COORDINATE_H_HH
#define SUPREME_STATION_COORDINATE_H_HH

#include "SUPREME_Coordinate.h"

class StationCrd
{
public:
	StationCrd();
	char m_StaName[50];                 // Station name
	char m_MarkerNumber[21];            // Marker number
	char m_RecType[21];                 // Receiver type
	char m_AntType[21];                 // Antenna type

	double m_AntOffset[3];              // North, East, Up

	Cart_Crd StaCoord;                  // Station Coordinate

	StationCrd *next;

	int GetStaCrd(const char* filename);
	int GetStaSnx(char* filename);
	int GetStaInfo(char* filename, gnsstime &t);
	int Display();

private:
	int flag;
};


#endif
