/* -------------------------------------------------------------------------
* SUPREME_ObsRNX.h : Observation Rinex
*
*           Copyright (C) by C. Zhao 2016-2018 All rights reserved
*           Contact zhaoC.@whigg.ac.cn
*
* Create : 2019.06.05
* ------------------------------------------------------------------------- */
#ifndef SUPREME_OBSERVATION_RINEX_H_HH
#define SUPREME_OBSERVATION_RINEX_H_HH

#include "SUPREME_GnssTime.h"
#include "SUPREME_Defvar.h"

/// Observation Type (Rinex2 and Rinex3)
/// Rinex 2 Obs Type
#define OBS_L1   0
#define OBS_L2   1
#define OBS_L3   2
#define OBS_L5   3
#define OBS_L6   4
#define OBS_L7   5
#define OBS_L8   6
#define OBS_C1   7
#define OBS_C2   8
#define OBS_C3   9
#define OBS_C5   10
#define OBS_C6   11
#define OBS_C7   12
#define OBS_C8   13
#define OBS_P1   14
#define OBS_P2   15
#define OBS_P3   16
#define OBS_P5   17
#define OBS_P6   18
#define OBS_P7   19
#define OBS_P8   20
#define OBS_D1   21
#define OBS_D2   22
#define OBS_D3   23
#define OBS_D5   24
#define OBS_D6   25
#define OBS_D7   26
#define OBS_D8   27
#define OBS_S1   28
#define OBS_S2   29
#define OBS_S3   30
#define OBS_S5   31
#define OBS_S6   32
#define OBS_S7   33
#define OBS_S8   34
/// Rinex 3 Obs Type
#define OBS_C1C  35
#define OBS_C1P  36
#define OBS_C1Y  37
#define OBS_C1W  38
#define OBS_C1X  39
#define OBS_C1I  40
#define OBS_C1Q  41

#define OBS_C1D  42  //COMPASS 

#define OBS_L1C  43
#define OBS_L1P  44
#define OBS_L1I  45
#define OBS_L1W  46
#define OBS_L1X  47

#define OBS_L1D  48  //COMPASS 

#define OBS_D1C  49
#define OBS_D1X  50
#define OBS_S1C  51
#define OBS_D1P  52
#define OBS_S1P  53
#define OBS_S1X  54

#define OBS_C2C  55
#define OBS_C2P  56
#define OBS_C2Y  57
#define OBS_C2W  58
#define OBS_C2X  59
#define OBS_C2I  60
#define OBS_C2Q  61
#define OBS_L2C  62
#define OBS_L2P  63
#define OBS_L2I  64
#define OBS_L2W  65
#define OBS_L2X  66
#define OBS_D2C  67
#define OBS_D2I  68
#define OBS_D2X  69
//#define OBS_D2X  66
//#define OBS_D2I  67  //S2C
#define OBS_D2P  70
#define OBS_S2P  71
#define OBS_S2X  72

#define OBS_C3C  73
#define OBS_C3P  74
#define OBS_C3Y  75
#define OBS_C3W  76
#define OBS_C3X  77
#define OBS_C3I  78
#define OBS_C3Q  79
#define OBS_L3C  80
#define OBS_L3P  81
#define OBS_L3I  82
#define OBS_L3W  83
#define OBS_L3X  84
#define OBS_D3C  85
#define OBS_D3X  86
#define OBS_S3C  87
#define OBS_D3P  88
#define OBS_S3P  89
#define OBS_S3X  90

#define OBS_C5C  91
#define OBS_C5P  92
#define OBS_C5Y  93
#define OBS_C5W  94
#define OBS_C5X  95
#define OBS_C5I  96
#define OBS_C5Q  97

#define OBS_C5D  98   //COMPASS 

#define OBS_L5C  99
#define OBS_L5P  100
#define OBS_L5I  101
#define OBS_L5W  102
#define OBS_L5X  103

#define OBS_L5D  104  //COMPASS 

#define OBS_D5C  105
#define OBS_D5X  106
#define OBS_S5C  107
#define OBS_D5P  108
#define OBS_S5P  109
#define OBS_S5X  110

#define OBS_C6C  111
#define OBS_C6P  112 
#define OBS_C6Y  113
#define OBS_C6W  114
#define OBS_C6X  115
#define OBS_C6I  116
#define OBS_C6Q  117
#define OBS_L6C  118
#define OBS_L6P  119
#define OBS_L6I  120
#define OBS_L6W  121
#define OBS_L6X  122
#define OBS_D6C  123
#define OBS_D6X  124
#define OBS_S6C  125
#define OBS_D6P  126
#define OBS_S6P  127
#define OBS_S6X  128

#define OBS_C7C  129
#define OBS_C7P  213
#define OBS_C7Y  131
#define OBS_C7W  132
#define OBS_C7X  212
#define OBS_C7I  211
#define OBS_C7Q  210

  //COMPASS

#define OBS_L7C  137
#define OBS_L7P  280

#define OBS_L7W  140
#define OBS_L7X  141

#define OBS_L7D  142  //COMPASS 

#define OBS_D7C  143
#define OBS_D7X  144
#define OBS_S7C  145
#define OBS_D7P  146
#define OBS_S7P  147
#define OBS_S7X  148

#define OBS_C8C  149
#define OBS_C8P  150
#define OBS_C8Y  151
#define OBS_C8W  152
#define OBS_C8X  153
#define OBS_C8I  215
#define OBS_C8Q  214

#define OBS_C8D  216  //COMPASS 

#define OBS_L8C  157
#define OBS_L8P  158
#define OBS_L8I  159
#define OBS_L8W  160
#define OBS_L8X  161

#define OBS_L8D  162

#define OBS_D8C  163
#define OBS_D8X  164
#define OBS_S8C  165
#define OBS_D8P  156
#define OBS_S8P  157
#define OBS_S8X  168



#define OBS_L5Q  169

#define OBS_D2W  170
#define OBS_L2L  171
#define OBS_C2L  172
#define OBS_C6A  173
#define OBS_L6A  174
#define OBS_C6B  175
#define OBS_L6B  176
#define OBS_L7I  178

#define OBS_C4A  179
#define OBS_L4A  180
#define OBS_C4B  181
#define OBS_L4B  182
#define OBS_C4X  183
#define OBS_L4X  184

#define OBS_L3Q  186


#define OBS_L7Q  189

#define OBS_L8Q  191
#define OBS_L2Q  192
#define OBS_L6Q  193
#define OBS_L7Z  194
#define OBS_C7Z  195
#define OBS_L6Z  196
#define OBS_C6Z  197
#define OBS_L6D  198
#define OBS_C6D  199
#define OBS_C7D  200

#define OBS_MAX_TYPE 500

#define MAX_HEADER_BUFFER_SIZE  4096

#define Obs_FileLineLength_Rinex2     100          // File length in rinex2
#define Obs_FileLineLength_Rinex3     400          // File length in rinex3

/// a: SatelliteData
/// b: ObsDataType
/// c: C1 P1 P2 L1 L2
#define GET_OBS_DATA(a,b,c)\
	(\
	(b->m_Reserved = b->m_TypeIndex[a.m_Prn / 100 * 100 + OBS_##c]), (b->m_Reserved<0 ? 0.0 : a.m_Data[b->m_Reserved].m_Obs)\
	)

#define GET_OBS_FLAG(a,b,c)\
	(\
	(b->m_Reserved = b->m_TypeIndex[a.m_Prn / 100 * 100 + OBS_##c]), (b->m_Reserved<0 ? 0 : a.m_Data[b->m_Reserved].m_Strength)\
	)

/// Observation data header information
struct ObsHeader
{
	char thebuffer[MAX_HEADER_BUFFER_SIZE];
};

/// Observation data type
struct ObsDataType
{
	char m_RecAnt[21];                   // Receiver Antenna Type
	double m_Version;                    // Rinex Version
	unsigned int m_TypeNum;              // Observation Data Type Number
	int m_Reserved;
	int m_TypeIndex[MAX_GNSS_PRN_NUM];   // Observation Data Type Index
};

/// Satellite data for X-Type
struct SatelliteDataValue
{
	int m_LLI;                           // Lost lock indicator
	int m_Strength;                      // signal strength
	double m_Obs;                        // Observation value
};

/// Satellite data of one sat
class SatelliteData
{
public:
	SatelliteData();
	~SatelliteData();

	unsigned int m_Prn;                  // SatNo
	int m_Week;                          // Week
	double m_Tow;                        // Second of Week
	double m_Jd;                         // Julian date

	SatelliteDataValue *m_Data;          // Observation value of the Satellite

	SatelliteData *m_Previous;           // Previous epoch data of this satellite
	SatelliteData *m_Next;               // Next epoch data of this satellite

	double Obs_GetP1();
	double Obs_GetP2();
	double Obs_GetL1();
	double Obs_GetL2();
	double Obs_GetC1();
	double Obs_GetC2();

	void Display_SatValue(unsigned int typenum);

	void Clear();

private:
	int l1_index, l2_index;
	int p1_index, p2_index;
	int c1_index, c2_index;
};

class w_SatelliteData
{
public:
	w_SatelliteData();
	unsigned int m_Prn;
	unsigned int m_EpochCount;           // Epoch number of m_Prn satellite
	SatelliteData *Sat_First_Data;       // First epoch data
	SatelliteData *Sat_Last_Data;        // Last epoch data
	w_SatelliteData *m_Next;
};

/// Observation Epoch Data
class ObsEpochData
{
public:
	ObsEpochData();
	~ObsEpochData();

	unsigned int m_SatCount;             // Observed satellite number in current epoch
	unsigned int m_EventFlag;            // Event flag

	gnsstime m_EpochTime;                // Observing epoch time

	ObsDataType *m_ObsType;              // Observation type
	SatelliteData *EpochData;            // All Satellite data of this epoch (current)

	ObsEpochData *m_Previous;            // Previous epoch data
	ObsEpochData *m_Next;                // Next epoch data

	void Display_EpochData();

	void Clear();
};

/// Observation data type
class ObsDataType_Rinex
{
public:
	ObsDataType_Rinex();
	~ObsDataType_Rinex();

	char m_SysType;                     // System type
	char *m_DataTypeArray;              // 2 characters for rinex2 data. 1 for '\0' "L1","L2"... 3 characters for rinex3 data
	unsigned int m_DataTypeCount;       // Number of observation type
	unsigned int m_ObsTypeNum;          // Observation type number of original rinex file
	unsigned int *m_Index;              // Index of observation type in rinex file and its size is m_DataTypeCount.

	void Clear();
};

/// Observation File Header Information
class ObsFileHeader
{
public:
	ObsFileHeader();
	ObsFileHeader(char* filename);

	char m_FileType;                            // Type of file: O or N...
	char m_SysType;                             // System type: "G" or "C" or "R" or "M"
	char m_MarkerName[61];                      // Marker Name.     60 characters for data, 1 for '\0'.
	char m_MarkerNumber[21];                    // Marker Number
	char m_AntennaSeria[21];                    // Anntenna Seria.  20 characters for data, 1 for '\0'.
	char m_AntennaType[21];                     // Anntenna Type.   20 characters for data, 1 for '\0'.
	char m_ReceiverSeria[21];                   // Receiver Seria   20 characters for data, 1 for '\0'.
	char m_ReceiverType[21];                    // Receiver Type    20 characters for data, 1 for '\0'.
	char m_ReceiverVers[21];                    // Receiver Version 20 characters for data, 1 for '\0'.
	char m_TimeSys[4];                          // 3 characters for data, 1 for '\0'.
	double m_RinexVersion;                      // Rinex version
	double m_Interval;                          // Observation Interval (sec)

	double m_StaPos_x, m_StaPos_y, m_StaPos_z;  // Station approximate position (X,Y,Z)
	double m_AntPos_H, m_AntPos_E, m_AntPos_N;  // Anntenna Position (E,N,H)

	gnsstime m_StartEpoch;                      // Start epoch time
	gnsstime m_EndEpoch;                        // End epoch time

	ObsDataType_Rinex m_DataTypeList;           // Rinex 2.0 observation type list
	ObsDataType_Rinex m_DataTypeList_GPS;       // Rinex 3.0 GPS observation type list
	ObsDataType_Rinex m_DataTypeList_GLO;       // Rinex 3.0 GLO observation type list
	ObsDataType_Rinex m_DataTypeList_GAL;       // Rinex 3.0 GAL observation type list
	ObsDataType_Rinex m_DataTypeList_BDS;       // Rinex 3.0 BDS observation type list

	int ReadObsFileHeader(char *filename);      // Read Observation File Header Information

	int Get_ObsType_List(FILE *fp);             // Get Rinex Observation Type List

	void DisplayHeader();

	void Clear();

private:
	int Get_ObsType_List_Rinex2(FILE *fp);      // Get Rinex 2.0 observation type
	int Get_ObsType_List_Rinex3(FILE *fp);      // Get Rinex 3.0 observation type
};

/// Real time observation data
class RealTimeObsData
{
public:
	RealTimeObsData();
	~RealTimeObsData();

	int flag;                            // data update indiator
	int week;                            // week
	double tow;                          // second of week

	ObsHeader m_ObsHeader;               // Observation header information
	ObsFileHeader m_Header;              // Header information
	ObsDataType m_DataType;              // Observation data type
	ObsDataType m_DataType_GPS;
	ObsDataType m_DataType_GLO;
	ObsDataType m_DataType_GAL;
	ObsDataType m_DataType_BDS;
	ObsEpochData m_ObsData;              // Observation data
};

/// Observation Data Class for Rinex O File
class ObsData
{
public:
	ObsData();
	ObsData(char* filename);
	~ObsData();

	unsigned int m_EpochCount;                   // Number of epoch time
	unsigned int m_SatCount;                     // Number of satellites

	ObsFileHeader m_Header;                      // Observation file header information
	ObsDataType m_DataType;                      // Observation data type	
	ObsEpochData *m_FirstEpoch;                  // First epoch observation data
	ObsEpochData *m_LastEpoch;                   // Last epoch observation data

	w_SatelliteData *First_Sat;                  // First Satellite Data
	w_SatelliteData *Last_Sat;                   // Last Satellite Data

	FILE* SkipObsRinexHeader(char *obsf);        // Skip observation file header

	w_SatelliteData* Search_Satellite_by_Prn(unsigned int prn); // Search satellite data by prn

	int Add_Satellite_to_Obs(SatelliteData *ptr_sat); // Add satellite data to obs data

	int Generate_ObsTypeIndex();                   // Generate observation type index in rinex

	int ReadObsEpochData(FILE* fp);                // Read observation epoch data

	int ReadObsFileBody(char *filename);           // Read observation rinex data body

	int ReadObsRinex(char *filename);              // Read observation rinex data

	int ReadObsRinex_EpochWise(FILE* fp);          // Read observation rinex data epoch wise


	int Clear(int i);
	int Clear();

private:
	int ReadObsEpochData_Rinex2(FILE* fp);//?
	int ReadObsEpochData_Rinex3(FILE* fp);//?

	int Generate_ObsTypeIndex2();
	int Generate_ObsTypeIndex3();

	int ReadObsFileBody_Rinex2(char *filename);//?
	int ReadObsFileBody_Rinex3(char *filename);//?

	int ReadObsRinex_EpochWise_2(FILE* fp);        // Read observation rinex2 data epoch wise
	int ReadObsRinex_EpochWise_3(FILE* fp);        // Read observation rinex3 data epoch wise
};



#endif
