#include <string.h>
#include <stdio.h>

#include "SUPREME_ObsRNX.h"
#include "SUPREME_CommonFunction.h"


#define CHECK_TYPE_INDEX(a,b,c)\
{\
	index = -1; \
if (m_DataType.m_Version<3)\
	{\
	n = m_Header.m_DataTypeList.m_DataTypeCount; \
for (i = 0; i<n; ++i)\
		{\
if (strcmp((char*)(#c), m_Header.m_DataTypeList.m_DataTypeArray + (i * 3)) == 0)\
			{\
			index = i; \
			break; \
			}\
		}\
if (i == n)\
		{\
		index = -1; \
		}\
	}\
	else\
	{\
	n = m_Header.m_DataTypeList_##a.m_DataTypeCount; \
for (i = 0; i<n; ++i)\
		{\
if (strcmp((char*)(#c), m_Header.m_DataTypeList_##a.m_DataTypeArray + (i * 4)) == 0)\
			{\
			index = i; \
			break; \
			}\
		}\
if (i == n)\
		{\
		index = -1; \
		}\
	}\
	m_DataType.m_TypeIndex[b + OBS_##c] = index; \
}
//a :卫星系统  b:卫星系统其实编号  c:观测值类型
/// Get Rinex 2 Observation Epoch Time
int Get_Obs_Epochtime_Rinex2(char *strline, unsigned int len, ObsEpochData *ptr)
{
	char year[3] = "", month[3] = "", day[3] = "";
	char hour[3] = "", minute[3] = "", second[12] = "";
	GetSubStr(strline, len, year, 3, 1, 2);
	int y = atoi(year);
	if (y >= 80)
	{
		y += 1900;
	}
	else
	{
		y += 2000;
	}
	ptr->m_EpochTime.m_Year = y;
	GetSubStr(strline, len, month, 3, 4, 2);
	ptr->m_EpochTime.m_Month = str_to_i(month);
	GetSubStr(strline, len, day, 3, 7, 2);
	ptr->m_EpochTime.m_Day = str_to_i(day);
	GetSubStr(strline, len, hour, 3, 10, 2);
	ptr->m_EpochTime.m_Hour = str_to_i(hour);
	GetSubStr(strline, len, minute, 3, 13, 2);
	ptr->m_EpochTime.m_Min = str_to_i(minute);
	GetSubStr(strline, len, second, 12, 15, 11);
	ptr->m_EpochTime.m_Sec = str_to_f(second);

	ptr->m_EpochTime._date2jd();
	ptr->m_EpochTime._date2gpst();

	return 1;
}

/// Get Rinex 3 Observation Epoch Time
int Get_Obs_Epochtime_Rinex3(char *strline, unsigned int len, ObsEpochData *ptr)
{
	char year[5] = "", month[3] = "", day[3] = "";
	char hour[3] = "", minute[3] = "", second[12] = "";
	GetSubStr(strline, len, year, 5, 2, 4);
	int y = atoi(year);
	ptr->m_EpochTime.m_Year = y;
	GetSubStr(strline, len, month, 3, 7, 2);
	ptr->m_EpochTime.m_Month = str_to_i(month);
	GetSubStr(strline, len, day, 3, 10, 2);
	ptr->m_EpochTime.m_Day = str_to_i(day);
	GetSubStr(strline, len, hour, 3, 13, 2);
	ptr->m_EpochTime.m_Hour = str_to_i(hour);
	GetSubStr(strline, len, minute, 3, 16, 2);
	ptr->m_EpochTime.m_Min = str_to_i(minute);
	GetSubStr(strline, len, second, 12, 19, 11);
	ptr->m_EpochTime.m_Sec = str_to_f(second);

	ptr->m_EpochTime._date2jd();
	ptr->m_EpochTime._date2gpst();

	return 1;
}



SatelliteData::SatelliteData()
{
	m_Prn = 0;
	m_Week = 0;	m_Tow = 0;
	m_Jd = 0;
	m_Data = NULL;
	m_Previous = NULL; m_Next = NULL;

	c1_index = c2_index = p1_index = p2_index = l1_index = l2_index = 0;
}

SatelliteData::~SatelliteData()
{
	Clear();
}

void SatelliteData::Display_SatValue(unsigned int typenum)
{
	unsigned int i = 0;
	for (; i < typenum; ++i)
	{
		printf("%18.3f%c", m_Data[i].m_Obs, ' ');
	}
	printf("\n");
}

void SatelliteData::Clear()
{
	m_Prn = 0;
	m_Week = 0;	m_Tow = 0; m_Jd = 0;
	m_Previous = NULL; m_Next = NULL;

	if (m_Data) free(m_Data);

	m_Data = NULL;

	c1_index = c2_index = p1_index = p2_index = l1_index = l2_index = 0;
}

/// w_Satellite constructor
w_SatelliteData::w_SatelliteData()
{
	m_Prn = 0; m_EpochCount = 0;
	Sat_First_Data = NULL; Sat_Last_Data = NULL;
	m_Next = NULL;
}

/// Observatioin epoch data constructor
ObsEpochData::ObsEpochData()
{
	m_SatCount = 0; m_EventFlag = 0;
	EpochData = NULL;
	m_Previous = NULL; m_Next = NULL;
	m_ObsType = NULL;
}

ObsEpochData::~ObsEpochData()
{
	Clear();
}

void ObsEpochData::Display_EpochData()
{
	printf("%d %d %d %d %d %f %f %d\n",
		m_EpochTime.m_Year, m_EpochTime.m_Month, m_EpochTime.m_Day,
		m_EpochTime.m_Hour, m_EpochTime.m_Min, m_EpochTime.m_Sec, m_EpochTime.m_Jd, m_SatCount);

	unsigned int i = 0, k = 0;
	for (; i < m_SatCount; ++i)
	{
		printf("%3d%c", EpochData[i].m_Prn, ' ');

		for (k = 0; k < m_ObsType->m_TypeNum; ++k)
		{
			printf("%14.3f%c", EpochData[i].m_Data[k].m_Obs, ' ');
		}
		printf("\n");
	}
}

void ObsEpochData::Clear()
{
	if (EpochData) free(EpochData);

	m_SatCount = 0; m_EventFlag = 0;
	EpochData = NULL;
	m_Previous = NULL; m_Next = NULL;
	m_ObsType = NULL;
}

ObsDataType_Rinex::ObsDataType_Rinex()
{
	m_SysType = 0;
	m_DataTypeArray = NULL;
	m_Index = NULL;
	m_DataTypeCount = m_ObsTypeNum = 0;
}

ObsDataType_Rinex::~ObsDataType_Rinex()
{
	Clear();
}

void ObsDataType_Rinex::Clear()
{
	m_SysType = 0;
	m_DataTypeCount = m_ObsTypeNum = 0;

	free(m_DataTypeArray);
	free(m_Index);

	m_DataTypeArray = NULL;
	m_Index = NULL;
}

/// Observation File Header constructor
ObsFileHeader::ObsFileHeader()
{
	m_FileType = 0;	m_SysType = 0;
	m_RinexVersion = 0.0;

	m_StaPos_x = m_StaPos_y = m_StaPos_z = 0.0;    // station approximate position
	m_AntPos_E = m_AntPos_H = m_AntPos_N = 0.0;    // anntenna position

	memset(m_MarkerName, 0, sizeof(char)* 61);
	memset(m_MarkerNumber, 0, sizeof(char)* 21);
	memset(m_AntennaSeria, 0, sizeof(char)* 21);
	memset(m_AntennaType, 0, sizeof(char)* 21);
	memset(m_ReceiverSeria, 0, sizeof(char)* 21);
	memset(m_ReceiverType, 0, sizeof(char)* 21);
	memset(m_ReceiverVers, 0, sizeof(char)* 21);
	memset(m_TimeSys, 0, sizeof(char)* 4);
}

ObsFileHeader::ObsFileHeader(char* filename)
{
	ReadObsFileHeader(filename);
}

/// Read rinex observation file headers
int ObsFileHeader::ReadObsFileHeader(char *filename)
{
	FILE *fp = NULL;
	fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("Error: Observation Rinex File (%s) Open Failed.\n", filename);
#ifdef _WIN32
		system("pause");
#endif
		return 0;
	}

	/// Get observation version and system
	unsigned int len = Obs_FileLineLength_Rinex3 + 1;
	char strline[Obs_FileLineLength_Rinex3 + 1] = "";

	char sel1[21] = "RINEX VERSION / TYPE";
	if (SetPointerPos(fp, sel1, strline, len))
	{
		char version[10] = "";
		GetSubStr(strline, len, version, 10, 0, 9);
		m_RinexVersion = str_to_f(version);
		m_FileType = strline[20];
		m_SysType = strline[40];
		if ((m_FileType != 'o') && (m_FileType != 'O'))
		{
			printf("Warning: The file %s is not observation file.\n", filename);
			fclose(fp);
			return 0;
		}
	}

	/// Get station marker name
	char sel2[21] = "MARKER NAME         ";
	if (SetPointerPos(fp, sel2, strline, len))
	{
		GetSubStr(strline, len, m_MarkerName, 61, 0, 60);
	}

	/// Get station marker number
	char sel3[21] = "MARKER NUMBER       ";
	if (SetPointerPos(fp, sel3, strline, len))
	{
		GetSubStr(strline, len, m_MarkerNumber, 21, 0, 20);
	}

	/// Get antenna type information
	char sel4[21] = "ANT # / TYPE        ";
	if (SetPointerPos(fp, sel4, strline, len))
	{
		GetSubStr(strline, len, m_AntennaSeria, 21, 0, 20);   // Receiver Antenna Seria
		GetSubStr(strline, len, m_AntennaType, 21, 20, 20);   // Receiver Antenna Type
	}

	/// Get receiver seris, type and version
	char sel5[21] = "REC # / TYPE / VERS ";
	if (SetPointerPos(fp, sel5, strline, len))
	{
		GetSubStr(strline, len, m_ReceiverSeria, 21, 0, 20);   // Receiver Antenna Seria
		GetSubStr(strline, len, m_ReceiverType, 21, 20, 20);   // Receiver Antenna Type
		GetSubStr(strline, len, m_ReceiverVers, 21, 40, 20);   // Receiver Antenna Version
	}

	/// Get approximate position (WGS84)
	char sel6[21] = "APPROX POSITION XYZ ";
	if (SetPointerPos(fp, sel6, strline, len))
	{
		char x[15] = "", y[15] = "", z[15] = "";
		GetSubStr(strline, len, x, 15, 0, 14);
		m_StaPos_x = str_to_f(x);
		GetSubStr(strline, len, y, 15, 14, 14);
		m_StaPos_y = str_to_f(y);
		GetSubStr(strline, len, z, 15, 28, 14);
		m_StaPos_z = str_to_f(z);
	}

	/// Get antenna offset
	char sel7[21] = "ANTENNA: DELTA H/E/N";
	if (SetPointerPos(fp, sel7, strline, len))
	{
		char H[15] = "", E[15] = "", N[15] = "";
		GetSubStr(strline, len, H, 15, 0, 14);
		m_AntPos_H = str_to_f(H);
		GetSubStr(strline, len, E, 15, 14, 14);
		m_AntPos_E = str_to_f(E);
		GetSubStr(strline, len, N, 15, 28, 14);
		m_AntPos_N = str_to_f(N);
	}

	/// Get Observation Type List
	if (Get_ObsType_List(fp) == 0) return 0;

	/// Get Time of First Epoch
	char sel8[21] = "TIME OF FIRST OBS   ";
	if (SetPointerPos(fp, sel8, strline, len))
	{
		char year[7] = "", month[7] = "", day[7] = "";
		char hour[7] = "", minute[7] = "", second[14] = "";
		GetSubStr(strline, len, year, 7, 0, 6);
		m_StartEpoch.m_Year = str_to_i(year);                 // Year
		GetSubStr(strline, len, month, 7, 6, 6);
		m_StartEpoch.m_Month = str_to_i(month);               // Month
		GetSubStr(strline, len, day, 7, 12, 6);
		m_StartEpoch.m_Day = str_to_i(day);                   // Day
		GetSubStr(strline, len, hour, 7, 18, 6);
		m_StartEpoch.m_Hour = str_to_i(hour);                 // Hour
		GetSubStr(strline, len, minute, 7, 24, 6);
		m_StartEpoch.m_Min = str_to_i(minute);                // Minute
		GetSubStr(strline, len, second, 14, 30, 13);
		m_StartEpoch.m_Sec = str_to_f(second);                 // Second

		m_StartEpoch._date2jd();
		m_StartEpoch._date2gpst();
		if (m_SysType == 'G')
		{
			char timesys[] = "GPS";
			GetSubStr(timesys, 4, m_TimeSys, 4, 0, 3);
		}
		else if (m_SysType == 'R')
		{
			char timesys[] = "GLO";
			GetSubStr(timesys, 4, m_TimeSys, 4, 0, 3);
		}
		else if (m_SysType == 'C')
		{
			char timesys[] = "BDT";
			GetSubStr(timesys, 4, m_TimeSys, 4, 0, 3);
		}
		else if (m_SysType == 'E')
		{
			char timesys[] = "GAL";
			GetSubStr(timesys, 4, m_TimeSys, 4, 0, 3);
		}
		else
		{
			GetSubStr(strline, len, m_TimeSys, 4, 48, 3);
		}
	}

	/// Get Time of Last Epoch
	char sel9[21] = "TIME OF LAST OBS    ";
	if (SetPointerPos(fp, sel9, strline, len))
	{
		char year[7] = "", month[7] = "", day[7] = "";
		char hour[7] = "", minute[7] = "", second[14] = "";
		GetSubStr(strline, len, year, 7, 0, 6);
		m_EndEpoch.m_Year = str_to_i(year);                            // Year
		GetSubStr(strline, len, month, 7, 6, 6);
		m_EndEpoch.m_Month = str_to_i(month);                          // Month
		GetSubStr(strline, len, day, 7, 12, 6);
		m_EndEpoch.m_Day = str_to_i(day);                              // Day
		GetSubStr(strline, len, hour, 7, 18, 6);
		m_EndEpoch.m_Hour = str_to_i(hour);                            // Hour
		GetSubStr(strline, len, minute, 7, 24, 6);
		m_EndEpoch.m_Min = str_to_i(minute);                           // Minute
		GetSubStr(strline, len, second, 14, 30, 13);
		m_EndEpoch.m_Sec = str_to_f(second);                           // Second

		m_EndEpoch._date2jd();
		m_EndEpoch._date2gpst();
	}

	char sel10[21] = "INTERVAL            ";
	if (SetPointerPos(fp, sel10, strline, len))
	{
		m_Interval = str2num(strline, 0, 10);
	}

	if (fp) fclose(fp);
	return 1;
}

int ObsFileHeader::Get_ObsType_List(FILE *fp)
{
	if (m_RinexVersion < 3.0)
		Get_ObsType_List_Rinex2(fp);
	else if (m_RinexVersion < 4.0)
		Get_ObsType_List_Rinex3(fp);
	else
		return 0;

	return 1;
}

/// Get observation type list of rinex 2.00
int ObsFileHeader::Get_ObsType_List_Rinex2(FILE *fp)
{
	m_DataTypeList.m_SysType = m_SysType;

	char sel[21] = "# / TYPES OF OBSERV ";
	unsigned int len = Obs_FileLineLength_Rinex2 + 1;
	char strline[Obs_FileLineLength_Rinex2 + 1] = "";
	SetPointerPos(fp, sel, strline, len);

	char typecount[7] = { ' ' };
	GetSubStr(strline, len, typecount, 7, 0, 6);
	unsigned int count = (unsigned int)(str_to_i(typecount));        // Number of observation type
	m_DataTypeList.m_ObsTypeNum = count;                             // Number of observation type in original file

	// Temporary list for saving observation type: L1\0P1\0L2\0P2\0C1\0C2\0
	char *temp_list = (char*)(malloc(count*sizeof(char)* 3));                         // Obs type list
	unsigned int *temp_index = (unsigned int*)(malloc(count*sizeof(unsigned int)));   // Obs type indexs

	unsigned int i = 0, j = 0, real_index = 0;
	for (; i<count; ++i)
	{
		char typeobs[3] = "";
		GetSubStr(strline, len, typeobs, 3, 10 + j * 6, 2);
		if (typeobs[0] == 'C' || typeobs[0] == 'P' || typeobs[0] == 'L')
		{
			//get_substr(strline,len,temp_list+real_index*3,3,10+j*6,2);
			temp_list[real_index * 3 + 0] = typeobs[0];
			temp_list[real_index * 3 + 1] = typeobs[1];
			temp_list[real_index * 3 + 2] = typeobs[2];
			temp_index[real_index] = i;
			++real_index;
		}
		if ((++j) == 9){ j = 0; }
		//if number of observation type greater than 9,we must read another line.
		if ((((i + 1) % 9) == 0) && ((i + 1)<count))
		{
			GetFileLine(fp, strline, len);
		}
	}

	m_DataTypeList.m_DataTypeCount = real_index;

	// Put observation type information into file header struct/class
	m_DataTypeList.m_DataTypeArray = (char*)(malloc(real_index*sizeof(char)* 3));
	m_DataTypeList.m_Index = (unsigned int*)(malloc(real_index*sizeof(unsigned int)));
	for (i = 0; i<real_index; ++i)
	{
		m_DataTypeList.m_DataTypeArray[i * 3] = temp_list[i * 3];
		m_DataTypeList.m_DataTypeArray[i * 3 + 1] = temp_list[i * 3 + 1];
		m_DataTypeList.m_DataTypeArray[i * 3 + 2] = '\0';

		m_DataTypeList.m_Index[i] = temp_index[i];
	}

	// Release memory
	if (temp_index) free(temp_index);
	if (temp_list) free(temp_list);

	return 1;
}

/// Get observation type list of rinex 3.00
int ObsFileHeader::Get_ObsType_List_Rinex3(FILE *fp)
{
	char sel[21] = "SYS / # / OBS TYPES ";
	unsigned int len = Obs_FileLineLength_Rinex3 + 1;
	char strline[Obs_FileLineLength_Rinex3 + 1] = { 0 };

	while (!strstr(strline, "END OF HEADER"))
	{
		memset(strline, 0, sizeof(char)*(Obs_FileLineLength_Rinex3 + 1));
		GetFileLine(fp, strline, len);
		if (!strstr(strline, sel)) continue;

		char sys_flag = 0;
		sys_flag = strline[0];                                           // System
		char typecount[7] = { ' ' };
		GetSubStr(strline, len, typecount, 7, 3, 3);
		unsigned int count = (unsigned int)(str_to_i(typecount));        // Number of observation type

		switch (sys_flag)
		{
			case 'G':{m_DataTypeList_GPS.m_DataTypeCount = count; m_DataTypeList_GPS.m_SysType = sys_flag; break; }
			case 'R':{m_DataTypeList_GLO.m_DataTypeCount = count; m_DataTypeList_GLO.m_SysType = sys_flag; break; }
			case 'E':{m_DataTypeList_GAL.m_DataTypeCount = count; m_DataTypeList_GAL.m_SysType = sys_flag; break; }
			case 'C':{m_DataTypeList_BDS.m_DataTypeCount = count; m_DataTypeList_BDS.m_SysType = sys_flag; break; }
			default:break;
		}

		// Temporary list for saving observation type: C1C\0L1C\0C2W\0L2L\0
		
		char*temp_list = (char*)(malloc(count*sizeof(char)* 4));                         // Obs type list
		unsigned int *temp_index = (unsigned int*)(malloc(count*sizeof(unsigned int)));   // Obs type indexs
		memset(temp_list, 0, sizeof(char)*count * 4);
		memset(temp_index, 0, sizeof(unsigned int)*count);

		unsigned int i = 0, j = 0, real_index = 0;
		for (; i < count; ++i)
		{
			char typeobs[4] = "";
			GetSubStr(strline, len, typeobs, 4, 7 + j * 4, 3);
			if (typeobs[0] == 'C' || typeobs[0] == 'L' || typeobs[0] == 'S' || typeobs[0] == 'D')
			{
				temp_list[real_index * 4 + 0] = typeobs[0];
				temp_list[real_index * 4 + 1] = typeobs[1];
				temp_list[real_index * 4 + 2] = typeobs[2];
				temp_list[real_index * 4 + 3] = typeobs[3];
				temp_index[real_index] = i;
				++real_index;
			}

			if ((++j) == 13){ j = 0; }

			if ((((i + 1) % 13) == 0) && ((i + 1) < count))
			{
				GetFileLine(fp, strline, len);
			}
		}

		switch (sys_flag)
		{
			case 'G':{m_DataTypeList_GPS.m_DataTypeCount = real_index; break; }
			case 'R':{m_DataTypeList_GLO.m_DataTypeCount = real_index; break; }
			case 'E':{m_DataTypeList_GAL.m_DataTypeCount = real_index; break; }
			case 'C':{m_DataTypeList_BDS.m_DataTypeCount = real_index; break; }
			default:break;
		}

		// Put observation type information into file header struct/class
		switch (sys_flag)
		{
			case 'G':
			{
				m_DataTypeList_GPS.m_DataTypeArray = (char*)(malloc(real_index*sizeof(char)* 4));
				m_DataTypeList_GPS.m_Index = (unsigned int*)(malloc(real_index*sizeof(unsigned int)));
				break;
			}
			case 'R':
			{
				m_DataTypeList_GLO.m_DataTypeArray = (char*)(malloc(real_index*sizeof(char)* 4));
				m_DataTypeList_GLO.m_Index = (unsigned int*)(malloc(real_index*sizeof(unsigned int)));
				break;
			}
			case 'E':
			{
				m_DataTypeList_GAL.m_DataTypeArray = (char*)(malloc(real_index*sizeof(char)* 4));
				m_DataTypeList_GAL.m_Index = (unsigned int*)(malloc(real_index*sizeof(unsigned int)));
				break;
			}
			case 'C':
			{
				m_DataTypeList_BDS.m_DataTypeArray = (char*)(malloc(real_index*sizeof(char)* 4));
				m_DataTypeList_BDS.m_Index = (unsigned int*)(malloc(real_index*sizeof(unsigned int)));
				break;
			}
			default:break;
		}

		for (i = 0; i < real_index; ++i)
		{
			switch (sys_flag)
			{
				case 'G':
				{
					m_DataTypeList_GPS.m_DataTypeArray[i * 4 + 0] = temp_list[i * 4 + 0];
					m_DataTypeList_GPS.m_DataTypeArray[i * 4 + 1] = temp_list[i * 4 + 1];
					m_DataTypeList_GPS.m_DataTypeArray[i * 4 + 2] = temp_list[i * 4 + 2];
					m_DataTypeList_GPS.m_DataTypeArray[i * 4 + 3] = '\0';

					m_DataTypeList_GPS.m_Index[i] = temp_index[i];
					break;
				}
				case 'R':
				{
					m_DataTypeList_GLO.m_DataTypeArray[i * 4 + 0] = temp_list[i * 4 + 0];
					m_DataTypeList_GLO.m_DataTypeArray[i * 4 + 1] = temp_list[i * 4 + 1];
					m_DataTypeList_GLO.m_DataTypeArray[i * 4 + 2] = temp_list[i * 4 + 2];
					m_DataTypeList_GLO.m_DataTypeArray[i * 4 + 3] = '\0';

					m_DataTypeList_GLO.m_Index[i] = temp_index[i];
					break;
				}
				case 'E':
				{
					m_DataTypeList_GAL.m_DataTypeArray[i * 4 + 0] = temp_list[i * 4 + 0];
					m_DataTypeList_GAL.m_DataTypeArray[i * 4 + 1] = temp_list[i * 4 + 1];
					m_DataTypeList_GAL.m_DataTypeArray[i * 4 + 2] = temp_list[i * 4 + 2];
					m_DataTypeList_GAL.m_DataTypeArray[i * 4 + 3] = '\0';

					m_DataTypeList_GAL.m_Index[i] = temp_index[i];
					break;
				}
				case 'C':
				{
					m_DataTypeList_BDS.m_DataTypeArray[i * 4 + 0] = temp_list[i * 4 + 0];
					m_DataTypeList_BDS.m_DataTypeArray[i * 4 + 1] = temp_list[i * 4 + 1];
					m_DataTypeList_BDS.m_DataTypeArray[i * 4 + 2] = temp_list[i * 4 + 2];
					m_DataTypeList_BDS.m_DataTypeArray[i * 4 + 3] = '\0';

					m_DataTypeList_BDS.m_Index[i] = temp_index[i];
					break;
				}
				default:break;
			}
		}

		if (temp_index) free(temp_index);
		if (temp_list) free(temp_list);
	}

	return 1;
}

/// Display observation file header information
void ObsFileHeader::DisplayHeader()
{
	printf("%16s %f\n", "Version: ", m_RinexVersion);
	printf("%16s %c\n", "File type: ", m_FileType);
	printf("%16s %c\n", "System type: ", m_SysType);
	printf("%16s %s\n", "Maker name: ", m_MarkerName);
	printf("%16s %s\n", "Antenna serial: ", m_AntennaSeria);
	printf("%16s %s\n", "Receiver type: ", m_ReceiverType);
	printf("%16s %s\n", "Antenna type: ", m_AntennaType);
	printf("%16s %f\n", "Station x: ", m_StaPos_x);
	printf("%16s %f\n", "Station y: ", m_StaPos_y);
	printf("%16s %f\n", "Station z: ", m_StaPos_z);
	printf("%16s %f\n", "Antenna H: ", m_AntPos_H);
	printf("%16s %f\n", "Antenna E: ", m_AntPos_E);
	printf("%16s %f\n", "Antenna N: ", m_AntPos_N);
	printf("%16s %s\n", "Time system: ", m_TimeSys);
	printf("%16s\n", "Type list: ");
	printf("%16s", "  ");
	unsigned int i = 0;
	for (; i<m_DataTypeList.m_DataTypeCount; ++i)
	{
		printf("%6s ", m_DataTypeList.m_DataTypeArray + i * 3);
	}
	printf("\n");
	printf("%16s %d %d %d %d %d %5.1f\n", "start epoch: ",
		m_StartEpoch.m_Year, m_StartEpoch.m_Month, m_StartEpoch.m_Day,
		m_StartEpoch.m_Hour, m_StartEpoch.m_Min, m_StartEpoch.m_Sec);
	printf("%16s %d %d %d %d %d %5.1f\n", "end   epoch: ",
		m_EndEpoch.m_Year, m_EndEpoch.m_Month, m_EndEpoch.m_Day,
		m_EndEpoch.m_Hour, m_EndEpoch.m_Min, m_EndEpoch.m_Sec);
}

/// Clear obs file head information
void ObsFileHeader::Clear()
{
	m_FileType = 0;	m_SysType = 0;
	m_RinexVersion = 0.0;

	m_StaPos_x = m_StaPos_y = m_StaPos_z = 0.0;    // station approximate position
	m_AntPos_E = m_AntPos_H = m_AntPos_N = 0.0;    // anntenna position

	memset(m_MarkerName,   0, sizeof(char)* 61);
	memset(m_MarkerNumber, 0, sizeof(char)* 21);
	memset(m_AntennaSeria, 0, sizeof(char)* 21);
	memset(m_AntennaType,  0, sizeof(char)* 21);
	memset(m_ReceiverSeria,0, sizeof(char)* 21);
	memset(m_ReceiverType, 0, sizeof(char)* 21);
	memset(m_ReceiverVers, 0, sizeof(char)* 21);
	memset(m_TimeSys, 0, sizeof(char)* 4);

	m_DataTypeList.Clear();
	m_DataTypeList_GPS.Clear();
	m_DataTypeList_BDS.Clear();
	m_DataTypeList_GLO.Clear();
	m_DataTypeList_GAL.Clear();
}

RealTimeObsData::RealTimeObsData()
{
	flag = 0;
	week = 0;
	tow = 0;

	memset((void*)(&m_ObsHeader), 0, sizeof(ObsHeader));
	memset((void*)(&m_DataType), 0, sizeof(ObsDataType));
}
RealTimeObsData::~RealTimeObsData()
{
	int satcount = 0;
	satcount = m_ObsData.m_SatCount;

	for (int k = 0; k < satcount; ++k)
	{
		if (m_ObsData.EpochData[k].m_Data)
			free(m_ObsData.EpochData[k].m_Data);
	}

	if (m_ObsData.EpochData) free(m_ObsData.EpochData);
}


ObsData::ObsData()
{
	memset((void*)this, 0, sizeof(ObsData));
}
ObsData::ObsData(char* filename)
{
	ReadObsRinex(filename);
}
ObsData::~ObsData()
{
	this->Clear();
}

/// Skip Observation file header
FILE*  ObsData::SkipObsRinexHeader(char *obsf)
{
	FILE *obsfp = NULL;
	obsfp = fopen(obsf, "r");
	if (obsfp == NULL)
	{
		return 0;
	}

	char sel[21] = "END OF HEADER       ";
	unsigned int len = Obs_FileLineLength_Rinex3 + 1;
	char strline[Obs_FileLineLength_Rinex3 + 1] = "";
	if (SetPointerPos(obsfp, sel, strline, len) == 0)
	{
		fclose(obsfp);
		return NULL;
	}

	m_EpochCount = 0;
	m_FirstEpoch = NULL; m_LastEpoch = NULL;

	for (int i = 0; i<OBS_MAX_TYPE; ++i)
	{
		m_DataType.m_TypeIndex[i] = -1;
	}

	return obsfp;
}

/// Search satellite data by prn
w_SatelliteData* ObsData::Search_Satellite_by_Prn(unsigned int prn)
{
	w_SatelliteData* p = First_Sat;
	w_SatelliteData* re = NULL;
	while (p != NULL)
	{
		if (p->m_Prn == prn)
		{
			re = p;
			break;
		}
		p = p->m_Next;
	}
	return re;
}

int ObsData::Add_Satellite_to_Obs(SatelliteData *ptr_sat)
{
	if (ptr_sat == NULL) return 0;

	w_SatelliteData *s_p = Search_Satellite_by_Prn(ptr_sat->m_Prn);

	if (s_p == NULL) // a new satellite
	{
		w_SatelliteData *temp_ptr = (w_SatelliteData *)(malloc(sizeof(w_SatelliteData)));
		memset((void*)(temp_ptr), 0, sizeof(w_SatelliteData));

		temp_ptr->m_Prn = ptr_sat->m_Prn;
		temp_ptr->Sat_First_Data = ptr_sat;
		temp_ptr->Sat_Last_Data = ptr_sat;
		temp_ptr->m_Next = NULL;
		temp_ptr->m_EpochCount++;

		if (First_Sat == NULL)
		{
			First_Sat = temp_ptr;
			Last_Sat = temp_ptr;
		}
		else
		{
			Last_Sat->m_Next = temp_ptr;
			Last_Sat = temp_ptr;
		}

		m_SatCount++;
	}
	else             // a new epoch data
	{
		s_p->Sat_Last_Data->m_Next = ptr_sat;
		ptr_sat->m_Previous = s_p->Sat_Last_Data;
		s_p->Sat_Last_Data = ptr_sat;
		s_p->m_EpochCount++;
	}

	return 1;
}

/// Generate Observation Type Index of Rinex2
int ObsData::Generate_ObsTypeIndex2()
{
	unsigned int n = 0;
	unsigned int i = 0;
	int index = -1;

	m_DataType.m_Version = m_Header.m_RinexVersion;

	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, L1);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, L2);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C1);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C2);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, P1);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, P2);

	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L1);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L2);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C1);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C2);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, P1);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, P2);

	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L1);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L2);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C1);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C2);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, P1);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, P2);

	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L1);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L2);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L7);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C1);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C2);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C7);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, P1);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, P2);

	return 1;
}

/// Generate Observation Type of Rinex3
int ObsData::Generate_ObsTypeIndex3()
{
	unsigned int n = 0;
	unsigned int i = 0;
	int index = -1;

	m_DataType.m_Version = m_Header.m_RinexVersion;

	// For GPS obs type generation
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C1C);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, L1C);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, D1C);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C1W);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, L1W);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C2W);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, L2W);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, D2W);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C1P);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C2P);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C2X);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, L2X);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, L5Q);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C5Q);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, L5I);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C5I);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, L2L);
	CHECK_TYPE_INDEX(GPS, MIN_GPS_SATNO - 1, C2L);


	// For GLO obs type generation
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C1C);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L1C);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C1P);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L1P);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C2P);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L2P);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C2C);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L2C);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C4A);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L4A);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C4B);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L4B);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C4X);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L4X);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C6A);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L6A);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C6B);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L6B);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C6X);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L6X);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C3I);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L3I);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C3Q);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L3Q);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, C3X);
	CHECK_TYPE_INDEX(GLO, MIN_GLO_SATNO - 1, L3X);


	// For GAL obs type generation
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C1X);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L1X);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C1C);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L1C);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C5I);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L5I);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C5X);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L5X);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C5Q);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L5Q);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C7I);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L7I);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C7Q);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L7Q);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C7X);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L7X);

	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C8I);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L8I);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C8Q);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L8Q);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C8X);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L8X);

	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C6B);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L6B);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C6C);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L6C);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, C6X);
	CHECK_TYPE_INDEX(GAL, MIN_GAL_SATNO - 1, L6X);


	// For BDS obs type generation
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C2I);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L2I);
	/*CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C2Q);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L2Q);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C2X);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L2X);*///B1
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C1P);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L1P);//B1c
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C1D);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L1D);//B1c
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C1X);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L1X);//B1c
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C5D);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L5D);//B2a
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C5P);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L5P);//B2a
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C5X);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L5X);//B2a
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C7D);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L7D);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C7P);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L7P);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C7Z);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L7Z); //B2b
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C7I);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L7I);//B2
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C7Q);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L7Q);//B2
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C7X);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L7X);//B2
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C6I);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L6I);
	/*CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C6Q);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L6Q);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C6X);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L6X);*///B3
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C6D);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L6D);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C6P);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L6P);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, C6Z);
	CHECK_TYPE_INDEX(BDS, MIN_BDS_SATNO - 1, L6Z);//B3A


	return 1;
}

/// Generate Observation Type Index
int ObsData::Generate_ObsTypeIndex()
{
	if (m_Header.m_RinexVersion < 3)
		Generate_ObsTypeIndex2();
	else if (m_Header.m_RinexVersion < 4)
		Generate_ObsTypeIndex3();
	else
		return 0;

	return 1;
}

/// Read obs epoch data in rinex2 version
int ObsData::ReadObsEpochData_Rinex2(FILE *fp)
{
	unsigned int len = Obs_FileLineLength_Rinex2 + 1;
	char strline[Obs_FileLineLength_Rinex2 + 1] = "";
	GetFileLine(fp, strline, len);

	if (is_eof(fp))	return 0;

	// New list node.
	ObsEpochData * temp_ptr = (ObsEpochData *)malloc(sizeof(ObsEpochData));
	memset((void*)(temp_ptr), 0, sizeof(ObsEpochData));
	temp_ptr->m_Previous = NULL;
	temp_ptr->m_Next = NULL;

	// Get Epoch Time
	Get_Obs_Epochtime_Rinex2(strline, len, temp_ptr);

	// Get event flag
	char eventf[2] = "", count[4] = "";
	GetSubStr(strline, len, eventf, 2, 28, 1);
	int eventflag = str_to_i(eventf);

	// Get satellite number of cunrrent epoch
	GetSubStr(strline, len, count, 4, 29, 3);
	temp_ptr->m_SatCount = (unsigned int)(str_to_i(count));
	temp_ptr->m_EventFlag = (unsigned int)(eventflag);

	if (eventflag <= 1)
	{
		// Allocate memory to save temporary satellite list
		int *prn = (int *)(malloc(sizeof(int)*temp_ptr->m_SatCount));
		unsigned int i = 0, j = 0;
		for (; i<temp_ptr->m_SatCount; ++i)
		{
			char prnstr[4] = "";//get prn string,for example: G20、R01、C12
			GetSubStr(strline, len, prnstr, 4, 32 + j * 3, 3);
			prn[i] = GetSatNo(prnstr, 4, prnstr[0]);//'G'
			if ((++j) == 12){ j = 0; }
			//if number of satellites greater than 12,we must read another line data.
			if (((i + 1) % 12) == 0 && (i + 1) != temp_ptr->m_SatCount)
			{
				GetFileLine(fp, strline, len);
			}
		}

		// Statistic valid satellite nubmer
		int real_index = 0, valid_count = 0;
		for (i = 0; i<temp_ptr->m_SatCount; ++i)
		{
			if (prn[i])
			{
				valid_count++;
			}
		}

		// Allocate memory to save all satellite data
		temp_ptr->EpochData = (SatelliteData*)(malloc(sizeof(SatelliteData)*valid_count));
		memset((void*)(temp_ptr->EpochData), 0, sizeof(SatelliteData)*valid_count);

		// read data of whole satellites

		// Observation type number when read rinex file
		int obstypenum = m_Header.m_DataTypeList.m_DataTypeCount;

		// n-是每颗卫星数据所占的行数
		unsigned int n = (m_Header.m_DataTypeList.m_ObsTypeNum - 1) / 5 + 1;
		char **tempstr = (char **)(malloc(sizeof(char *)*n));
		for (i = 0; i<n; ++i)
		{
			tempstr[i] = (char *)(malloc(sizeof(char)*(Obs_FileLineLength_Rinex2 + 1)));
		}
		// Loop for each satellite record
		for (i = 0; i<temp_ptr->m_SatCount; ++i)
		{
			// Put each satellite data to data block
			unsigned int y = 0;
			for (y = 0; y<n; ++y)
			{
				GetFileLine(fp, tempstr[y], len);
			}

			if (prn[i] == 0)
			{
				continue;
			}

			// Allocate memory for each satellite observation value
			temp_ptr->EpochData[real_index].m_Data = (SatelliteDataValue*)(malloc(sizeof(SatelliteDataValue)*obstypenum));
			memset((void*)(temp_ptr->EpochData[real_index].m_Data), 0, sizeof(SatelliteDataValue)*obstypenum);

			temp_ptr->EpochData[real_index].m_Tow = temp_ptr->m_EpochTime.m_SecondofWeek;
			temp_ptr->EpochData[real_index].m_Week = temp_ptr->m_EpochTime.m_Week;
			temp_ptr->EpochData[real_index].m_Jd = temp_ptr->m_EpochTime.m_Jd;

			int q = 0;
			for (; q<obstypenum; ++q)
			{
				int type_t = m_Header.m_DataTypeList.m_Index[q];
				int k = type_t % 5, m = type_t / 5;
				char data[15] = "", lli[2] = "", strength[2] = "";
				GetSubStr(tempstr[m], len, data, 15, k * 16, 14);
				GetSubStr(tempstr[m], len, lli, 2, (k + 1) * 16 - 2, 1);
				GetSubStr(tempstr[m], len, strength, 2, (k + 1) * 16 - 1, 1);
				temp_ptr->EpochData[real_index].m_Data[q].m_Obs = str_to_f(data);
				temp_ptr->EpochData[real_index].m_Data[q].m_LLI = str_to_i(lli);
				temp_ptr->EpochData[real_index].m_Data[q].m_Strength = str_to_i(strength);
			}
			temp_ptr->EpochData[real_index].m_Prn = prn[i];
			Add_Satellite_to_Obs(&temp_ptr->EpochData[real_index]);
			++real_index;
		}
		temp_ptr->m_SatCount = real_index;

		if (valid_count == 0)
		{
			free(temp_ptr->EpochData);
			free(temp_ptr);
		}
		else
		{
			if (m_FirstEpoch == NULL)
			{
				m_FirstEpoch = temp_ptr;
				m_LastEpoch = temp_ptr;
			}
			else
			{
				temp_ptr->m_Previous = m_LastEpoch;
				m_LastEpoch->m_Next = temp_ptr;
				m_LastEpoch = temp_ptr;
			}
			++m_EpochCount;
		}

		// Release temporary memory
		free(prn);
		for (i = 0; i<n; ++i)
		{
			free(tempstr[i]);
		}
		free(tempstr);
	}
	else
	{
		// skip these lines
		unsigned int i = 0;
		for (; i<temp_ptr->m_SatCount; ++i)
		{
			GetFileLine(fp, strline, len);
		}
		free(temp_ptr);
	}

	return 1;
}

/// Read Rinex 3.0 Epoch Observation Data
int ObsData::ReadObsEpochData_Rinex3(FILE *fp)
{
	unsigned int len = Obs_FileLineLength_Rinex3 + 1;
	char strline[Obs_FileLineLength_Rinex3 + 1] = "";

	while (strline[0] != '>' && !is_eof(fp))
	{
		GetFileLine(fp, strline, len);
	}

	if (is_eof(fp))	return 0;

	// New list node.
	ObsEpochData * temp_ptr = (ObsEpochData *)malloc(sizeof(ObsEpochData));
	memset((void*)(temp_ptr), 0, sizeof(ObsEpochData));
	temp_ptr->m_Previous = NULL;
	temp_ptr->m_Next = NULL;

	// Get Rinex 3.0 Epoch Time
	Get_Obs_Epochtime_Rinex3(strline, len, temp_ptr);

	// Get event flag
	char epochf[2] = "", count[4] = "";
	GetSubStr(strline, len, epochf, 2, 31, 1);
	int epochflag = str_to_i(epochf);

	// Get satellite number of cunrrent epoch
	GetSubStr(strline, len, count, 4, 32, 3);
	temp_ptr->m_SatCount = (unsigned int)(str_to_i(count));
	temp_ptr->m_EventFlag = (unsigned int)(epochflag);

	if (epochflag <= 1)
	{
		unsigned int i = 0, j = 0;
		// Allocate memory to save temporary satellite list
		int *prn = (int *)(malloc(sizeof(int)*temp_ptr->m_SatCount));
		// Statistic valid satellite nubmer
		int real_index = 0, valid_count = 0;

		// Allocate memory to save all satellite data
		temp_ptr->EpochData = (SatelliteData*)(malloc(sizeof(SatelliteData)*temp_ptr->m_SatCount));
		memset((void*)(temp_ptr->EpochData), 0, sizeof(SatelliteData)*temp_ptr->m_SatCount);

		// Observation type number when read rinex file
		int obstypenum_gps = m_Header.m_DataTypeList_GPS.m_DataTypeCount;
		int obstypenum_glo = m_Header.m_DataTypeList_GLO.m_DataTypeCount;
		int obstypenum_gal = m_Header.m_DataTypeList_GAL.m_DataTypeCount;
		int obstypenum_bds = m_Header.m_DataTypeList_BDS.m_DataTypeCount;

		for (; i<temp_ptr->m_SatCount; ++i)
		{
			char prnstr[4] = "";//get prn string,for example: G20、R01、C12
			GetFileLine(fp, strline, len);
			GetSubStr(strline, len, prnstr, 4, j * 3, 3);
			int prn_tmp = GetSatNo(prnstr, 4, prnstr[0]);
			if (!prn_tmp) continue;
			prn[i] = prn_tmp;
			valid_count++;

			int obstypenum = 0;
			// Allocate memory for each satellite observation value
			switch (prnstr[0])
			{
				case 'G':
				{
					obstypenum = obstypenum_gps;
					break;
				}
				case 'R':
				{
					obstypenum = obstypenum_glo;
					break;
				}
				case 'E':
				{
					obstypenum = obstypenum_gal;
					break;
				}
				case 'C':
				{
					obstypenum = obstypenum_bds;
					break;
				}
				default:break;
			}

			temp_ptr->EpochData[real_index].m_Data = (SatelliteDataValue*)(malloc(sizeof(SatelliteDataValue)*obstypenum));
			memset((void*)(temp_ptr->EpochData[real_index].m_Data), 0, sizeof(SatelliteDataValue)*obstypenum);

			temp_ptr->EpochData[real_index].m_Tow = temp_ptr->m_EpochTime.m_SecondofWeek;
			temp_ptr->EpochData[real_index].m_Week = temp_ptr->m_EpochTime.m_Week;
			temp_ptr->EpochData[real_index].m_Jd = temp_ptr->m_EpochTime.m_Jd;

			int q = 0;
			for (; q<obstypenum; ++q)
			{
				char data[15] = "", lli[2] = "", strength[2] = "";
				GetSubStr(strline, len, data, 15, 3 + 16 * q, 14);
				GetSubStr(strline, len, lli, 2, 17 + q * 16, 1);
				GetSubStr(strline, len, strength, 2, 18 + q * 16 - 1, 1);
				temp_ptr->EpochData[real_index].m_Data[q].m_Obs = str_to_f(data);
				temp_ptr->EpochData[real_index].m_Data[q].m_LLI = str_to_i(lli);
				temp_ptr->EpochData[real_index].m_Data[q].m_Strength = str_to_i(strength);
			}
			temp_ptr->EpochData[real_index].m_Prn = prn[i];

			Add_Satellite_to_Obs(&temp_ptr->EpochData[real_index]);

			real_index++;
		}
		temp_ptr->m_SatCount = real_index;

		if (valid_count == 0)
		{
			free(temp_ptr->EpochData);
			free(temp_ptr);
		}
		else
		{
			if (m_FirstEpoch == NULL)
			{
				m_FirstEpoch = temp_ptr;
				m_LastEpoch = temp_ptr;
			}
			else
			{
				temp_ptr->m_Previous = m_LastEpoch;
				m_LastEpoch->m_Next = temp_ptr;
				m_LastEpoch = temp_ptr;
			}
			++m_EpochCount;
		}

		// Release temporary memory
		free(prn);
	}
	else
	{
		// skip these lines
		unsigned int i = 0;
		for (; i<temp_ptr->m_SatCount; ++i)
		{
			GetFileLine(fp, strline, len);
		}
		free(temp_ptr);
	}

	return 1;
}

/// Read Observation Epoch Data
int ObsData::ReadObsEpochData(FILE *fp)
{
	if (m_Header.m_RinexVersion < 3)
		ReadObsEpochData_Rinex2(fp);
	else if (m_Header.m_RinexVersion < 4)
		ReadObsEpochData_Rinex3(fp);
	else
		return 0;

	return 1;
}

// Read rinex2 obs file body
int ObsData::ReadObsFileBody_Rinex2(char *filename)
{
	FILE *fp = NULL;
	fp = fopen(filename, "r");
	if (fp == NULL)
	{
		return 0;
	}

	char sel[21] = "END OF HEADER       ";                //Skip Header.
	unsigned int len = Obs_FileLineLength_Rinex2 + 1;
	char strline[Obs_FileLineLength_Rinex2 + 1] = "";
	if (SetPointerPos(fp, sel, strline, len) == 0)
	{
		fclose(fp);
		return 0;
	}

	m_EpochCount = 0;
	m_FirstEpoch = NULL;
	m_LastEpoch = NULL;

	while (ReadObsEpochData_Rinex2(fp))
	{
		int i = 0;
		for (i = 0; i<OBS_MAX_TYPE; ++i)
		{
			m_DataType.m_TypeIndex[i] = -1;
		}

		ObsEpochData *ptr = m_FirstEpoch;

		while (ptr != NULL)
		{
			ptr->m_ObsType = &m_DataType;
			ptr = ptr->m_Next;
		}

		//m_FirstEpoch->Display_EpochData();
	}
	fclose(fp);

	Generate_ObsTypeIndex2();

	return 1;
}

/// Read rinex3 file body
int ObsData::ReadObsFileBody_Rinex3(char *filename)
{
	FILE *fp = NULL;
	fp = fopen(filename, "r");
	if (fp == NULL)
	{
		return 0;
	}

	char sel[21] = "END OF HEADER       ";                //Skip Header.
	unsigned int len = Obs_FileLineLength_Rinex3 + 1;
	char strline[Obs_FileLineLength_Rinex3 + 1] = "";
	if (SetPointerPos(fp, sel, strline, len) == 0)
	{
		fclose(fp);
		return 0;
	}

	m_EpochCount = 0;
	m_FirstEpoch = NULL; m_LastEpoch = NULL;

	while (ReadObsEpochData_Rinex3(fp))//正式读观测文件  
	{
		int i = 0;
		for (i = 0; i<OBS_MAX_TYPE; ++i)
		{
			m_DataType.m_TypeIndex[i] = -1;
		}

		ObsEpochData *ptr = m_FirstEpoch;

		while (ptr != NULL)
		{
			ptr->m_ObsType = &m_DataType;
			ptr = ptr->m_Next;
		}

	}
	fclose(fp);

	Generate_ObsTypeIndex3();

	return 1;
}

/// Read obs file body
int ObsData::ReadObsFileBody(char *filename)
{
	if (m_Header.m_RinexVersion < 3)
		ReadObsFileBody_Rinex2(filename);
	else if (m_Header.m_RinexVersion < 4)
		ReadObsFileBody_Rinex3(filename);
	else
		return 0;

	return 1;
}

/// Read Observation Rinex File
int ObsData::ReadObsRinex(char* filename)
{
	// Read obs file header
	if (m_Header.ReadObsFileHeader(filename) == 0)
		return 0;

	// Read obs file body
	if (ReadObsFileBody(filename) == 0)
		return 0;

	return 1;
}

int ObsData::ReadObsRinex_EpochWise_2(FILE* fp)
{
	if (fp == NULL)
		return 0;
	/* Get observation data epoch wise ----------------------------------------------- */
	unsigned int len = Obs_FileLineLength_Rinex2 + 1;
	char strline[Obs_FileLineLength_Rinex2 + 1] = "";
	GetFileLine(fp, strline, len);

	if (is_eof(fp))
		return 0;

	// New list node.
	ObsEpochData * temp_ptr = (ObsEpochData *)malloc(sizeof(ObsEpochData));
	memset((void*)(temp_ptr), 0, sizeof(ObsEpochData));
	temp_ptr->m_Previous = NULL;
	temp_ptr->m_Next = NULL;

	// Get Epoch Time
	Get_Obs_Epochtime_Rinex2(strline, len, temp_ptr);

	// Get event flag
	char eventf[2] = "", count[4] = "";
	GetSubStr(strline, len, eventf, 2, 28, 1);
	int eventflag = str_to_i(eventf);

	// Get satellite number of cunrrent epoch
	GetSubStr(strline, len, count, 4, 29, 3);
	temp_ptr->m_SatCount = (unsigned int)(str_to_i(count));
	temp_ptr->m_EventFlag = (unsigned int)(eventflag);

	if (eventflag <= 1)
	{
		// Allocate memory to save temporary satellite list
		int *prn = (int *)(malloc(sizeof(int)*temp_ptr->m_SatCount));
		unsigned int i = 0, j = 0;
		for (; i<temp_ptr->m_SatCount; ++i)
		{
			char prnstr[4] = "";//get prn string,for example: G20、R01、C12
			GetSubStr(strline, len, prnstr, 4, 32 + j * 3, 3);
			prn[i] = GetSatNo(prnstr, 4, 'G');
			if (prn[i] == 0)
				prn[i] = GetSatNo(prnstr, 4, 'R');
			if ((++j) == 12){ j = 0; }
			//if number of satellites greater than 12,we must read another line data.
			if (((i + 1) % 12) == 0 && (i + 1) != temp_ptr->m_SatCount)
			{
				GetFileLine(fp, strline, len);
			}
		}

		// Statistic valid satellite nubmer
		int real_index = 0, valid_count = 0;
		for (i = 0; i<temp_ptr->m_SatCount; ++i)
		{
			if (prn[i])
			{
				valid_count++;
			}
		}

		// Allocate memory to save all satellite data
		temp_ptr->EpochData = (SatelliteData*)(malloc(sizeof(SatelliteData)*valid_count));
		memset((void*)(temp_ptr->EpochData), 0, sizeof(SatelliteData)*valid_count);

		// read data of whole satellites

		// Observation type number when read rinex file
		int obstypenum = m_Header.m_DataTypeList.m_DataTypeCount;

		// n-是每颗卫星数据所占的行数
		unsigned int n = (m_Header.m_DataTypeList.m_ObsTypeNum - 1) / 5 + 1;
		char **tempstr = (char **)(malloc(sizeof(char *)*n));
		for (i = 0; i<n; ++i)
		{
			tempstr[i] = (char *)(malloc(sizeof(char)*(Obs_FileLineLength_Rinex2 + 1)));
		}
		// Loop for each satellite record
		for (i = 0; i<temp_ptr->m_SatCount; ++i)
		{
			// Put each satellite data to data block
			unsigned int y = 0;
			for (y = 0; y<n; ++y)
			{
				GetFileLine(fp, tempstr[y], len);
			}
			if (prn[i] == 0)
			{
				continue;
			}

			// Allocate memory for each satellite observation value
			temp_ptr->EpochData[real_index].m_Data = (SatelliteDataValue*)(malloc(sizeof(SatelliteDataValue)*obstypenum));
			memset((void*)(temp_ptr->EpochData[real_index].m_Data), 0, sizeof(SatelliteDataValue)*obstypenum);

			temp_ptr->EpochData[real_index].m_Tow = temp_ptr->m_EpochTime.m_SecondofWeek;
			temp_ptr->EpochData[real_index].m_Week = temp_ptr->m_EpochTime.m_Week;
			temp_ptr->EpochData[real_index].m_Jd = temp_ptr->m_EpochTime.m_Jd;

			int q = 0;
			for (; q<obstypenum; ++q)
			{
				int type_t = m_Header.m_DataTypeList.m_Index[q];
				int k = type_t % 5, m = type_t / 5;
				char data[15] = "", lli[2] = "", strength[2] = "";
				GetSubStr(tempstr[m], len, data, 15, k * 16, 14);
				GetSubStr(tempstr[m], len, lli, 2, (k + 1) * 16 - 2, 1);
				GetSubStr(tempstr[m], len, strength, 2, (k + 1) * 16 - 1, 1);
				temp_ptr->EpochData[real_index].m_Data[q].m_Obs = str_to_f(data);
				temp_ptr->EpochData[real_index].m_Data[q].m_LLI = str_to_i(lli);
				temp_ptr->EpochData[real_index].m_Data[q].m_Strength = str_to_i(strength);
			}
			temp_ptr->EpochData[real_index].m_Prn = prn[i];
			Add_Satellite_to_Obs(&temp_ptr->EpochData[real_index]);
			++real_index;
		}
		temp_ptr->m_SatCount = real_index;

		if (valid_count == 0)
		{
			free(temp_ptr->EpochData);
			free(temp_ptr);
			return 0;
		}

		// Release temporary memory
		free(prn);
		for (i = 0; i<n; ++i)
		{
			free(tempstr[i]);
		}
		free(tempstr);
	}
	else
	{
		// skip these lines
		unsigned int i = 0;
		for (; i<temp_ptr->m_SatCount; ++i)
		{
			GetFileLine(fp, strline, len);
		}
		free(temp_ptr);

		return 0;
	}

	/* ------------------------------- */
	m_FirstEpoch = temp_ptr;
	m_DataType.m_TypeNum = m_Header.m_DataTypeList.m_DataTypeCount;
	m_FirstEpoch->m_ObsType = &m_DataType;
	Generate_ObsTypeIndex2();
	/* ------------------------------- */

	return 1;
}

int ObsData::ReadObsRinex_EpochWise_3(FILE* fp)
{
	if (fp == NULL)	return 0;

	unsigned int len = Obs_FileLineLength_Rinex3 + 1;
	char strline[Obs_FileLineLength_Rinex3 + 1] = "";

	if (is_eof(fp))	return 0;

	m_EpochCount = 0;
	m_FirstEpoch = NULL; m_LastEpoch = NULL;

	if (ReadObsEpochData_Rinex3(fp))
	{
		int i = 0;
		for (i = 0; i<OBS_MAX_TYPE; ++i)
		{
			m_DataType.m_TypeIndex[i] = -1;
		}

		ObsEpochData *ptr = m_FirstEpoch;

		while (ptr != NULL)
		{
			ptr->m_ObsType = &m_DataType;
			ptr = ptr->m_Next;
		}
	}

	Generate_ObsTypeIndex3();

	return 1;
}

/// Read observation epoch-wise
int ObsData::ReadObsRinex_EpochWise(FILE* fp)
{
	int re = 0;

	if (m_Header.m_RinexVersion < 3.0)
		re = this->ReadObsRinex_EpochWise_2(fp);
	else if (m_Header.m_RinexVersion < 4.0)
		re = this->ReadObsRinex_EpochWise_3(fp);

	return re;
}


int ObsData::Clear()
{
	m_Header.Clear();

	ObsEpochData *ptr = m_FirstEpoch;
	while (ptr != NULL)
	{
		ObsEpochData *p = ptr->m_Next;

		unsigned int i = 0;
		for (; i < ptr->m_SatCount; ++i)
		{
			free(ptr->EpochData[i].m_Data);
		}

		ptr->m_SatCount = 0;
		free(ptr->EpochData);
		free(ptr);

		ptr = p;
	}

	m_FirstEpoch = NULL; m_LastEpoch = NULL;
	m_EpochCount = NULL;

	w_SatelliteData *s_p = First_Sat;
	while (s_p != NULL)
	{
		w_SatelliteData *s_p_t = s_p->m_Next;
		s_p->m_Prn = 0;
		s_p->m_EpochCount = 0;
		s_p->Sat_First_Data = NULL;
		s_p->Sat_Last_Data = NULL;

		free(s_p);
		s_p = s_p_t;
	}
	First_Sat = NULL;
	Last_Sat = NULL;
	m_SatCount = 0;

	return 1;
}
