#include <string.h>

#include "SUPREME_STACRD.h"
#include "SUPREME_CommonFunction.h"

StationCrd::StationCrd()
{
	flag = 0;
	next = NULL;
}

/// Get station information (Receiver type, Antenna type, Antenna offset...)
int StationCrd::GetStaInfo(char *filename, gnsstime &t)
{
	char strline[500] = { 0 };
	char sta_name[5] = { 0 };
	char *name;

	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("Error: STA information file open failed.\n");
		return -1;
	}

	flag = 0;
	while (!feof(fp))
	{
		GetFileLine(fp, strline, 500);
		char strflag[50] = { 0 };
		memcpy(strflag, strline + 9, sizeof(char)* 20);

		if (strstr(strflag, "STATION INFORMATION"))
		{
			while (!feof(fp))
			{
				GetFileLine(fp, strline, 500);
				name = GetSubStr(strline, 0, 4);

				if (strcmp(name, m_StaName) == 0)
				{
					gnsstime t1, t2;
					t1.m_Year = (int)str2num(strline, 27, 4);
					t1.m_Month = (int)str2num(strline, 32, 2);
					t1.m_Day = (int)str2num(strline, 35, 2);
					t1.m_Hour = (int)str2num(strline, 38, 2);
					t1.m_Min = (int)str2num(strline, 41, 2);
					t1.m_Sec = str2num(strline, 44, 2);
					t1._date2gpst();

					t2.m_Year = (int)str2num(strline, 48, 4);
					t2.m_Month = (int)str2num(strline, 53, 2);
					t2.m_Day = (int)str2num(strline, 56, 2);
					t2.m_Hour = (int)str2num(strline, 59, 2);
					t2.m_Min = (int)str2num(strline, 62, 2);
					t2.m_Sec = str2num(strline, 65, 2);
					t2._date2gpst();

					if (t > t2 && t2.m_Year != 0) continue;	if (t < t1 && t1.m_Year != 0) continue;

					memcpy(m_RecType, strline + 69, sizeof(char)* 20);   // Receiver Type
					memcpy(m_AntType, strline + 121, sizeof(char)* 20);  // Antenna Type
					m_AntOffset[0] = str2num(strline, 173, 8);           // North
					m_AntOffset[1] = str2num(strline, 183, 8);           // East
					m_AntOffset[2] = str2num(strline, 193, 8);           // Up

					flag = 1;
					break;
				}
			}
		}
		if (flag == 1) break;
	}

	fclose(fp);

	if (flag == 0)
	{
		printf("Warning: No %s station information in %s file.\n", m_StaName, filename);
		return 0;
	}

	return 1;
}

/// Get station coordinate from SNX file
int StationCrd::GetStaSnx(char* filename)
{
	char strline[200] = { 0 };
	char sta_name[5] = { 0 };
	char *name;

	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("Error: SNX file open failed.\n");
		return -1;
	}

	memcpy(sta_name, m_StaName, 4);

	while (!feof(fp))
	{
		GetFileLine(fp, strline, 200);
		char strflag[50] = { 0 };

		memcpy(strflag, strline, sizeof(char)* 18);

		if (strcmp(strflag, "+SOLUTION/ESTIMATE") == 0)
		{
			while (!feof(fp))
			{
				GetFileLine(fp, strline, 200);
				name = GetSubStr(strline, 14, 4);

				if (strcmp(name, sta_name) == 0)
				{
					StaCoord.setX(str2num(strline, 46, 22));
					GetFileLine(fp, strline, 200);
					StaCoord.setY(str2num(strline, 46, 22));
					GetFileLine(fp, strline, 200);
					StaCoord.setZ(str2num(strline, 46, 22));

					flag = 1;
					break;
				}
				else continue;
			}
		}

		if (flag == 1)
			break;

		if (strcmp(strflag, "-SOLUTION/ESTIMATE") == 0)
			break;
	}

	fclose(fp);

	if (flag == 0)
	{
		printf("Warning: No %s coordinate information in %s file.\n", m_StaName, filename);
		return 0;
	}

	next = NULL;

	return 1;
}

/// Get station coordinate from CRD file
int StationCrd::GetStaCrd(const char* filename)
{
	char strline[200] = { 0 };
	char sta_name[5] = { 0 };
	char *name;

	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("Error: CRD file open failed.\n");
		return -1;
	}

	memcpy(sta_name, m_StaName, 4);

	while (!feof(fp))
	{
		GetFileLine(fp, strline, 200);

		name = GetSubStr(strline, 5, 4);
		//name = substr(strline, 6, 5);

		if (strcmp(name, sta_name) == 0)
		{
			StaCoord.setX(str_to_f(GetSubStr(strline, 21, 15)));
			StaCoord.setY(str_to_f(GetSubStr(strline, 36, 15)));
			StaCoord.setZ(str_to_f(GetSubStr(strline, 51, 15)));

			/*StaCoord.setX(str_to_f(get_substr(strline, 23, 15)));
			StaCoord.setY(str_to_f(get_substr(strline, 38, 15)));
			StaCoord.setZ(str_to_f(get_substr(strline, 53, 15)));*/

			flag = 1;
			break;
		}
		else continue;
	}

	fclose(fp);

	if (flag == 0)
	{
		printf("Warning: No %s coordinate information in %s file.\n", m_StaName, filename);
		return 0;
	}

	next = NULL;

	return 1;
}

/// Display station coordinate
int StationCrd::Display()
{
	if (flag)
		printf("%s: %16.5f%16.5f%16.5f\n", m_StaName, StaCoord._X, StaCoord._Y, StaCoord._Z);
	else
	{
		//printf("Get crd failed.\n");
	}

	return 1;
}