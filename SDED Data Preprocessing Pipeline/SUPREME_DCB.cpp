#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SUPREME_DCB.h"
#include "SUPREME_CommonFunction.h"

#define gnss_dcb_line 100

CodeDCB::CodeDCB()
{
	string a[28] = {"C1C  C1W", "C2C  C2W", "C2W  C2S", "C2W  C2L", "C2W  C2X","C1C  C2W","C1C  C5Q","C1C  C5X","C1W  C2W","C1C  C1P","C2C  C2P","C1C  C2C"," C1C  C2P","C1P  C2P","C1C  C6C","C1C  C7Q","C1C  C8Q","C1X  C5X","C1X  C7X","C1X  C8X","C1P  C5P","C1D  C5D","C1X  C6I","C1P  C6I","C1D  C6I","C2I  C6I","C1X  C7Z","C2I  C7I"};
	for (int i = 0; i < 28; i++)
	{
		DCB_mgex_type[i] = a[i];
	}
	//
	for (int i = 0; i < 28; i++)
	{
		a[i].clear();
	}
}
CodeDCB::~CodeDCB()
{
	unsigned int i = 0;
	DCB_Item *ptr1   = m_DCB_First;
	DCB_Item* ptr111 = m_mgex_DCB_First1;
	while (ptr1 != NULL)
	{
		DCB_Item *p = ptr1->m_Next;
		free(ptr1);
		ptr1 = p;
	}
	while (ptr111 != NULL)
	{
		DCB_Item* p11 = ptr111->m_Next;
		free(ptr111);
		ptr111 = p11;
	}
	if (m_IndexTable) free(m_IndexTable);
	if (m_IndexTable_mgex1) free(m_IndexTable_mgex1);
	m_Size = 0;
	m_TableSize = 0;
	m_Size_mgex1 = 0;
	m_TableSize_mgex1 = 0;
	for (int i = 0; i < 28; i++)
	{
		DCB_mgex_type[i].clear();
	}

}

DCB_Item* dcb_mgex_find_data(CodeDCB* dcb, unsigned int prn)
{
	DCB_Item* p; p = dcb->m_mgex_DCB_First1;
	if (dcb->m_mgex_DCB_First1 == NULL )
	{
		return NULL;
	}

	if (prn == dcb->m_mgex_DCB_First1->m_Prn) 
	{
		return dcb->m_mgex_DCB_First1;
	}
	else
	{
		while (p!=NULL)
		{
			if (prn == p->m_Prn) 
			{
				return p; 
			}
			else
			{
				p = p->m_Next;
			}
		}
	}
	return NULL;
	
}
DCB_Item *dcb_find_data(CodeDCB *dcb, unsigned int prn)
{
	if (dcb->m_DCB_First == NULL || dcb->m_IndexTable == NULL)
	{
		return NULL;
	}

	int i = 0, j = dcb->m_TableSize - 1;
	int d = abs(j - i);

	while (d>1)
	{
		unsigned int k = d / 2 + i;
		if (prn<dcb->m_IndexTable[k]->m_Prn)
		{
			j = k;
		}
		else
		{
			i = k;
		}
		d = abs(j - i);
	}

	if (prn == dcb->m_IndexTable[i]->m_Prn)
	{
		return dcb->m_IndexTable[i];
	}
	else if (prn == dcb->m_IndexTable[j]->m_Prn)
	{
		return dcb->m_IndexTable[j];
	}
	else
	{
		return NULL;
	}
}
int dcb_mgex_add_item(unsigned int prn, double value, char *flag, CodeDCB* dcb, DCB_Item** tail,string type[])
{
	if (prn == 0)
	{
		return 0;
	}
	DCB_Item* ptr1;
	int find_result = 0;
	
	ptr1 = dcb_mgex_find_data(dcb, prn);

	if (ptr1 != NULL) 
	{
		find_result = 1;
	}
	else
	{
		ptr1 = (DCB_Item*)(malloc(sizeof(DCB_Item)));
		memset((void*)(ptr1), 0, sizeof(DCB_Item));
	}
	for (int i = 0; i < 28; i++)
	{
		if (strstr(flag, type[i].c_str()))
		{
			ptr1->DCB_mgex[i]= value; ptr1->m_Prn = prn;
			break;
		}
	}
	
	if (dcb->m_Size_mgex1 == 0)
	{
		dcb->m_mgex_DCB_First1 = ptr1;
		*tail = ptr1;
		if (find_result != 1) {
		++dcb->m_Size_mgex1;
		}
	}
	else if(find_result != 1)
	{
		if (ptr1->m_Prn > (*tail)->m_Prn)
		{
			(*tail)->m_Next = ptr1;
			*tail = ptr1;
			if (find_result != 1) {
				++dcb->m_Size_mgex1;
			}
		}
		else
		{
			DCB_Item* p = dcb->m_mgex_DCB_First1;
			DCB_Item* q = dcb->m_mgex_DCB_First1;
			while (p != NULL)
			{
				if (ptr1->m_Prn < p->m_Prn)
				{
					break;
				}
				q = p;
				p = p->m_Next;
			}
			if (p == dcb->m_mgex_DCB_First1)
			{
				ptr1->m_Next = dcb->m_mgex_DCB_First1;
				dcb->m_mgex_DCB_First1 = ptr1;
				if (find_result != 1) {
					++dcb->m_Size_mgex;
				}
			}
			else
			{
				ptr1->m_Next = q->m_Next;
				q->m_Next = ptr1;
				if (find_result != 1) {
					++dcb->m_Size_mgex1;
				}
			}
		}
	}
	return 1;
}

int dcb_add_item(unsigned int prn, double value, char flag, CodeDCB* dcb, DCB_Item** tail)
{
	if (prn == 0)
	{
		return 0;
	}

	int find_result = 0;

	DCB_Item* ptr = dcb_find_data(dcb, prn);
	if (ptr != NULL)
	{
		find_result = 1;
	}
	else
	{
		ptr = (DCB_Item*)(malloc(sizeof(DCB_Item))); //申请内存了
		memset((void*)(ptr), 0, sizeof(DCB_Item));
	}

	if (flag == 'C')
	{
		ptr->m_P1C1 = value;
	}
	else if (flag == 'P')
	{
		ptr->m_P1P2 = value;
	}
	else if (flag == 'B') 
	{
		ptr->m_P2C2 = value;
	}
	if (find_result == 1)
	{
		return 1;
	}

	ptr->m_Prn = prn;

	if (dcb->m_Size == 0)
	{
		dcb->m_DCB_First = ptr;//  接受链表头
		*tail = ptr;           //  把指针指到与链表头地址一致
		++dcb->m_Size;
	}
	else
	{
		if (ptr->m_Prn > (*tail)->m_Prn)
		{
			(*tail)->m_Next = ptr;
			*tail = ptr;
			++dcb->m_Size;
		}
		else
		{
			DCB_Item* p = dcb->m_DCB_First;
			DCB_Item* q = dcb->m_DCB_First;
			while (p != NULL)
			{
				if (ptr->m_Prn < p->m_Prn)
				{
					break;
				}
				q = p;
				p = p->m_Next;
			}
			if (p == dcb->m_DCB_First)
			{
				ptr->m_Next = dcb->m_DCB_First;
				dcb->m_DCB_First = ptr;
				++dcb->m_Size;
			}
			else
			{
				ptr->m_Next = q->m_Next;
				q->m_Next = ptr;
				++dcb->m_Size;
			}
		}
	}

	return 1;
}


/// Read Dcb data
int CodeDCB::Read_DCB(char* filename)
{
	FILE* fp = fopen(filename, "r");

	if (fp == NULL)
	{
		printf("Error:DCB file can not be opened.\n");
		return 0;
	}

	char strline[gnss_dcb_line + 1] = "";
	unsigned int len = gnss_dcb_line + 1;
	unsigned int i = 0, j = 0;

	for (i = 0; i < 4; ++i)
	{
		GetFileLine(fp, strline, len);
		if (is_eof(fp) == 1)// come to end of file.
		{
			printf("Error:DCB file is not integrate.\n");
			return 0;
		}
	}
	char str_flag[6] = "";
	GetSubStr(strline, len, str_flag, 6, 14, 5);

	char file_type = 'C';

	if (strcmp(str_flag, "P1-P2") == 0)
	{
		file_type = 'P';
	}
	if (strcmp(str_flag, "P2-C2") == 0) 
	{
		file_type = 'B';
	}
	// skip 3 lines.
	for (i = 0; i < 3; ++i)
	{
		GetFileLine(fp, strline, len);
		if (is_eof(fp) == 1)// come to end of file.
		{
			printf("Error:DCB file is not integrate.\n");
			return 0;
		}
	}

	DCB_Item* tail = NULL;

	if (m_Size != 0 && m_DCB_First != NULL && m_IndexTable != NULL)
	{
		tail = m_IndexTable[m_Size - 1];
	}

	int l = 0;

	while (is_eof(fp) == 0)
	{
		char prn[4] = "";
		char val[10] = "";

		GetFileLine(fp, strline, len);

		GetSubStr(strline, len, prn, 4, 0, 3);
		GetSubStr(strline, len, val, 10, 26, 9);

		unsigned int prnu = GetSatNo(prn, 4, 0);
		double data = str_to_f(val);

		dcb_add_item(prnu, data, file_type, this, &tail);
	}
	// create index table.
	DCB_Item* p = m_DCB_First;
	if (p != NULL)
	{
		if (m_IndexTable != NULL)
		{
			free(m_IndexTable);
		}
		m_IndexTable = (DCB_Item**)(malloc(m_Size * sizeof(DCB_Item*)));
		memset((void*)(m_IndexTable), 0, m_Size * sizeof(DCB_Item*));
	}

	i = 0;
	while (p != NULL)
	{
		m_IndexTable[i] = p;
		++i;
		p = p->m_Next;
	}

	m_TableSize = m_Size;

	fclose(fp);
	return 1;
}

DCB_Item* CodeDCB::Find_DCB_Data(unsigned int prn)
{
	if (m_DCB_First == NULL || m_IndexTable == NULL)
		return NULL;

	int i = 0, j = m_TableSize - 1;
	int d = abs(j - i);

	while (d > 1)
	{
		unsigned int k = d / 2 + i;
		if (prn < m_IndexTable[k]->m_Prn)
			j = k;
		else
			i = k;
		d = abs(j - i);
	}

	if (prn == m_IndexTable[i]->m_Prn)
		return m_IndexTable[i];
	else if (prn == m_IndexTable[j]->m_Prn)
		return m_IndexTable[j];
	else
		return NULL;
}

/// Get p1c1 dcb
double CodeDCB::Get_DCB_P1C1(unsigned int prn)
{
	DCB_Item* p = Find_DCB_Data(prn);
	if (p == NULL)
		return 0.0;
	else
		return p->m_P1C1;
}

/// Get p1p2 dcb
double CodeDCB::Get_DCB_P1P2(unsigned int prn)
{
	DCB_Item* p = Find_DCB_Data(prn);
	if (p == NULL) return 0.0;
	else
		return p->m_P1P2;
}
double CodeDCB::Get_DCB_P2C2(unsigned int prn)
{
	DCB_Item* p = Find_DCB_Data(prn);
	if (p == NULL) return 0.0;
	else
		return p->m_P2C2;
}
double CodeDCB::Get_DCB_B1B3(unsigned int prn)
{
	DCB_Item* p = Find_mgex_DCB_Data(prn);
	if (p == NULL) return 0.0;
	else
		return p->m_B1B3;
}



DCB_Item* CodeDCB::Find_mgex_DCB_Data(unsigned int prn)
{
	if (m_mgex_DCB_First1 == NULL )
		return NULL;

	int i = 0, j = m_TableSize_mgex1 - 1;
	int d = abs(j - i);

	while (d > 1)
	{
		unsigned int k = d / 2 + i;
		if (prn < m_IndexTable_mgex1[k]->m_Prn)
			j = k;
		else
			i = k;
		d = abs(j - i);
	}

	if (prn == m_IndexTable_mgex1[i]->m_Prn)
		return m_IndexTable_mgex1[i];
	else if (prn == m_IndexTable_mgex1[j]->m_Prn)
		return m_IndexTable_mgex1[j];
	else
		return NULL;

}

int CodeDCB::Read_mgex_DCB(char* filename)
{
	FILE* fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("Error:mgex_DCB file can not be opened.\n");
		return 0;
	}
	char strline1[gnss_dcb_line + 10] = "";
	unsigned int len = gnss_dcb_line + 10;
	unsigned int i = 0, j = 0;

	while (!strstr(strline1, "+BIAS/SOLUTION"))
	{
		GetFileLine(fp, strline1, len);
	}
	DCB_Item* tail1 = NULL; //char file_type;
	if (m_Size_mgex1 != 0 && m_mgex_DCB_First1 != NULL && m_IndexTable_mgex1 != NULL)
	{
		tail1 = m_IndexTable_mgex1[m_Size_mgex - 1];
	}
	string type[28];
	for (int kk = 0; kk < 28; kk++) { type[kk] = DCB_mgex_type[kk]; }
	int l = 0;
	while (is_eof(fp) == 0)
	{
		GetFileLine(fp, strline1, len);
		if (strncmp(strline1 + 1, "DSB", 3) || strncmp(strline1 + 15, "    ", 4)) continue;
		char prn[4] = "";
		char val[10] = "";
		GetSubStr(strline1, len, prn, 4, 11, 3);
		GetSubStr(strline1, len, val, 10, 80, 9);

		unsigned int prnu = GetSatNo(prn, 4, 0);
		double data = str_to_f(val);
		
		char file_type[9]="";
		GetSubStr(strline1, len, file_type, 9, 25, 8);
		
		dcb_mgex_add_item(prnu, data, file_type, this, &tail1,type);
		
	}
	DCB_Item* p = m_mgex_DCB_First1;//p没有申请内存所以出了函数就直接消除了；
	if (p != NULL)
	{
		if (m_IndexTable_mgex1 != NULL)
		{
			free(m_IndexTable_mgex1);
		}

		m_IndexTable_mgex1 = (DCB_Item**)(malloc(m_Size_mgex1 * sizeof(DCB_Item*)));
		memset((void*)(m_IndexTable_mgex1), 0, m_Size_mgex1 * sizeof(DCB_Item*));
	}

	i = 0;
	while (p != NULL)
	{
		m_IndexTable_mgex1[i] = p;
		++i;
		p = p->m_Next;
	}

	m_TableSize_mgex1 = m_Size_mgex1;

	fclose(fp);
	return 1;
}
double find_Dcb_value(string type, DCB_Item* PtrDcb)
{
	string DCB_mgex_type[28] = { "C1C  C1W", "C2C  C2W", "C2W  C2S",
		"C2W  C2L", "C2W  C2X","C1C  C2W","C1C  C5Q","C1C  C5X",
		"C1W  C2W","C1C  C1P","C2C  C2P","C1C  C2C"," C1C  C2P",
		"C1P  C2P","C1C  C6C","C1C  C7Q","C1C  C8Q","C1X  C5X",
		"C1X  C7X","C1X  C8X","C1P  C5P","C1D  C5D","C1X  C6I",
		"C1P  C6I","C1D  C6I","C2I  C6I","C1X  C7Z","C2I  C7I" };
	for (int i = 0; i < 28; i++)
	{
		if (!(type.compare(DCB_mgex_type[i]))) 
		{
			return PtrDcb->DCB_mgex[i] * 1E-9 * 299792458.0;
		}
	}
	return 0;
}

