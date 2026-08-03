/* -------------------------------------------------------------------------
* SUPREME_DCB.h : SUPREME for DCB Input and Application
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2019.06.05
* ------------------------------------------------------------------------- */

#ifndef SUPREME_DCB_H_HH
#define SUPREME_DCB_H_HH
#include<string>
using namespace std;
class DCB_Item
{
public:
	unsigned int m_Prn;
	DCB_Item* m_Next;

	double m_P1C1;
	double m_P1P2;
	double m_P2C2;
	double m_B1B3;
	double DCB_mgex[28];
	
};

class CodeDCB
{
public:
	CodeDCB();
	~CodeDCB();
	
	string DCB_mgex_type[28]; 
	unsigned int m_Size;
	unsigned int m_TableSize;
	unsigned int m_Size_mgex;
	unsigned int m_TableSize_mgex;

	unsigned int m_Size_mgex1;
	unsigned int m_TableSize_mgex1;

	DCB_Item *m_DCB_First;
	DCB_Item **m_IndexTable;

	DCB_Item* m_mgex_DCB_First1;
	DCB_Item** m_IndexTable_mgex1;


	int Read_DCB(char *filename);
	DCB_Item* dcb_mgex_find_data1(CodeDCB* dcb, unsigned int prn);
	DCB_Item* dcb_mgex_find_data(CodeDCB* dcb, unsigned int prn);
	DCB_Item* Find_DCB_Data(unsigned int prn);

	DCB_Item* Find_mgex_DCB_Data(unsigned int prn);

	int Read_mgex_DCB(char* filename);

	double Get_DCB_P1C1(unsigned int prn);

	double Get_DCB_P1P2(unsigned int prn);

	double Get_DCB_B1B3(unsigned int prn);
	
	double Get_DCB_P2C2(unsigned int prn);
	
};

double find_Dcb_value(string type, DCB_Item* PtrDcb);




#endif
