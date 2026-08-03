#pragma once
#include<string>
#include<vector>
#include "SUPREME_GnssTime.h"
#include"SUPREME_CommonFunction.h"
using namespace std;

class ifcb_
{
public:
	//ifcb_();
	//~ifcb_();
	unsigned int Prn[200];
	double ifcb_value[200];  //这用指针就是不好使 妈的


};
class IFCB_Value
{
public:
	gnsstime time;
	ifcb_* IFCB_;
	
	int IFCB_size;
	IFCB_Value* m_Next;
};

class IFCB
{
public:
	IFCB();
	~IFCB();
	int IFsize;
	int IFindexsize;
	IFCB_Value *IF;
	IFCB_Value** IFCB_Index;
	int read_ifcb(char* filename);
	ifcb_* add(IFCB_Value** PTR, ifcb_* pp, char* prn, char* value, gnsstime time_cur,bool flag,int k);
	double get_g_ifcb(gnsstime t, unsigned int prn);
};

