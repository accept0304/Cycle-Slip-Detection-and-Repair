#include"SUPREME_Ifcb.h"
#include<iostream>
IFCB::IFCB()
{

}
IFCB::~IFCB()
{
	IFCB_Value* ptr;
	ptr = IF; //ptr1 = time;
	while (ptr != NULL)
	{
		IFCB_Value* p = ptr->m_Next;
		free(ptr);
		ptr = p;
	}

}


int IFCB::read_ifcb(char* filename)
{
	char strline[100] = ""; int len = 100; ifcb_* p = NULL; 
	gnsstime time_pre;
	gnsstime time_cur;
	int year = 0, month = 0, day = 0, hour = 0, min = 0, sec = 0;
	char Year[5] = ""; char Month[3] = ""; char Day[3] = ""; char Hour[3] = "";
	char Min[3] = ""; char Sec[3] = "";
	char prn[4]=""; char value[8] = "";
	bool fuhao = true; 
	FILE* fp = fopen(filename, "r");
	vector<unsigned int>Prn_; vector<double>value_;
	if (fp == NULL)
	{
		printf("Error:IFCB file can not be opened.\n");
		return 0;
	}
	
	int k = 0,j = 0;
	p = new ifcb_;
	memset((ifcb_*)p, 0.0, sizeof(ifcb_));
	IFCB_Value* PTR=NULL ;
	while (is_eof(fp) == 0)
	{
		time_pre = time_cur;
		GetFileLine(fp, strline, len);
		GetSubStr(strline, len, prn, 4, 3, 3);
		GetSubStr(strline, len, Year, 5, 8, 4);  // 长度  位置  长度
		GetSubStr(strline, len, Month, 3, 14, 2);
		GetSubStr(strline, len, Day, 3, 16, 2);
		GetSubStr(strline, len, Hour, 3, 19, 2);
		GetSubStr(strline, len, Min, 3, 22, 2);
		GetSubStr(strline, len, Sec, 3, 25, 2);
		GetSubStr(strline, len, value, 8, 43, 7);
		
		time_cur.m_Year = str_to_i(Year); time_cur.m_Month = str_to_i(Month);
		time_cur.m_Day = str_to_i(Day); time_cur.m_Hour = str_to_i(Hour);
		time_cur.m_Min = str_to_i(Min); time_cur.m_Sec = str_to_f(Sec);
		time_cur._date2jd();
		time_cur._date2gpst();
		
		if (time_cur == time_pre ||k==0)
		{
			unsigned int prn_1 = GetSatNo(prn, 4, 0);
			p->Prn[k]= prn_1;
			p->ifcb_value[k] = str_to_f(value);
			
			k++;
			continue; 
		}
		if ((time_cur > time_pre) &&(!fuhao))
		{ 
			p=add(&PTR, p, prn, value, time_cur, fuhao,k);
			j++;
			k = 1;
		}
		else 
		{
			p=add(&PTR, p, prn, value, time_cur, fuhao,k);
			fuhao = false;
			j++;
			k = 1;
		}
	}
	PTR->m_Next = NULL;
	IFCB_Value* P_in = IF;
	if (P_in != NULL) 
	{
		IFCB_Index= (IFCB_Value**)(malloc(j * sizeof(IFCB_Value*)));
	}
	int i = 0;
	while (P_in != NULL)
	{
		IFCB_Index[i] = P_in;
		++i;
		P_in = P_in->m_Next;
	}
	IFindexsize = j;
	return 1;
}

ifcb_* IFCB::add(IFCB_Value **PTR, ifcb_ *pp,char *prn,char *value,gnsstime time_cur,bool flag,int k)
{

	IFCB_Value* p1 = (IFCB_Value*)(malloc(sizeof(IFCB_Value)));
	p1->IFCB_ = pp; p1->IFCB_size = k;
	p1->time = time_cur;
	if (!flag) 
	{
		(*PTR)->m_Next = p1;
		*PTR = p1;
	}
	else
	{
		IF = p1;
		*PTR = p1;
	}
	pp = (ifcb_*)(malloc(sizeof(ifcb_)));
	pp->Prn[0]=GetSatNo(prn, 4, 0);
	pp->ifcb_value[0]=(str_to_f(value));
	return pp;
}

double IFCB::get_g_ifcb(gnsstime t, unsigned int prn)
{
	if (IF == NULL || IFCB_Index == NULL) 
	{
		return 0;
	}
	int i = 0, j = IFindexsize - 1;
	int d = abs(j - i);
	while (d > 1)
	{
		unsigned int k = d / 2 + i;
		if (t < IFCB_Index[k]->time)
			j = k;
		else
			i = k;
		d = abs(j - i);
	}
	if (t == IFCB_Index[i]->time) 
	{
		for (int ii = 0; ii < IFCB_Index[i]->IFCB_size; ii++)
		{
			if (prn== IFCB_Index[i]->IFCB_->Prn[ii])
			{
				return  IFCB_Index[i]->IFCB_->ifcb_value[ii];
			}
		}
	}
	else if(t == IFCB_Index[j]->time)
	{
		for (int ii = 0; ii < IFCB_Index[j]->IFCB_size; ii++)
		{
			if (prn == IFCB_Index[j]->IFCB_->Prn[ii])
			{
				return  IFCB_Index[j]->IFCB_->ifcb_value[ii];
			}
		}
	}
	else
	{
		return 0.0;
	}
}

//ifcb_::ifcb_()
//{
//	std::cout << "shuchu" << endl;
//	int satn = 200;
//	Prn = (unsigned int*)(malloc(sizeof(unsigned int) * satn));
//	memset((unsigned int*)Prn, 0, sizeof(unsigned int) * satn);
//	ifcb_value= (double*)(malloc(sizeof(double) * satn));
//	memset((double*)ifcb_value, 0.0, sizeof(double) * satn);
//}
//
//ifcb_::~ifcb_()
//{
//	if (Prn)free(Prn);
//	if (ifcb_value)free(ifcb_value);
//}
