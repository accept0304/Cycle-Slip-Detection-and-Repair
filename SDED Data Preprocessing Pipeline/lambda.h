#pragma once
#include <iostream>
//#include"SUPREME_CommonFunction.h"
#include <math.h>




#define ROUND(x)    (int)floor((x)+0.5)
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)
#define SGN(x)      ((x)<=0.0?-1.0:1.0)

double* zeros(int n, int m);
double* mat(int n, int m);
int lambda(

	int n,
	int m,
	const double* a,
	const double* Q,
	double* F,
	double* s,
	double Pf,
	bool& pass);