/* -------------------------------------------------------------------------
* SUPREME_Matrix.h : Matrix in Multi-GNSS Solutions (Using Eigen Matrix Lib)
*
*           Copyright (C) by C. Zhao All rights reserved
*           Contact caszcb@163.com
*
* Create : 2018.12.09
* ------------------------------------------------------------------------- */
#ifndef SUPREME_MATRIX_H_HH
#define SUPREME_MATRIX_H_HH

#ifndef __cplusplus
#error Error: Must use C++ for the type matrix.
#endif

#include <Dense>
using namespace Eigen;
using namespace std;

/// Class Matrix Defination
class MatrixT
{
public:
	MatrixT();
	MatrixT(const MatrixT& m);
	MatrixT(unsigned int r);
	MatrixT(unsigned int row, unsigned int col);
	~MatrixT();

	unsigned int m_row, m_col;

	int GetRow()const { return m_row; }
	int GetCol()const { return m_col; }
	int Initialize(unsigned int row, unsigned int col);
	int Resize(unsigned int row, unsigned int col);

	double GetValue(unsigned int row, unsigned int col);
	int SetValue(unsigned int row, unsigned int col, double value);

	double GetMax(int& index);
	double GetMin(int& index);
	double GetMaxAbs(int& index);
	double GetMinAbs(int& index);
	int GetRank();

	MatrixT Cholesky();

	MatrixT Trans();
	MatrixT Inv();
	MatrixT InvCholesky();

	MatrixT Eye(int A);
	double& operator()(unsigned int row, unsigned int col);
	double operator()(unsigned int row, unsigned int col)const;

	void operator=(const MatrixT& mat);
	MatrixT operator+(const MatrixT& mat);
	MatrixT operator-(const MatrixT& mat);
	MatrixT operator*(const MatrixT& mat);
	MatrixT operator*(const double elem);
	MatrixT operator/(const double elem);

	void Output_to_File_no(const char* path, unsigned int flag);

	friend ostream& operator<<(ostream&, MatrixT&);
	void Clear();
	void display();
	void Output_to_File(const char* path, unsigned int flag);

	friend double Trop_UNB3(double blh[3], double doy, double elev, double& dry_delay, double& wet_delay, int flag);
	MatrixXd m_mat;
private:

};

#endif
