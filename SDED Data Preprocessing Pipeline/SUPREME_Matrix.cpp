#include <iostream>
#include "SUPREME_Matrix.h"
#pragma warning(disable:4996)

MatrixT::MatrixT()
{
	m_row = m_col = 0;
	m_mat.resize(m_row, m_col);
	m_mat.setZero();
}

MatrixT::MatrixT(const MatrixT& mat)
{
	m_row = mat.GetRow();
	m_col = mat.GetCol();

	m_mat = mat.m_mat;
}

MatrixT::MatrixT(unsigned int r)
{
	m_row = m_col = r;
	m_mat.resize(m_row, m_col);
	m_mat.setZero();
}

MatrixT::MatrixT(unsigned int row, unsigned int col)
{
	m_row = row;
	m_col = col;
	m_mat.resize(m_row, m_col);
	m_mat.setZero();
}

MatrixT::~MatrixT()
{
	m_row = m_col = 0;
	m_mat.resize(m_row, m_col);
	m_mat.setZero();
}

int MatrixT::Initialize(unsigned int row, unsigned int col)
{
	m_row = row; m_col = col;
	m_mat.resize(m_row, m_col);
	m_mat.setZero();

	return 1;
}

int MatrixT::Resize(unsigned int row, unsigned int col)
{
	m_row = row; m_col = col;
	m_mat.resize(m_row, m_col);
	m_mat.setZero();

	return 1;
}

double MatrixT::GetValue(unsigned int row, unsigned int col)
{
	return (m_mat(row - 1, col - 1));
}

int MatrixT::SetValue(unsigned int row, unsigned int col, double value)
{
	m_mat(row - 1, col - 1) = value;
	return 1;
}

double MatrixT::GetMax(int& index)
{
	double maxvalue = m_mat(0, 0);
	index = 0;
	for (unsigned int i = 0; i < m_row; ++i)
	{
		for (unsigned int j = 0; j < m_col; ++j)
		{
			if (m_mat(i, j) > maxvalue)
			{
				maxvalue = m_mat(i, j);
				index = m_col * i + j;
			}
		}
	}

	return maxvalue;
}

double MatrixT::GetMin(int& index)
{
	double minvalue = m_mat(0, 0);
	index = 0;
	for (unsigned int i = 0; i < m_row; ++i)
	{
		for (unsigned int j = 0; j < m_col; ++j)
		{
			if (m_mat(i, j) < minvalue)
			{
				minvalue = m_mat(i, j);
				index = m_col * i + j;
			}
		}
	}

	return minvalue;
}

double MatrixT::GetMaxAbs(int& index)
{
	double maxvalue = abs(m_mat(0, 0));
	index = 0;
	for (unsigned int i = 0; i < m_row; ++i)
	{
		for (unsigned int j = 0; j < m_col; ++j)
		{
			if (abs(m_mat(i, j)) > maxvalue)
			{
				maxvalue = abs(m_mat(i, j));
				index = m_col * i + j;
			}
		}
	}

	return maxvalue;
}

double MatrixT::GetMinAbs(int& index)
{
	double minvalue = abs(m_mat(0, 0));
	index = 0;
	for (unsigned int i = 0; i < m_row; ++i)
	{
		for (unsigned int j = 0; j < m_col; ++j)
		{
			if (abs(m_mat(i, j)) < minvalue)
			{
				minvalue = abs(m_mat(i, j));
				index = m_col * i + j;
			}
		}
	}

	return minvalue;
}

int MatrixT::GetRank()
{
	int rank = 0;

	return rank;
}

MatrixT MatrixT::Cholesky()
{
	MatrixT chol(m_row, m_col);

	chol.m_mat = m_mat.llt().matrixL();

	return chol;
}

MatrixT MatrixT::Trans()
{
	MatrixT trans_tmp(m_col, m_row);
	trans_tmp.m_mat = m_mat.transpose();

	return trans_tmp;
}

MatrixT MatrixT::Inv()
{
	MatrixT inv_tmp(m_row, m_col);
	inv_tmp.m_mat = m_mat.inverse();

	return inv_tmp;
}

MatrixT MatrixT::InvCholesky()
{
	MatrixT invchol_tmp(m_row, m_col);

	/*LLT<MatrixXd> lltOfA(m_mat);
	MatrixXd R = lltOfA.matrixL();*/

	//MatrixXd R = m_mat.llt().matrixL();

	invchol_tmp.m_mat = R.inverse().transpose() * R.inverse();

	return invchol_tmp;
}

MatrixT MatrixT::Eye(int A)
{
	m_row = A;
	m_col = A;

	m_mat.resize(m_row, m_col);
	m_mat.setZero();
	for (int i = 0; i < A; i++)
	{
		for (int j = 0; j < A; j++)
		{
			if (i == j)
			{
				m_mat(i, j) = 1.0;
			}
		}
	}
	return *this;
}

double& MatrixT:: operator () (unsigned int row, unsigned int col)
{
	if (row > m_row || col > m_col)
	{
		printf("Error: Matrix index out of range.\n");
	}

	return m_mat(row - 1, col - 1);
}

// subscript operator to get individual elements ()
double MatrixT::operator () (unsigned int row, unsigned int col)const
{
	if (row > m_row || col > m_col)
		printf("Error: Matrix index out of range.\n");

	return m_mat(row - 1, col - 1);
}

MatrixT MatrixT::operator + (const MatrixT& mat)
{
	MatrixT add_tmp(m_row, m_col);

	add_tmp.m_mat = m_mat + mat.m_mat;

	return add_tmp;
}

MatrixT MatrixT::operator - (const MatrixT& mat)
{
	MatrixT minus_tmp(m_row, m_col);

	minus_tmp.m_mat = m_mat - mat.m_mat;

	return minus_tmp;
}

MatrixT MatrixT::operator * (const MatrixT& mat)
{
	MatrixT multi_tmp(m_row, mat.m_col);

	if (m_col != mat.m_row)
	{
		printf("Error: Column and row is not equal in matrix multiply.\n");

#ifdef _WIN32
		system("pause");
		exit(0);
#endif
	}

	multi_tmp.m_mat = m_mat * mat.m_mat;

	return multi_tmp;
}

MatrixT MatrixT::operator * (const double elem)
{
	MatrixT multi_tmp(m_row, m_col);

	multi_tmp.m_mat = m_mat * elem;

	return multi_tmp;
}

MatrixT MatrixT::operator / (const double elem)
{
	MatrixT multi_tmp(m_row, m_col);

	multi_tmp.m_mat = m_mat / elem;

	return multi_tmp;
}


void MatrixT::operator = (const MatrixT& mat)
{
	m_row = mat.m_row;
	m_col = mat.m_col;

	m_mat = mat.m_mat;
}

void MatrixT::Clear()
{
	m_row = m_col = 0;
	m_mat.resize(0, 0);
}

/// Output matrix to screen
void MatrixT::display()
{
	cout << endl << m_mat << endl << endl;
}

/// Output matrix to file
void MatrixT::Output_to_File(const char* path, unsigned int flag)
{
	FILE* fp = NULL;
	if (flag)	fp = fopen(path, "a+");
	else		fp = fopen(path, "w");

	if (fp)
	{
		fprintf(fp, "\n");
		for (unsigned int i = 0; i < m_row; ++i)
		{
			for (unsigned int j = 0; j < m_col; ++j)
				fprintf(fp, " %15.9f ", m_mat(i, j));
			fprintf(fp, "\n");
		}
	}

	fprintf(fp, "\n");
	if (fp) fclose(fp);
}
void MatrixT::Output_to_File_no(const char* path, unsigned int flag)
{
	FILE* fp = NULL;
	if (flag)	fp = fopen(path, "a+");
	else		fp = fopen(path, "w");

	if (fp)
	{

		for (unsigned int i = 0; i < m_row; ++i)
		{
			for (unsigned int j = 0; j < m_col; ++j)
				fprintf(fp, " %15.9f ", m_mat(i, j));
			fprintf(fp, "\n");
		}
	}

	if (fp) fclose(fp);
}