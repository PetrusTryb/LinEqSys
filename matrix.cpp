#include "matrix.h"

Matrix::Matrix()
{
	rows = 0;
	cols = 0;
	data = nullptr;
}

Matrix::Matrix(int rows, int cols)
{
	this->rows = rows;
	this->cols = cols;
	data = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		data[i] = new double[cols];
		for (int j = 0; j < cols; j++)
		{
			data[i][j] = 0;
		}
	}
}

Matrix::Matrix(const Matrix &m)
{
	this->rows = m.rows;
	this->cols = m.cols;
	data = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		data[i] = new double[cols];
		for (int j = 0; j < cols; j++)
		{
			data[i][j] = m.data[i][j];
		}
	}
}

Matrix::~Matrix()
{
	for (int i = 0; i < rows; i++)
	{
		delete[] data[i];
	}
	delete[] data;
}

void Matrix::fillDiagonal(int diagonal, double value)
{
	int startRow = 0;
	int startCol = 0;
	if (diagonal > 0)
		startCol = diagonal;
	else
		startRow = -diagonal;

	for (int row = startRow, col = startCol; row < rows && col < cols; row++, col++)
	{
		data[row][col] = value;
	}
}

void Matrix::print()
{
	for (int i = 0; i < rows; i++)
	{
		printf("Row %d: ", i);
		for (int j = 0; j < cols; j++)
		{
			printf("%lf ", data[i][j]);
		}
		printf("\n");
	}
}

int Matrix::getRows()
{
	return rows;
}

int Matrix::getCols()
{
	return cols;
}

double Matrix::getElement(int i, int j)
{
	return data[i][j];
}

void Matrix::setElement(int i, int j, double value)
{
	data[i][j] = value;
}

Matrix &Matrix::operator=(const Matrix &m)
{
	if (this == &m)
		return *this;

	for (int i = 0; i < rows; i++)
	{
		delete[] data[i];
	}
	delete[] data;

	this->rows = m.rows;
	this->cols = m.cols;
	data = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		data[i] = new double[cols];
		for (int j = 0; j < cols; j++)
		{
			data[i][j] = m.data[i][j];
		}
	}
	return *this;
}

Matrix Matrix::operator+(const Matrix &m)
{
	ASSERT(rows == m.rows && cols == m.cols, "Matrix dimensions must match\n");
	Matrix result(rows, cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			result.data[i][j] = data[i][j] + m.data[i][j];
		}
	}
	return result;
}

Matrix Matrix::operator-(const Matrix &m)
{
	ASSERT(rows == m.rows && cols == m.cols, "Matrix dimensions must match\n");
	Matrix result(rows, cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			result.data[i][j] = data[i][j] - m.data[i][j];
		}
	}
	return result;
}

Matrix Matrix::operator*(const Matrix &m)
{
	ASSERT(cols == m.rows, "Matrix dimensions must match\n");
	Matrix result(rows, m.cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < m.cols; j++)
		{
			for (int k = 0; k < cols; k++)
			{
				result.data[i][j] += data[i][k] * m.data[k][j];
			}
		}
	}
	return result;
}

LU Matrix::decomposeLU()
{
	ASSERT(rows == cols, "Matrix must be square\n");
	int N = rows;
	LU LU_arr;
	LU_arr.L = new Matrix(N,N);
	LU_arr.U = new Matrix(N,N);
	
	LU_arr.L->fillDiagonal(0,1);

	for (int i = 0; i < N; i++)
	{
		for (int j = i; j < N; j++)
		{
			double sum = 0;
			for (int k = 0; k < i; k++)
			{
				sum += LU_arr.L->getElement(i,k) * LU_arr.U->getElement(k,j);
			}
			LU_arr.U->setElement(i,j,data[i][j]-sum);
		}
		for (int j = i; j < N; j++)
		{
			double sum = 0;
			for (int k = 0; k < i; k++)
			{
				sum += LU_arr.L->getElement(j,k) * LU_arr.U->getElement(k,i);
			}
			LU_arr.L->setElement(j,i,(data[j][i]-sum)/LU_arr.U->getElement(i,i));
		}
	}
	return LU_arr;
}
