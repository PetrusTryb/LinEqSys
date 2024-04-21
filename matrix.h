#pragma once
#include <stdio.h>
#include <stdlib.h>
#define ASSERT(cond, msg) if (!(cond)) { printf(msg); exit(1); }

typedef struct LU LU;

class Matrix
{
 public:
	Matrix();
	Matrix(int rows, int cols=1);
	Matrix(const Matrix& m);
	~Matrix();
	void fillDiagonal(int diagonal, double value);
	void print();
	int getRows();
	int getCols();
	double getElement(int i, int j=0);
	void setElement(int i, int j, double value);
	Matrix& operator=(const Matrix& m);
	Matrix operator+(const Matrix& m);
	Matrix operator-(const Matrix& m);
	Matrix operator*(const Matrix& m);
	LU decomposeLU();

 private:
	int rows;
	int cols;
	double** data;
};

struct LU
{
	Matrix* L;
	Matrix* U;
};