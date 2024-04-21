#include "matrix.h"
#include <math.h>
#include <chrono>
#include <vector>
using namespace std;
using namespace std::chrono;

vector<double> solve_Jacobi(Matrix& A, Matrix& b, double err_norm, double failure_threshold = 1e5)
{
	ASSERT(A.getRows() == A.getCols(), "Matrix A must be square\n");
	ASSERT(A.getRows() == b.getRows(), "Matrix A and b must have the same number of rows\n");
	ASSERT(b.getCols() == 1, "Matrix b must have only one column\n");
	int N = A.getRows();
	Matrix x(N);
	vector<double> error(1);
	error[0] = 1;
	int iter = 0;

	auto start = high_resolution_clock::now();
	while(error[iter]>=err_norm){
		Matrix x_new(N);
		for (int i = 0; i < N; i++){
			double sum = 0;
			for (int j = 0; j < N; j++){
				if (i != j){
					sum += A.getElement(i, j) * x.getElement(j);
				}
			}
			x_new.setElement(i, 0, (b.getElement(i) - sum) / A.getElement(i, i));
		}
		iter++;
		x = x_new;
		Matrix error_vec = (A*x) - b;
		error.push_back(0);
		for (int i = 0; i < N; i++){
			error[iter] += error_vec.getElement(i) * error_vec.getElement(i);
		}
		error[iter] = sqrt(error[iter]);
		if(error[iter]>failure_threshold){
			printf(":\"Jacobi\", \"%d\", \"-\", \"%e\",\n", iter, error[iter]);
			return error;
		}
	}
	auto end = high_resolution_clock::now();
	duration<double> elapsed = end - start;
	printf(":\"Jacobi\", \"%d\", \"%.2f\", \"%e\",\n", iter, elapsed.count(), error[iter]);
	error[0] = elapsed.count();
	return error;
}

vector<double> solve_GaussSeidel(Matrix& A, Matrix& b, double err_norm, double failure_threshold = 1e5){
	ASSERT(A.getRows() == A.getCols(), "Matrix A must be square\n");
	ASSERT(A.getRows() == b.getRows(), "Matrix A and b must have the same number of rows\n");
	ASSERT(b.getCols() == 1, "Matrix b must have only one column\n");

	int N = A.getRows();
	Matrix x(N);
	vector<double> error(1);
	error[0] = 1;
	int iter = 0;

	auto start = high_resolution_clock::now();
	while(error[iter]>=err_norm){
		Matrix x_new(N);
		for(int i = 0; i<N; i++){
			double sum1 = 0;
			double sum2 = 0;
			for(int j = 0; j<i; j++){
				sum1 += A.getElement(i,j)*x_new.getElement(j);
			}
			for(int j = i+1; j<N; j++){
				sum2 += A.getElement(i,j)*x.getElement(j);
			}
			x_new.setElement(i,0,(b.getElement(i)-sum1-sum2)/A.getElement(i,i));
		}
		iter++;
		x = x_new;
		Matrix error_vec = (A*x) - b;
		error.push_back(0);
		for (int i = 0; i < N; i++){
			error[iter] += error_vec.getElement(i) * error_vec.getElement(i);
		}
		error[iter] = sqrt(error[iter]);
		if(error[iter]>failure_threshold){
			printf(":\"Gauss-Seidel\", \"%d\", \"-\", \"%e\",\n", iter, error[iter]);
			return error;
		}
	}
	auto end = high_resolution_clock::now();
	duration<double> elapsed = end - start;
	printf(":\"Gauss-Seidel\", \"%d\", \"%.2f\", \"%e\",\n", iter, elapsed.count(), error[iter]);
	error[0] = elapsed.count();
	return error;
}

double solve_LU(Matrix& A, Matrix& b){
	ASSERT(A.getRows() == A.getCols(), "Matrix A must be square\n");
	ASSERT(A.getRows() == b.getRows(), "Matrix A and b must have the same number of rows\n");
	ASSERT(b.getCols() == 1, "Matrix b must have only one column\n");

	int N = A.getRows();
	auto start = high_resolution_clock::now();
	LU lu = A.decomposeLU();

	Matrix y(N);
	Matrix x(N);

	for(int i = 0; i<N; i++){
		double sum = 0;
		for(int j = 0; j<i; j++){
			sum += lu.L->getElement(i,j)*y.getElement(j);
		}
		y.setElement(i,0,(b.getElement(i)-sum)/lu.L->getElement(i,i));
	}

	for(int i = N-1; i>=0; i--){
		double sum = 0;
		for(int j = i+1; j<N; j++){
			sum += lu.U->getElement(i,j)*x.getElement(j);
		}
		x.setElement(i,0,(y.getElement(i)-sum)/lu.U->getElement(i,i));
	}

	Matrix error_vec = (A*x) - b;
	double error = 0;
	for (int i = 0; i < N; i++){
		error += error_vec.getElement(i) * error_vec.getElement(i);
	}
	error = sqrt(error);
	auto end = high_resolution_clock::now();
	duration<double> elapsed = end - start;
	printf(":\"Faktoryzacja LU\", \"nie dotyczy\", \"%.2f\", \"%e\",\n", elapsed.count(), error);
	delete lu.L;
	delete lu.U;
	return elapsed.count();
}

int main()
{
	const int index[] = {1,9,3,5,5,7};
	//const int index[] = {1,9,3,3,2,8};
	//const int index[] = {1,9,3,0,4,4};
	//Zadanie A
	int a1 = 5+index[3];
	int a2 = -1;
	int a3 = -1;
	int N = 900+index[4]*10+index[5];
	Matrix A(N,N);
	A.fillDiagonal(0, a1);
	A.fillDiagonal(1, a2);
	A.fillDiagonal(-1, a2);
	A.fillDiagonal(2, a3);
	A.fillDiagonal(-2, a3);
	Matrix b(N);
	for (int i = 0; i < N; i++)
	{
		b.setElement(i,0,sin((i+1)*(index[2]+1)));
		if(i<5||i>N-6)
			printf("%d: %.3f;\n",i, b.getElement(i));
		if(i==6)
			printf("dots.v;\n");
	}
	//Zadanie B
	auto jacobi1_err = solve_Jacobi(A, b, 1e-9);
	auto gs1_err = solve_GaussSeidel(A, b, 1e-9);
	FILE* output_file = fopen("C:\\Users\\Ptryb\\LinEqSys\\B_jacobi.txt", "w");
	for(int i = 1; i<jacobi1_err.size(); i++){
		fprintf(output_file, "%e\n", jacobi1_err[i]);
	}
	fclose(output_file);
	output_file = fopen("C:\\Users\\Ptryb\\LinEqSys\\B_gauss_seidel.txt", "w");
	for(int i = 1; i<gs1_err.size(); i++){
		fprintf(output_file, "%e\n", gs1_err[i]);
	}
	fclose(output_file);
	
	solve_LU(A, b);
	//Zadanie C
	a1 = 3;
	Matrix A2(N,N);
	A2.fillDiagonal(0, a1);
	A2.fillDiagonal(1, a2);
	A2.fillDiagonal(-1, a2);
	A2.fillDiagonal(2, a3);
	A2.fillDiagonal(-2, a3);
	auto jacobi2_err = solve_Jacobi(A2, b, 1e-9);
	auto gs2_err = solve_GaussSeidel(A2, b, 1e-9);
	output_file = fopen("C:\\Users\\Ptryb\\LinEqSys\\C_jacobi.txt", "w");
	for(int i = 1; i<jacobi2_err.size(); i++){
		fprintf(output_file, "%e\n", jacobi2_err[i]);
	}
	fclose(output_file);
	output_file = fopen("C:\\Users\\Ptryb\\LinEqSys\\C_gauss_seidel.txt", "w");
	for(int i = 1; i<gs2_err.size(); i++){
		fprintf(output_file, "%e\n", gs2_err[i]);
	}
	fclose(output_file);
	//Zadanie D
	solve_LU(A2, b);
	//Zadanie E
	auto output_j = fopen("C:\\Users\\Ptryb\\LinEqSys\\E_jacobi.txt", "w");
	auto output_gs = fopen("C:\\Users\\Ptryb\\LinEqSys\\E_gauss_seidel.txt", "w");
	auto output_lu = fopen("C:\\Users\\Ptryb\\LinEqSys\\E_lu.txt", "w");
	int sizes[] = {100,500,1000,2000,3000,5000};
	a1 = 5+index[3];
	for(int i = 0; i<6; i++){
		a1 = 5+index[3];
		N = sizes[i];
		Matrix A3(N,N);
		A3.fillDiagonal(0, a1);
		A3.fillDiagonal(1, a2);
		A3.fillDiagonal(-1, a2);
		A3.fillDiagonal(2, a3);
		A3.fillDiagonal(-2, a3);
		Matrix b3(N);
		for (int i = 0; i < N; i++)
		{
			b3.setElement(i,0,sin((i+1)*(index[2]+1)));
		}
		double jacobi_time = solve_Jacobi(A3, b3, 1e-9)[0];
		double gs_time = solve_GaussSeidel(A3, b3, 1e-9)[0];
		double lu_time = solve_LU(A3, b3);
		fprintf(output_j, "%e\n", jacobi_time);
		fprintf(output_gs, "%e\n", gs_time);
		fprintf(output_lu, "%e\n", lu_time);
	}
	fclose(output_j);
	fclose(output_gs);
	fclose(output_lu);
	return 0;
}