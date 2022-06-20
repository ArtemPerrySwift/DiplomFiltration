#pragma once
#include "sparsematrix.h"

using namespace matrix;
typedef double real;
namespace slae
{
	template<class MatrixType, typename = std::enable_if_t<std::is_base_of<SparseMatrix, MatrixType>::value>>
	struct SLAE
	{
		MatrixType A;
		double* b;
		void init(MatrixType A, double* b)
		{
			this->A = A;
			this->b = b;
		}

		void init(MatrixType A)
		{
			b = new double[A.n];
			init(A, b);
		}

		void setOneVariableSolve(int iVar, double varMean);
	};

	/*
	struct SLAE
	{
		SparseMatrixAsym A;
		SparseMatrixAsym LU;

		double* b;

		int count_LOS(real* x, int maxiter, real eps);
		int count_LOS_Simple(real* x, int maxiter, real eps);
		void init(SparseMatrixAsym A);
		void init(SparseMatrixAsym A, double* b);

		void setOneVariableSolve(int iVar, double varMean);
	private:

		real* r;;
		real* z; 
		real* p; 
		real* f; 
		real* buf_v; 
		real* buf_v1;

		int calc_Lx(SparseMatrixAsym& LU, real* x, real* f, bool is_diag_1);

		int calc_Ux(SparseMatrixAsym& LU, real* x, real* f, bool is_diag_1);	
	};
	*/
	int solveSLAU3(double A[3][3], double b[3], double x[3]);
	
}
