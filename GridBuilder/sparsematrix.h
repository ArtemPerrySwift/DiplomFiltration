#pragma once
//#include "mesh.h"
#include "gridbuilder.h"

namespace matrix
{
	struct Matrix
	{
		double* operator * (double* v);
		virtual void mult(double* v, double* res) = 0;
		virtual void fillMatrix(double mean) = 0;
		int n;
	};
	struct SparseMatrix : public Matrix
	{
		int* ig;
		int* jg;
		double* di;	
		SparseMatrix();
		//SparseMatrixPortrait(SparseMatrixPortrait& portarait);	
		void buildPortrait(CoordStorage coordStorage, FinitElementStore finitElementStore);
		void mult(double* v, double* res) override;
		virtual void allocateMemoryForElems();
	protected:
		//void copy();
		void copyPortrait(SparseMatrix& M);
		bool copyElems(SparseMatrix& M);
		void changePortrait(SparseMatrix& M);
		virtual void multMatElemAndVecElem(int i, int j, int ij, double* v, double*res) = 0;
		
		virtual void freeMemory();
		//virtual void printInFullFormat() = 0;
	};


	struct SparseMatrixSym : SparseMatrix
	{
		double* gg;

		SparseMatrixSym();

		void copy(SparseMatrixSym& M);

		//operator SparseMatrixAsym() const;

		void fillMatrix(double mean) override;
		void allocateMemoryForElems() override;
	private:
		bool copyElems(SparseMatrixSym& M);
		void multMatElemAndVecElem(int i, int j, int ij, double* v, double* res) override;

		void freeMemory() override;
	};

	struct SparseMatrixAsym :  SparseMatrix
	{
		double* ggl;
		double* ggu;

		SparseMatrixAsym();

		void copy(SparseMatrixAsym& M);
		void copy(SparseMatrixSym& M);

		int decomp_mat_LU(SparseMatrixAsym& LU);

		void printFullMatrix();

		void fillMatrix(double mean) override;
		void allocateMemoryForElems() override;
	private:
		bool copyElems(SparseMatrixAsym& M);
		bool copyElems(SparseMatrixSym& M);
		void multMatElemAndVecElem(int i, int j, int ij, double* v, double* res) override;
		
		void freeMemory() override;
	};


	//template<class MatrixType, typename = std::enable_if_t<std::is_base_of<Matrix, MatrixType>::value>>
	//struct DecompMatrix
	//{
		//virtual void decompMatrix(MatrixType M) = 0;
		//virtual void solveLy(double* f, double* y) = 0;
		//virtual void solveUx(double* y, double* x) = 0;
	//};

	//template<class MatrixType, typename = std::enable_if_t<std::is_base_of<Matrix, MatrixType>::value>>
	struct DecompSparseMatrixLDLt : SparseMatrixSym
	{
	public:
		void decompMatrix(SparseMatrixSym&M);
		virtual void solveLy(double* f, double* y);
		virtual void solveUx(double* y, double* x);
		void solveDz(double* f, double* z);
	};
}

