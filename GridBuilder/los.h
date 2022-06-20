#pragma once
#include "sparsematrix.h"
#include "array.h"
#include <iostream>
#include <iomanip>

using namespace matrix;
class LOS
{
public:
	virtual int solve(Matrix& A, double* b, double* x, int maxiter, double eps);
	LOS();
protected:
	int n;
	double* r;
	double* z;
	double* p;
	double* buf_v;

	bool changeMemory(unsigned int n);
	virtual void deleteMemory();
	virtual bool allocateMemory(unsigned int n);
	//double* buf_v1;
};

class LOS_precond : public LOS
{
public:
	DecompSparseMatrixLDLt LLt;
	int solve(SparseMatrixSym& A, double* b, double* x, int maxiter, double eps);
private:
	double* buf_v1;
	void deleteMemory() override;
	bool allocateMemory(unsigned int n) override;
};
