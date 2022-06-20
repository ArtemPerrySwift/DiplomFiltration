#pragma once
#include "slae.h"
//#include "mesh.h"
#include "FinalElemetsBasic.h"
#include "gridbuilder.h"
using namespace slae;
//using namespace meshspace;

class GlobalMatAssembler
{
public:
	void assembleGlobalMatrix(SLAE<SparseMatrixSym>& A, CalculationArea& mesh, double* q);

protected:
	static const int N_KNOTS_1D = FinalElemBase::N_KNOTS_1D;
	static const int N_KNOTS_2D = FinalElemBase::N_KNOTS_2D;
	static const int N_KNOTS_3D = FinalElemBase::N_KNOTS_3D;

	double localG[N_KNOTS_3D][N_KNOTS_3D], localM[N_KNOTS_3D][N_KNOTS_3D], localB[N_KNOTS_3D];

	void addLocalToGlobalAll(SLAE<SparseMatrixSym>& slae, double G_local[N_KNOTS_3D][N_KNOTS_3D], double M_local[N_KNOTS_3D][N_KNOTS_3D], double localB[N_KNOTS_3D], int L[N_KNOTS_3D]);

	virtual void assembleLocal(FinalElemBase finalElemBase) = 0;
};

/*
class NewtAssembler : public GlobalMatAssembler
{
	void assembleLocal(FinalElemBase finalElemBase) override;
};
*/

class SimpleAssembler :public GlobalMatAssembler
{
	void assembleLocal(FinalElemBase finalElemBase) override;
};

