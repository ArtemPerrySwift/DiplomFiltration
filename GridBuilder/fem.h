#pragma once
#include "slae.h"
//#include "mesh.h"
#include "FinalElemetsBasic.h"
#include "GlobalMatAssembler.h"
#include "finitelemculcer.h"
#define NON_LINEAR false

//using namespace meshspace;
using namespace slae;

class FEM
{
	SLAE<SparseMatrixSym> slae;
	//Mesh mesh;
	//NewtAssembler newtAssembler;
	CalculationArea calculationArea;
	SimpleAssembler simpleAssembler;


	double* descr;
	double* bForCheck;
	void addFirstConditions();
	void addSecondConditions();
	//void addEmptyConditions();
	void LinearTask();
	void buildMatrixPortrait();
#if NON_LINEAR
	void NonLinearTask();
#endif

	/// <summary>
	/// Расчёт невязки
	/// </summary>
	double countDescr();
	double countChange();
	double* qPrev;
	FinitElemCulcer finitElemCulcer;
public:
	double* q;
	unsigned int nBufLookingArea;
	unsigned int nBufLookingCenter;
	//double solution(double x, double y, double& A, double& Bx, double& By);
	void init(CalculationArea& calculationArea);
	void printSolution(std::ofstream& out);
	double* getSolutWeights();
	void getSolutWeights(double*& qRes);
	double countSolution(Coord p);
};

