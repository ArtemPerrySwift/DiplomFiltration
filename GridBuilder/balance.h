#pragma once
#include "gridbuilder.h"
#include "slae.h"
#include "los.h"
#include "newtoptim.h"

class Balance : NewtOptim
{
	FaceStore faceStore;
	FlowStore flowStore;
	FinitElementStore finitElementStore;
	slae::SLAE<SparseMatrixSym> slae;
	LOS_precond los_prec;

	double* betta;
	double alpha;
	double* dQ;
	double* diBuf;
	double* bNew;
	int nFace;
	double epsBalance;

	void buildPortrait();
	void fillSlae();
	double functMin(double mean) override;
	bool changeMemory(unsigned int n);
	virtual void deleteMemory();
	virtual bool allocateMemory(unsigned int n);
	double findAlpha();
	double reculcBetta();

public:
	Balance();

	void init(CalculationArea calculationArea);
	void balanceFlows();
	bool setEpsBalanse(double eps);
	//bool checkBalance();
	
	
};