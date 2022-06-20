#pragma once
#include "gridbuilder.h"
#include "slae.h"

class Balance
{
	FaceStore faceStore;
	FlowStore flowStore;
	FinitElementStore finitElementStore;
	slae::SLAE<SparseMatrixSym> slae;
	
	double* betta;
	double alpha;
	double* dQ;

	double epsBalance;
	void buildPortrait();
	void fillSlae();

public:
	void init(CalculationArea calculationArea);
	void balanceFlows();
	bool checkBalance();
	void setEpsBalanse();
	void reculcBetta();
	void findAlpha();
	void trying();
};