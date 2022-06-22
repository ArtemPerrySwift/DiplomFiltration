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
	BorderFacesStore faces2CondStore;
	slae::SLAE<SparseMatrixSym> slae;
	LOS los;

	double* betta;
	double alpha;
	double* dQ;
	double* diBuf;
	double* bNew;
	int nFaces;
	double epsBalance;
	int nElems;

	void buildPortrait();
	void fillSlae();
	double functMin(double mean) override;
	bool changeMemory(unsigned int nFaces, unsigned int nElems);
	virtual void deleteMemory();
	virtual bool allocateMemory(unsigned int n);
	double findAlpha();
	double reculcBetta();
	double calcDisbalanse();
	void add2Cond();
	void refreshFlowSigns();
public:
	Balance();

	void init(CalculationArea calculationArea);
	void balanceFlows();
	bool setEpsBalanse(double eps);
	//bool checkBalance();
	
	
};