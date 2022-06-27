#pragma once
#include "gridbuilder.h"
#include "FinalElemetsBasic.h"

/// <summary>
/// Класс расчёта потоков через грани
/// </summary>
class FlowCulcer
{
	CalculationArea calculationArea;
	PhaseStorage outPhaseStorage;
	double* q;
	int qN;
	double* flows;
	int nFlows;
public:
	FlowCulcer();
	bool isInnerN(int iLocalFace);
	void init(CalculationArea calculationArea, double* q);
	void calcFlows();
	double* getflows();
	void setPressure(double* q);
	void printFlows(std::ofstream &out);

};
