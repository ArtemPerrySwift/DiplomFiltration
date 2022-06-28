#pragma once
#include "gridbuilder.h"
#include <vector>
#include <set>

struct PhaseOut
{
	int iFinElem;
	int iPhase;

	/*bool operator < (const PhaseOut op)
	{
		if (iFinElem == op.iFinElem)
			return iPhase < op.iPhase;
		return iFinElem < op.iFinElem;
	}
	*/
};


bool operator <(const PhaseOut& left, const PhaseOut& right);

class SatCulcer
{
private:
	static const double S_MIN;
	static const double S_MAX;
	double dtMax;
	CalculationArea calculationArea;
	double* flows;
	double* phaseSatGen;
	double genPoreVol;
	std::vector<PhaseOut> phaseSatSmin;
	std::vector<PhaseOut> phaseSatSmax;
	std::set<PhaseOut> phaseOut;

	double choose_dt();
	void pushOutPhases();
	double calcOutQSum(FinitElement finitElement);
	//void calcOutVol();
	double calcPhaseCoeff(Phase phase, double coeffSum);
	double calcCoeffSum(PhaseStorage phaseStorage);
	double calcAllVol(double dt);
	void calcNewSat();
	void calcBorderVol(BorderFacesPhases borderFacesStore, double dt);
	void calcGenSat();
	void calcGenPoreVol();
	void calcGenVol();
public:
	void init(CalculationArea calculationArea, double dtMax);
	void printSat(std::ofstream &out);
	double reculcSat();
	double getGenPhaseSat(int iPhase);
};
