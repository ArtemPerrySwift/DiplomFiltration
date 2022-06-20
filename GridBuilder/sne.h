#pragma once
#include "coord.h"
using namespace ecoord;
class SNE3D
{
public:
	SNE3D();
	bool countSolution(ECoord initAppr, ECoord& ans);
	bool setMaxiter(int maxiter);
	bool setErr(double err);
	bool setBettaEps(double bettaEps);
	bool setEpsDif(double epsDif);
	bool setMinEpsDif(double minEpsDif);
	bool setMinMean(double minMean);
	bool isNumericalDif;
private:
	static const int DIMMEN = 3;
	double JacMatrix[DIMMEN][DIMMEN];
	double F[DIMMEN];
	double ans[DIMMEN];
	int maxiter;
	double err;
	double bettaEps;
	double epsDif; // Шаг чисенного расчёта производной
	double minEpsDif;
	double minMean;

	void fillJacobMatrix(double J[DIMMEN][DIMMEN], ECoord p);
	void countF(double F[DIMMEN], ECoord p);
	ECoord getIterSolutUsingBetta(ECoord p, double& betta, double& normF);

	virtual double equ1(ECoord p) = 0;
	virtual double equ2(ECoord p) = 0;
	virtual double equ3(ECoord p) = 0;

	virtual double dEqu1dksi(ECoord p);
	virtual double dEqu2dksi(ECoord p);
	virtual double dEqu3dksi(ECoord p);

	virtual double dEqu1dnu(ECoord p);
	virtual double dEqu2dnu(ECoord p);
	virtual double dEqu3dnu(ECoord p);

	virtual double dEqu1detta(ECoord p);
	virtual double dEqu2detta(ECoord p);
	virtual double dEqu3detta(ECoord p);
};
