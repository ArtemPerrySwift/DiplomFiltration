#pragma once
#include "gridbuilder.h"
#include "gaussintegr.h"
#include "basicfunct3d.h"

class FinalElemBase: public GausseIntegr3D, public BasicFunct3D
{
public:
	static const int N_KNOTS_1D = 2;
	static const int N_KNOTS_2D = 4;
	static const int N_KNOTS_3D = 8;

	int globalInd[N_KNOTS_3D];

	void setFinalElement(FinitElement finalElement, CoordStorage knotsStorage);
	void buildLocalM(double M[N_KNOTS_3D][N_KNOTS_3D], double gamma);
	void buildLocalG(double G[N_KNOTS_3D][N_KNOTS_3D]);
	void buildLocalB(double B[N_KNOTS_3D]);
	void fillF(double F[N_KNOTS_3D]);

	double calcSqare();
	
	double funcIntegrG(Coord transform);

private:
	static enum INTEG_CHOICE{INTEG_G, INTEG_SQARE};
	int integChoice;
	int iGradPhi, jGradPhi;

	FinitElement finitElement;
	bool isLinear;

	double str1, str2;
	double nu1, nu2;
	double K;

	double functIntegSqare(Coord p);

	virtual double integrFunct(Coord p) override;
};

