#pragma once
#include "coord.h"
#include "gaussintegr.h"
#include "gridbuilder.h"

class FaceBasic
{
public:
	static const int DIM = 2;
	static const int N_BASIC = 4;
	void init(Face face, CoordStorage coordStorage);
	void init(Coord coords[N_BASIC]);
	double getDetJ(double ksi, double nu);
protected:	
	Coord coordFace[N_BASIC];
	static double phi2D(double ksi, double nu, int i);

	static double dphi_dksi2D(double ksi, double nu, int i);
	static double dphi_dnu2D(double ksi, double nu, int i);

	static int u(int i);
	static int v(int i);

	Coord dcordFunct_dksi2D(double ksi, double nu);
	Coord dcordFunct_dnu2D(double ksi, double nu);

	Coord cordFunct2D(double ksi, double nu);

	//double integrateGausse(Coord left, Coord right);
};

class FaceBasic2Cond: public FaceBasic, public GausseIntegr2D
{
	double LocalC[N_BASIC][N_BASIC];
	double LocalEtta[N_BASIC];

	virtual double integrFunct(Coord p) override;
	int iPhi2d, jPhi2D;
	int buildLocalM2D();
	void fillLocalEtta();
public:
	void calcAddLocalB(double localAddB[N_BASIC]);
};

class FaceBasicFlow: public FaceBasic
{
	static const int DIM_XYZ = 3;
	bool isInnerNorm;
	double lambda;

	double localQ[N_BASIC];
	void countGrad_i2D(double grad_i[2], int i, double ksi, double nu);

	void buildJ2D(double JacMatrix[DIM][2], double ksi, double nu);
public:
	void init(Face face, CoordStorageXYZ coordStorage, double* q, bool isInnerNorm, double lambda);
	double getFlow();
};

