#pragma once
#include "coord.h"

class GausseIntegr
{
protected:
	static const int POW_GAUSSE = 5;
	static const double R_GAUSSE[POW_GAUSSE];
	static const double W_GAUSSE[POW_GAUSSE];
public:
	virtual double integrateGausse(Coord left, Coord right) = 0;
	virtual double integrFunct(Coord p) = 0;
};

class GausseIntegr3D: GausseIntegr
{
public:
	double integrateGausse(Coord left, Coord right);
	virtual double integrFunct(Coord p) = 0;
};

class GausseIntegr2D : GausseIntegr
{
public:
	double integrateGausse(Coord left, Coord right);
	virtual double integrFunct(Coord p) = 0;
};