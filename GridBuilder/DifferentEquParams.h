#pragma once
#include "phase.h"

class DifferentEquParams
{
public:
	static double f(double x, double y, double z);

	static double u1(double x, double y, double z, int functNum);

	static double du_dn(double x, double y, double z, int functNum);

	static double lambda(double K, Phase* phase, int nPhases);

	static double km(double str);
};

