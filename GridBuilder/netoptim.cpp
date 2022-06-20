#include "newtoptim.h"
#include <cmath>

double NewtOptim::calcMin(double initBetta)
{
	bettaPrev = initBetta;
	meanPrev = functMin(bettaPrev);
	if (abs(meanPrev) < eps) return bettaPrev;

	bettaCur = bettaCur == 0 ? 1 : bettaPrev * 1.2;
	meanCur = functMin(bettaCur);
	double diff;
	for (int i = 0; i < maxiter && abs(meanCur) < eps; i++)
	{
		diff = diffFuntcMin(bettaCur);
		if (diff == 0)
			return bettaCur;
		bettaPrev = bettaCur;
		bettaCur -= meanCur / diff;
		meanPrev = meanCur;
		meanCur = functMin(bettaCur);
	}

	return bettaCur;
}

double NewtOptim::diffFuntcMin(double mean)
{
	if (abs((meanCur - meanPrev) / meanCur) < eps)
	{
		bool bettaCurChanged = false;

		if (minPossibleMeanEps > abs(bettaCur))
		{
			bettaCur = minPossibleMeanEps;
			bettaCurChanged = true;
		}
		while (abs((meanCur - meanPrev) / meanCur) < eps && abs((bettaCur - bettaPrev)/bettaCur) > eps)
		{
			bettaPrev = 1.0 / 2.0 * (bettaCur - bettaPrev);
			meanPrev = functMin(bettaPrev);
		}

		if (abs((bettaCur - bettaPrev) / bettaCur) <= eps)
		{
			if (bettaCurChanged) bettaCur = 0;
			return 0;
		}
		if (bettaCurChanged) bettaCur = 0;
	}
	return (meanCur - meanPrev) / (bettaCur - bettaPrev);
}