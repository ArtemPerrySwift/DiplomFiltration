#pragma once
#include <cfloat>
class NewtOptim
{
public:
	double calcMin(double initMean);
	unsigned int maxiter;
protected:
	double eps;
	virtual double functMin(double mean) = 0; // Минимизируемая функция
	virtual double diffFuntcMin(double mean); // Производная минимизируемой функции
private:
	double minPossibleMeanEps;
	double bettaCur, bettaPrev;
	double meanCur, meanPrev;
	bool fl;
};
