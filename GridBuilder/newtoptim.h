#pragma once
#include <cfloat>
class NewtOptim
{
public:
	double calcMin(double initMean);
	unsigned int maxiter;
protected:
	double eps;
	virtual double functMin(double mean) = 0; // �������������� �������
	virtual double diffFuntcMin(double mean); // ����������� �������������� �������
private:
	double minPossibleMeanEps;
	double bettaCur, bettaPrev;
	double meanCur, meanPrev;
	bool fl;
};
