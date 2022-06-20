#include "gaussintegr.h"
#include <iostream>

const double GausseIntegr::R_GAUSSE[POW_GAUSSE] = { -0.90617984593866399279763,-0.5384693101056830910363, 0, 0.5384693101056830910363,0.90617984593866399279763 };
const double GausseIntegr::W_GAUSSE[POW_GAUSSE] = { 0.23692688505618908751426, 0.4786286704993664680412, 0.5688888888888888888889 , 0.4786286704993664680412, 0.23692688505618908751426 };

double GausseIntegr3D::integrateGausse(Coord left, Coord right)
{
	Coord avr = (left + right) / 2.0;
	Coord hCoordHalf = (right - left) / 2.0;
	Coord r_gausse;
	double sum = 0;

	int l, m, k;
	double mean;
	//std::cout << "\t Itegral is started" << std::endl;
	for (l = 0; l < POW_GAUSSE; l++)
	{
		for (m = 0; m < POW_GAUSSE; m++)
		{
			for (k = 0; k < POW_GAUSSE; k++)
			{
				//std::cout << "\t Itegral is counting 1" << std::endl;
				r_gausse = { R_GAUSSE[l], R_GAUSSE[m],  R_GAUSSE[k] };
				//std::cout << "\t Itegral is counting 2" << std::endl;
				mean = integrFunct(avr + r_gausse * hCoordHalf);
				//std::cout << "\t Itegral is counting 3" << std::endl;
				sum += W_GAUSSE[l] * W_GAUSSE[m] * W_GAUSSE[k] * mean;
				//std::cout << "\t Itegral is counting 4" << std::endl;
			}

		}
	}
	//std::cout << "\t Itegral is ending" << std::endl;
	return hCoordHalf.x * hCoordHalf.y * hCoordHalf.z * sum;
}

double GausseIntegr2D::integrateGausse(Coord left, Coord right)
{
	Coord avr = (left + right) / 2.0;
	Coord hCoordHalf = (right - left) / 2.0;
	Coord r_gausse;
	double sum = 0;

	int l, m, k;
	for (l = 0; l < POW_GAUSSE; l++)
	{
		for (m = 0; m < POW_GAUSSE; m++)
		{
			r_gausse = { R_GAUSSE[l], R_GAUSSE[m], R_GAUSSE[l] };
			sum += W_GAUSSE[l] * W_GAUSSE[m] * integrFunct(avr + r_gausse * hCoordHalf);
		}
	}
	return hCoordHalf.x * hCoordHalf.y * sum;
}