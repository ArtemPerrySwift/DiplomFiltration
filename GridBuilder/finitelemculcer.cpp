#include "finitelemculcer.h"
#include "array.h"
double FinitElemCulcer::culcPhi(ECoord p)
{
	double phi = 0;
	for (int i = 0; i < N; i++)
		phi += phi_i(p, i)*q[i];

	return phi;
}

ECoord FinitElemCulcer::culcGrad(ECoord p)
{
	ECoord grad(0, 0, 0);
	for (int i = 0; i < N; i++)
		grad += countGrad_i(p, i) * q[i];
	return grad;
}

void FinitElemCulcer::init(FinitElement finitElement, CoordStorage coordsStore, double* q)
{
	this->finitElement = finitElement;

	coordsStore.copyElemsByInd(finitElement.ver, coords, VER_NUM);
	arrayspace::copyElemsByInd(finitElement.ver, this->q, q, VER_NUM);
}

FinitElemCulcer::FinitElemCulcer(): SNE3D()
{
	isNumericalDif = false;
}

double FinitElemCulcer::equ1(ECoord p)
{
	return pointSolut.x - coordFunct_x(p);
}
double FinitElemCulcer::equ2(ECoord p)
{
	return pointSolut.y - coordFunct_y(p);
}
double FinitElemCulcer::equ3(ECoord p)
{
	return pointSolut.z - coordFunct_z(p);
}

double FinitElemCulcer::dEqu1dksi(ECoord p)
{
	Coord p1;
	p1 = dcordFunct_dksi(p);
	return -p1.x;
}

double FinitElemCulcer::dEqu2dksi(ECoord p)
{
	Coord p1;
	p1 = dcordFunct_dksi(p);
	return -p1.y;
}

double FinitElemCulcer::dEqu3dksi(ECoord p)
{
	Coord p1;
	p1 = dcordFunct_dksi(p);
	return -p1.z;
}


double FinitElemCulcer::dEqu1dnu(ECoord p)
{
	Coord p1;
	p1 = dcordFunct_dnu(p);
	return -p1.x;
}

double FinitElemCulcer::dEqu2dnu(ECoord p)
{
	Coord p1;
	p1 = dcordFunct_dnu(p);
	return -p1.y;
}

double FinitElemCulcer::dEqu3dnu(ECoord p)
{
	Coord p1;
	p1 = dcordFunct_dnu(p);
	return -p1.z;
}


double FinitElemCulcer::dEqu1detta(ECoord p)
{
	Coord p1;
	p1 = dcordFunct_detta(p);
	return -p1.x;
}

double FinitElemCulcer::dEqu2detta(ECoord p)
{
	Coord p1;
	p1 = dcordFunct_detta(p);
	return -p1.y;
}
double FinitElemCulcer::dEqu3detta(ECoord p)
{
	Coord p1;
	p1 = dcordFunct_detta(p);
	return -p1.z;
}


double FinitElemCulcer::countFunct(Coord p, bool& isInFinitElem)
{
	pointSolut = p;
	ECoord initMean(0, 0, 0);
	ECoord ans;
	if (countSolution(initMean, ans))
	{
		if (ans.ksi < 0.0 || ans.ksi > 1.0 || ans.nu < 0.0 || ans.nu > 1.0 || ans.etta < 0.0 || ans.etta > 1.0)
		{
			isInFinitElem = false;
			return 0.0;
		}
		isInFinitElem = true;
		return culcPhi(ans);
	}
	isInFinitElem = false;
	return 0.0;
	
}

