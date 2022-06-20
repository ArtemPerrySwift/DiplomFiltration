#include "sne.h"
#include "programlog.h"
#include "array.h"
#include "slae.h"

SNE3D::SNE3D()
{
	maxiter = 1000;
	err = 1e-14;
	bettaEps = 1e-14;
	epsDif = 1e-7;
	minMean = 1e-14;
	isNumericalDif = true;
}
void SNE3D::fillJacobMatrix(double JacMatrix[DIMMEN][DIMMEN], ECoord p)
{
	JacMatrix[0][0] = dEqu1dksi(p);
	JacMatrix[0][1] = dEqu1dnu(p);
	JacMatrix[0][2] = dEqu1detta(p);

	JacMatrix[1][0] = dEqu2dksi(p);
	JacMatrix[1][1] = dEqu2dnu(p);
	JacMatrix[1][2] = dEqu2detta(p);

	JacMatrix[2][0] = dEqu3dksi(p);
	JacMatrix[2][1] = dEqu3dnu(p);
	JacMatrix[2][2] = dEqu3detta(p);
}

void SNE3D::countF(double F[DIMMEN], ECoord p)
{
	F[0] = -equ1(p);
	F[1] = -equ2(p);
	F[2] = -equ3(p);
}

bool SNE3D::setMaxiter(int maxiter)
{
	if (maxiter < 1)
	{
		programlog::writeErr("SNE3D - Attempting to set the maximum number of iterations parameter to less than 1");
		return false;
	}
	this->maxiter = maxiter;
	return true;
}
bool SNE3D::setErr(double err)
{
	if (err <= 0.0)
	{
		programlog::writeErr("SNE3D - An attempt to assign a value equal to 0 or less than 0 to the parameter of the permissible error in solving a system of nonlinear equations");
		return false;
	}
	this->err = err;
	return true;
}
bool SNE3D::setBettaEps(double bettaEps)
{
	if (bettaEps <= 0.0)
	{
		programlog::writeErr("SNE3D - An attempt to assign to the parameter of the minimum allowable change in the solution of a system of nonlinear equations at each iteration a value equal to 0 or less than 0");
		return false;
	}
	this->bettaEps = bettaEps;
	return true;
}
bool SNE3D::setEpsDif(double epsDif)
{
	if (epsDif <= 0.0)
	{
		programlog::writeErr("SNE3D - An attempt to assign a value equal to 0 or less than 0 to the step parameter of the numerical calculation of derivatives");
		return false;
	}
	this->epsDif = epsDif;
	return true;
}

bool SNE3D::setMinEpsDif(double minEpsDif)
{
	if (minEpsDif <= 0.0)
	{
		programlog::writeErr("SNE3D - Attempt to assign a value equal to 0 or less than 0 to the parameter of the minimum step of the numerical calculation of derivatives");
		return false;
	}
	this->epsDif = epsDif;
	return true;
}

bool SNE3D::setMinMean(double minMean)
{
	if (minMean <= 0.0)
	{
		programlog::writeErr("SNE3D - Attempt to set the minimum number parameter not considered 0 to a value equal to 0 or less than 0");
		return false;
	}
	this->minMean = minMean;
	return true;
}

double SNE3D::dEqu1dksi(ECoord p)
{
	ECoord p1 = p;
	double d = (abs(p.ksi) < minMean) ? minMean : p.ksi * epsDif;
	p1.ksi += d;
	return (equ1(p1) - equ1(p)) / d;
}

double SNE3D::dEqu1dnu(ECoord p)
{
	ECoord p1 = p;
	double d = (abs(p.nu) < minMean) ? minMean : p.nu * epsDif;
	p1.nu += d;
	return (equ1(p1) - equ1(p)) / d;
}

double SNE3D::dEqu1detta(ECoord p)
{
	ECoord p1 = p;
	double d = (abs(p.etta) < minMean) ? minMean : p.etta * epsDif;
	p1.etta += d;
	return (equ1(p1) - equ1(p)) / d;
}

double SNE3D::dEqu2dksi(ECoord p)
{
	ECoord p1 = p;
	double d = (abs(p.ksi) < minMean) ? minMean : p.ksi * epsDif;
	p1.ksi += d;
	return (equ2(p1) - equ2(p)) / d;
}

double SNE3D::dEqu2dnu(ECoord p)
{
	ECoord p1 = p;
	double d = (abs(p.nu) < minMean) ? minMean : p.nu * epsDif;
	p1.nu += d;
	return (equ2(p1) - equ2(p)) / d;
}

double SNE3D::dEqu2detta(ECoord p)
{
	ECoord p1 = p;
	double d = (abs(p.etta) < minMean) ? minMean : p.etta * epsDif;
	p1.etta += d;
	return (equ2(p1) - equ2(p)) / d;
}

double SNE3D::dEqu3dksi(ECoord p)
{
	ECoord p1 = p;
	double d = (abs(p.ksi) < minMean) ? minMean : p.ksi * epsDif;
	p1.ksi += d;
	return (equ3(p1) - equ3(p)) / d;
}

double SNE3D::dEqu3dnu(ECoord p)
{
	ECoord p1 = p;
	double d = (abs(p.nu) < minMean) ? minMean : p.nu * epsDif;
	p1.nu += d;
	return (equ3(p1) - equ3(p)) / d;
}

double SNE3D::dEqu3detta(ECoord p)
{
	ECoord p1 = p;
	double d = (abs(p.etta) < minMean) ? minMean : p.etta * epsDif;
	p1.etta += d;
	return (equ3(p1) - equ3(p)) / d;
}


bool SNE3D::countSolution(ECoord initAppr, ECoord& answer)
{
	countF(F, initAppr);
	double curErr = arrayspace::scal(F, F, DIMMEN);
	double betta = 1;
	ECoord p = initAppr;
	int i;
	if(isNumericalDif)
	{
		double begEpsDif = epsDif;
		for (i = 0; i < maxiter && curErr > err && betta > bettaEps && epsDif > minEpsDif; i++)
		{
			betta = 1;
			fillJacobMatrix(JacMatrix, p);
			if (slae::solveSLAU3(JacMatrix, F, ans) == -1)
			{	
				while (epsDif > minEpsDif && slae::solveSLAU3(JacMatrix, F, ans) == -1)
				{
					epsDif /= 2.0;
					fillJacobMatrix(JacMatrix, p);
				}

				if (epsDif > minEpsDif)
				{
					epsDif = begEpsDif;
					p = getIterSolutUsingBetta(p, betta, curErr);
				}
				else
				{
					programlog::writeErr("SLAE made in SNU cannot be solved");
					return false;
				}
					
					
			}
			else
				p = getIterSolutUsingBetta(p, betta, curErr);
		}
	}
	else
	{
		for (i = 0; i < maxiter && curErr > err && betta > bettaEps; i++)
		{
			betta = 1;
			fillJacobMatrix(JacMatrix, p);
			if (slae::solveSLAU3(JacMatrix, F, ans) == -1)
			{
				programlog::writeErr("SLAE made in SNU cannot be solved");
				return false;
			}
				
			p = getIterSolutUsingBetta(p, betta, curErr);
		}
	}

	if (i == maxiter)
	{
		programlog::writeErr("Solving of SNE reached the limit of iterations");
		return false;
	}

	if (betta <= bettaEps)
	{
		//programlog::writeErr("Solving of SNE stocked");
		return false;
	}

	answer = p;
	return true;
}

ECoord SNE3D::getIterSolutUsingBetta(ECoord p,double& betta, double& normF)
{
	ECoord p1;
	ECoord dp;
	dp = ans;
	double bettaBuf = betta;
	double normF0 = normF;
	p1 = p + dp * bettaBuf;
	countF(F, p1);
	double normCurF = arrayspace::scal(F, F, DIMMEN);
	while (normCurF > normF0 && bettaBuf > bettaEps)
	{
		bettaBuf /= 2;
		p1 = p + dp * bettaBuf;
		countF(F, p1);
		normCurF = arrayspace::scal(F, F, DIMMEN);
	}

	if (normCurF < normF)
	{
		normF = normCurF;
		betta = bettaBuf;
		return p1;
	}
	else
	{
		betta = 0;
		return p;
	}
}

