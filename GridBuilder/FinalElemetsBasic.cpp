#include "FinalElemetsBasic.h"
#include "slae.h"
#include "fullmatrix.h"
#include "vectgeom.h"
#include "DifferentEquParams.h"
#include "array.h"
#include <iostream>


void FinalElemBase::setFinalElement(FinitElement finalElement, CoordStorage knotsStorage)
{
	this->finitElement = finalElement;
	int i;
	for (i = 0; i < N_KNOTS_3D; i++)
		globalInd[i] = finalElement.ver[i];

	for (i = 0; i < N_KNOTS_3D; i++)
		coords[i] = knotsStorage.coords[globalInd[i]];

	/*Пока без понятия как там рополагаются x и y так что это место довольно опасное, если работать не будет, то  возможно ошибка здесь*/
		
}

void FinalElemBase::buildLocalM(double M[N_KNOTS_3D][N_KNOTS_3D], double gamma)
{
	int i, j;
	Coord left(0, 0, 0);
	Coord right(1, 1, 1);
	for (i = 0; i < N_KNOTS_3D; i++)
	{
		iPhi = i;
		for (j = i; j < N_KNOTS_3D; j++)
		{
			jPhi = j;
			//M[i][j] = gamma*LOCAL_MATRIX_M_1[i % 2][j % 2] * LOCAL_MATRIX_M_1[i / 2][j / 2] * hx * hy / 36.0;
			M[i][j] = gamma * integrateGausse(left, right);;
			M[j][i] = M[i][j];
		}
	}
	/*
	M[0][0] = gamma*4 * hx * hy / 36.0;
	M[1][0] = M[0][1] = gamma *2 * hx * hy / 36.0;
	M[2][0] = M[0][2] = gamma * 2 * hx * hy / 36.0;
	M[3][0] = M[0][3] = gamma * 1 * hx * hy / 36.0;

	M[1][1] = gamma * 4 * hx * hy / 36.0;
	M[1][2] = gamma * 1 * hx * hy / 36.0;
	M[3][1] = M[1][3] = gamma * 2 * hx * hy / 36.0;

	M[2][2] = gamma *4 * hx * hy / 36.0;
	M[3][2] = M[2][3] = gamma * 2 * hx * hy / 36.0;

	M[3][3] = gamma * 4 * hx * hy / 36.0;
	*/
}

double FinalElemBase::integrFunct(Coord p)
{
	double mean;
	switch (integChoice)
	{
	case INTEG_G:
		return funcIntegrG(p);
	case INTEG_M:
		return funcIntegrM(p);
	}
	return  funcIntegrG(p);
}

double FinalElemBase::functIntegSqare(Coord eCoord)
{
	ECoord p;
	p = eCoord;

	buildJ(J, p);
	double DetJ = countDetMatrix3(J);

	return DetJ;
}

double FinalElemBase::calcSqare()
{
	integChoice = INTEG_M;

	Coord left(0, 0, 0);
	Coord right(1, 1, 1);

	return integrateGausse(left, right);
}

double FinalElemBase::funcIntegrG(Coord eCord)
{
	ECoord p;
	p = eCord;
	ECoord gradi, gradj;
	double grad_i[DIM], grad_j[DIM];
	double x1[DIM], x2[DIM];
	buildJ(J, p);

	//std::cout << "\t \t \t Расчёт градиента i" << std::endl;
	gradi = countGrad_i(p, iGradPhi);
	//std::cout << "\t \t \t Расчёт градиента j" << std::endl;
	gradj = countGrad_i(p, jGradPhi);

	//std::cout << gradi.ksi << " " << gradi.nu << " " << gradi.etta << std::endl;
	//std::cout << "\t \t \t Расчёт определителя матрицы Якоби" << std::endl;
	double DetJ = countDetMatrix3(J);
	//std::cout << "\t \t \t Расчёт определителя матрицы Якоби" << std::endl;
	gradi.transIntoMas(grad_i);
	gradj.transIntoMas(grad_j);

	//std::cout << "\t \t \tРешение СЛАУ 3x3 (Обратная матрица умноженная на вектор)" << std::endl;
	slae::solveSLAU3(J, grad_i, x1);
	slae::solveSLAU3(J, grad_j, x2);
	
	//std::cout << "\t \t \tItegrFunct End " << ans << std::endl;
	return  arrayspace::scal(x1, x2, DIM) * abs(DetJ);
}

double FinalElemBase::funcIntegrM(Coord eCord)
{
	ECoord p;
	p = eCord;
	double x1[DIM], x2[DIM];
	buildJ(J, p);

	//std::cout << gradi.ksi << " " << gradi.nu << " " << gradi.etta << std::endl;
	//std::cout << "\t \t \t Расчёт определителя матрицы Якоби" << std::endl;
	double DetJ = countDetMatrix3(J);
	//std::cout << "\t \t \t Расчёт определителя матрицы Якоби" << std::endl;
	//std::cout << "\t \t \tItegrFunct End " << ans << std::endl;
	return  phi_i(p, iGradPhi) * phi_i(p, iGradPhi) *arrayspace::scal(x1, x2, DIM) * abs(DetJ);
}

void FinalElemBase::buildLocalG(double G[N_KNOTS_3D][N_KNOTS_3D])
{
	int i, j;
	Coord left(0, 0, 0);
	Coord right(1, 1, 1);
	integChoice = INTEG_G;
	//integrFunct = &funcIntegrG;
	for (i = 0; i < N; i++)
	{
		iGradPhi = i;
		for (j = i; j < N; j++)
		{
			jGradPhi = j;
			//std::cout << "\t \t Расчёт интеграла" << std::endl;
			G[j][i] = G[i][j] = DifferentEquParams::lambda(finitElement.K, finitElement.phaseStorage.phases, finitElement.phaseStorage.count)*integrateGausse(left, right);
			//std::cout << "\t \t Расчёт интеграла закончен" << std::endl;
		}		
	}
}


void FinalElemBase::buildLocalB(double B[N_KNOTS_3D])
{
	double M[N_KNOTS_3D][N_KNOTS_3D];
	double F[N_KNOTS_3D];
	buildLocalM(M, 1.0);
	fillF(F);
	multMatrixVect(M, F, B);
}

void FinalElemBase::fillF(double F[N_KNOTS_3D])
{	
	for (int i = 0; i < N_KNOTS_3D; i++)
	    F[i] = DifferentEquParams::f(coords[i].x, coords[i].y, coords[i].z);
}

#if UNLINEAR_TASK
double FinalElemBase::count_dG_dq(int i, int r, int j)
{
	return dlambda_du(q[j], xKnot[j % 2], yKnot[j / 2]) * countGausse2ForG(i, r, j);
}

void FinalElemBase::addNewt(double G[N_KNOTS_3D][N_KNOTS_3D], double B[N_KNOTS_3D])
{
	double sum;
	double sumB;
	int i, j, r;
	for (i = 0; i < N_KNOTS_3D; i++)
	{
		sumB = 0;
		for (j = 0; j < N_KNOTS_3D; j++)
		{
			sum = 0;
			for (r = 0; r < N_KNOTS_3D; r++)
			{
				sum += count_dG_dq(i, r, j)*q[r];
			}
			G[i][j] += sum;
			sumB += sum * q[j];
		}
		B[i] += sumB;
	}
}

void FinalElemBase::buildNewt(double G[N_KNOTS_3D][N_KNOTS_3D], double B[N_KNOTS_3D])
{
	buildLocalG(G);
	buildLocalB(B);
	//addNewt(G, B);
}
#endif
