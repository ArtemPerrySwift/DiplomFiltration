#include "facebasic.h"
#include "DifferentEquParams.h"
#include "fullmatrix.h"
#include "slae.h"
#include "vectgeom.h"

int FaceBasic::u(int i) { return i % 2; }
int FaceBasic::v(int i) { return (i / 2) % 2; }

double FaceBasic::phi2D(double ksi, double nu, int i)
{
	double res = 1;
	res *= (u(i) ? ksi : 1 - ksi);
	res *= (v(i) ? nu : 1 - nu);

	return res;
}

double FaceBasic::dphi_dksi2D(double ksi, double nu, int i)
{
	double res = 1;
	res *= (u(i) ? 1 : -1);
	res *= (v(i) ? nu : 1 - nu);

	return res;
}

double FaceBasic::dphi_dnu2D(double ksi, double nu, int i)
{
	double res = 1;
	res *= (u(i) ? ksi : 1 - ksi);
	res *= (v(i) ? 1 : -1);

	return res;
}

Coord FaceBasic::dcordFunct_dksi2D(double ksi, double nu)
{
	Coord p;
	for (int i = 0; i < N_BASIC; i++)
		p += coordFace[i] * dphi_dksi2D(ksi, nu, i);

	return p;
}

Coord FaceBasic::dcordFunct_dnu2D(double ksi, double nu)
{
	Coord p;
	for (int i = 0; i < N_BASIC; i++)
		p += coordFace[i] * dphi_dnu2D(ksi, nu, i);

	return p;
}

Coord FaceBasic::cordFunct2D(double ksi, double nu)
{
	Coord p;
	for (int i = 0; i < N_BASIC; i++)
		p += coordFace[i] * phi2D(ksi, nu, i);

	return p;
}

double  FaceBasic::getDetJ(double ksi, double nu)
{
	Coord pdksi = dcordFunct_dksi2D(ksi, nu);
	Coord pdnu = dcordFunct_dnu2D(ksi, nu);
	//J[0][0] = p.x;
	//J[0][1] = p.y;
	//J[0][2] = p.z;
	double detBuf = 0, sum = 0;

	detBuf = pdksi.x * pdnu.y - pdksi.y * pdnu.x;
	sum += detBuf * detBuf;

	detBuf = pdksi.y * pdnu.z - pdksi.z * pdnu.y;
	sum += detBuf * detBuf;

	detBuf = pdksi.x * pdnu.z - pdksi.z * pdnu.x;
	sum += detBuf * detBuf;

	return sqrt(sum);
}


double FaceBasic2Cond::integrFunct(Coord p)
{
	double ksi = p.x,
		   nu = p.y;

	double DetJ = getDetJ(ksi, nu);
	return phi2D(ksi, nu, iPhi2d) * phi2D(ksi, nu, jPhi2D) * DetJ;
}

void FaceBasic::init(Face face, CoordStorage coordStorage)
{
	int* iKnots = face.knots;
	Coord* coords = coordStorage.coords;
	for (int i = 0; i < N_BASIC; i++)
		coordFace[i] = coords[iKnots[i]];
}

void FaceBasic::init(Coord coords[N_BASIC])
{
	for (int i = 0; i < N_BASIC; i++)
		coordFace[i] = coords[i];
}


int FaceBasic2Cond::buildLocalM2D()
{
	int i, j;
	Coord left(0, 0, 0);
	Coord right(1, 1, 1);
	for (i = 0; i < N_BASIC; i++)
	{
		iPhi2d = i;
		for (j = i; j < N_BASIC; j++)
		{
			jPhi2D = j;
			LocalC[j][i] = LocalC[i][j] = integrateGausse(left, right);
		}			
	}
	return 0;
}

void FaceBasic2Cond::fillLocalEtta()
{
	for (int i = 0; i < N_BASIC; i++)
		LocalEtta[i] = DifferentEquParams::du_dn(coordFace[i].x, coordFace[i].y, coordFace[i].z);
}

void FaceBasic2Cond::calcAddLocalB(double localAddB[N_BASIC])
{
	buildLocalM2D();
	fillLocalEtta();
	multMatrixVect(LocalC, LocalEtta, localAddB);
}

void FaceBasicFlow::countGrad_i2D(double grad_i[2], int i, double ksi, double nu)
{
	grad_i[0] = dphi_dksi2D(ksi, nu, i);
	grad_i[1] = dphi_dnu2D(ksi, nu, i);
}

void FaceBasicFlow::init(Face face, CoordStorageXYZ coordStorage, double* q, bool isInnerNorm, double lambda)
{
	this->isInnerNorm = isInnerNorm;
	this->lambda = lambda;
	int* iKnots = face.knots;
	Coord* coords = coordStorage.coords;
	int iGlobal;
	for (int i = 0; i < N_BASIC; i++)
	{
		iGlobal = iKnots[i];
		coordFace[i] = coords[iGlobal];
		localQ[i] = q[iGlobal];
	}		
}

void FaceBasicFlow::buildJ2D(double J[3][DIM], double ksi, double nu)
{
	Coord pdksi = dcordFunct_dksi2D(ksi, nu);
	Coord pdnu = dcordFunct_dnu2D(ksi, nu);

	J[0][0] = pdksi.x;
	J[0][1] = pdnu.x;

	J[1][0] = pdksi.y;
	J[1][1] = pdnu.y;

	J[2][0] = pdksi.z;
	J[2][1] = pdnu.z;
}
/*
double FaceBasicFlow::integrFunct(Coord p)
{
	double ksi = p.x,
		   nu = p.y;
	Coord p1;
	double grad_i2D[DIM];
	double grad_sum2D[DIM];

	for (int j = 0; j < 2; j++)
		grad_sum2D[j] = 0;

	//Рассчёт суммы градиентов
	for (int i = 0; i < N_BASIC; i++)
	{
		countGrad_i2D(grad_i2D, i, ksi, nu);
		for (int j = 0; j < DIM; j++)
			grad_sum2D[j] += localQ[i] * grad_i2D[j];
	}

	double J2D[DIM_XYZ][2];
	buildJ2D(J2D, ksi, nu);

	double slae3[DIM_XYZ][DIM_XYZ];
	multMatrix3x2(J2D, J2D, slae3);

	double b3[DIM_XYZ];
	multMatrixVect(J2D, grad_sum2D, b3);

	double matBuf[DIM_XYZ][DIM_XYZ], bufv1[DIM_XYZ], slaeSol[DIM_XYZ];
	slae::solveSLAU3(slae3, matBuf, b3, bufv1, slaeSol);

	p1 = slaeSol;
	double detJ2D = getDetJ(ksi, nu);

	Vect n = findNormVect(coordFace[0], coordFace[2], coordFace[3]);
	if (isInnerNorm) n *= -1;
	p1 *= n;

	return (p1.x + p1.y + p1.z) * detJ2D;
}
*/
double FaceBasicFlow::getFlow()
{
	Coord left(0, 0, 0);
	Coord right(1, 1, 1);

	return lambda; //* integrateGausse(left, right);
}
