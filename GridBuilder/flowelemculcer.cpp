#include "flowelemculcer.h"
#include "array.h"
#include "slae.h"
#include "vectgeom.h"
#include "DifferentEquParams.h"

void FlowElemCulcer::init(FinitElement finitElem, CoordStorage coordsStore, double* q)
{
	this->finitElem = finitElem;

	int globInd[VER_NUM];
	arrayspace::copy(globInd, finitElem.ver, VER_NUM);

	coordsStore.copyElemsByInd(globInd, locCoords, VER_NUM);
	arrayspace::copyElemsByInd(globInd, locQ, q, VER_NUM);

}

void FlowElemCulcer::prepDataForGetFaceFlow(int iFace)
{
	finitElem.getFaceLocalNum(iFace, locFaceInd);
	finitElem.getOpposFaceLocalNum(iFace, locOpFaceInd);
	Coord coordBuf;
	double qBuf;

	for (int i = 0, j = VER_NUM_FACE; i < VER_NUM_FACE; i++, j++)
	{
		coords[i] = locCoords[locFaceInd[i]];
		coords[j] = locCoords[locOpFaceInd[i]];

		q[i] = locQ[locFaceInd[i]];
		q[j] = locQ[locOpFaceInd[i]];
	}
}

double FlowElemCulcer::integrFunct(Coord eCord)
{
	ECoord p;
	p = eCord;
	p.etta = 0;

	ECoord grad = culcGrad(p);

	buildJ(J, p);
	
	double gradB[DIM];
	grad.transIntoMas(gradB);

	double x[DIM];
	slae::solveSLAU3(J, gradB, x);

	Coord gradXYZ;
	gradXYZ = x;

	Vect n = findNormVect(coords[0], coords[1], coords[2]);
	gradXYZ *= n;

	faceBasic.init(coords);
	double detJ2D = faceBasic.getDetJ(p.ksi, p.nu);
	
	return (gradXYZ.x + gradXYZ.y + gradXYZ.z) * detJ2D;
}

bool FlowElemCulcer::isInnerN(int iLocalFace)
{
	switch (iLocalFace)
	{
	case 0:
		return true;
	case 1:
		return false;
	case 2:
		return true;
	case 3:
		return true;
	case 4:
		return false;
	case 5:
		return false;
	}
}

void FlowElemCulcer::getFlows(double flows[FACES_NUM], signed char flowsSigns[FACES_NUM])
{
	double flow;
	Coord left(0, 0, 0), right(1, 1, 1);

	for (int i = 0; i < FACES_NUM; i++)
	{
		prepDataForGetFaceFlow(i);
		flow = -integrateGausse(left, right);
		
		flow *= isInnerN(i) ? -1 : 1;

		flowsSigns[i] = flow < 0 ? -1: 1;
		flows[i] = abs(flow);
	}
}
