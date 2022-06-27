#pragma once
#include "gridbuilder.h"
#include "finitelemculcer.h"
#include "gaussintegr.h"
#include "facebasic.h"

class FlowElemCulcer : private FinitElemCulcer, public GausseIntegr2D
{
public:
	void init(FinitElement finitElem, CoordStorage coordsStore, double* q) override;
	void getFlows(double flows[FACES_NUM], signed char flowsSigns[FACES_NUM]);
private:
	FaceBasic faceBasic;
	FinitElement finitElem;
	Coord locCoords[VER_NUM];
	FaceStore faceStore;
	double locQ[VER_NUM];

	int locFaceInd[VER_NUM_FACE];
	int locOpFaceInd[VER_NUM_FACE];

	void prepDataForGetFaceFlow(int iFace);
	virtual double integrFunct(Coord p) override;
	bool isInnerN(int iLocalFace);
};