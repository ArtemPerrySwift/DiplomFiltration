#pragma once
#include "basicfunct3d.h"
#include "gridbuilder.h"
#include "sne.h"

class FinitElemCulcer: public BasicFunct3D, private SNE3D
{
public:
	double culcPhi(ECoord p);
	ECoord culcGrad(ECoord p);
	virtual void init(FinitElement finitElem, CoordStorage coordsStore, double* q);
	double countFunct(Coord p, bool& isInFinitElem);
	FinitElemCulcer();
protected:
	double q[N];
	FinitElement finitElement;

private:
	Coord pointSolut;
	double equ1(ECoord p) override;
	double equ2(ECoord p) override;
	double equ3(ECoord p) override;

	double dEqu1dksi(ECoord p) override;
	double dEqu2dksi(ECoord p) override;
	double dEqu3dksi(ECoord p) override;

	double dEqu1dnu(ECoord p) override;
	double dEqu2dnu(ECoord p) override;
	double dEqu3dnu(ECoord p) override;

	double dEqu1detta(ECoord p) override;
	double dEqu2detta(ECoord p) override;
	double dEqu3detta(ECoord p) override;
};
