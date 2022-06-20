#pragma once
#include "coord.h"
#include "ecoord.h"

using namespace ecoord;

class BasicFunct3D
{
public:
	static const int N = 8;
	static const int DIM = 3;

	static double phi_i(ECoord p, unsigned char i);

	static double dphi_dksi_i(ECoord p, unsigned char i);
	static double dphi_dnu_i(ECoord p, unsigned char i);
	static double dphi_detta_i(ECoord p, unsigned char i);

	static double coordFunct_i(ECoord p, unsigned char i);

	static double dcoordFunct_dksi_i(ECoord p, unsigned char i);
	static double dcoordFunct_dnu_i(ECoord p, unsigned char i);
	static double dcoordFunct_detta_i(ECoord p, unsigned char i);

	Coord dcordFunct_dksi(ECoord p);
	Coord dcordFunct_dnu(ECoord p);
	Coord dcordFunct_detta(ECoord p);

	Coord cordFunct(ECoord p);
	double coordFunct_x(ECoord p);
	double coordFunct_y(ECoord p);
	double coordFunct_z(ECoord p);

	static ECoord countGrad_i(ECoord p, unsigned char i);

	void buildJ(double J[DIM][DIM], ECoord p);

	virtual void init(Coord coords[N]);

protected:
	Coord coords[N]; // ¬ершины отображаемого шестигранника
	double J[DIM][DIM];

	static int u(unsigned char i);
	static int v(unsigned char i);
	static int g(unsigned char i);
};
