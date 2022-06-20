#include "basicfunct3d.h"
#include <iostream>
int BasicFunct3D::u(unsigned char i) { return i % 2; }
int BasicFunct3D::v(unsigned char i) { return (i / 2) % 2; }
int BasicFunct3D::g(unsigned char i) { return i / 4; }


ECoord BasicFunct3D::countGrad_i(ECoord p, unsigned char i)
{
	ECoord grad;
	grad.ksi = dphi_dksi_i(p, i);
	grad.nu = dphi_dnu_i(p, i);
	grad.etta = dphi_detta_i(p, i);

	return grad;
}

double BasicFunct3D::phi_i(ECoord p, unsigned char i)
{
	double res = 1;
	res *= (u(i) ? p.ksi : 1 - p.ksi);
	res *= (v(i) ? p.nu : 1 - p.nu);
	res *= (g(i) ? p.etta : 1 - p.etta);

	return res;
}

double BasicFunct3D::coordFunct_i(ECoord p, unsigned char i)
{
	double res = 1;
	res *= (u(i) ? p.ksi : 1 - p.ksi);
	res *= (v(i) ? p.nu : 1 - p.nu);
	res *= (g(i) ? p.etta : 1 - p.etta);

	return res;
}

double BasicFunct3D::dcoordFunct_dksi_i(ECoord p, unsigned char i)
{
	double res = 1;
	res *= (u(i) ? 1 : -1);
	res *= (v(i) ? p.nu : 1 - p.nu);
	res *= (g(i) ? p.etta : 1 - p.etta);

	return res;
}

double BasicFunct3D::dcoordFunct_dnu_i(ECoord p, unsigned char i)
{
	double res = 1;
	res *= (u(i) ? p.ksi : 1 - p.ksi);
	res *= (v(i) ? 1 : -1);
	res *= (g(i) ? p.etta : 1 - p.etta);

	return res;
}

double BasicFunct3D::dcoordFunct_detta_i(ECoord p, unsigned char i)
{
	double res = 1;
	res *= (u(i) ? p.ksi : 1 - p.ksi);
	res *= (v(i) ? p.nu : 1 - p.nu);
	res *= (g(i) ? 1 : -1);
	return res;
}

double BasicFunct3D::dphi_dksi_i(ECoord p, unsigned char i)
{
	double res = 1;
	res *= (u(i) ? 1 : -1);
	res *= (v(i) ? p.nu : 1 - p.nu);
	res *= (g(i) ? p.etta : 1 - p.etta);

	return res;
}

double BasicFunct3D::dphi_dnu_i(ECoord p, unsigned char i)
{
	double res = 1;
	res *= (u(i) ? p.ksi : 1 - p.ksi);
	res *= (v(i) ? 1 : -1);
	res *= (g(i) ? p.etta : 1 - p.etta);

	return res;
}

double BasicFunct3D::dphi_detta_i(ECoord p, unsigned char i)
{
	double res = 1;
	res *= (u(i) ? p.ksi : 1 - p.ksi);
	res *= (v(i) ? p.nu : 1 - p.nu);
	res *= (g(i) ? 1 : -1);
	return res;
}

Coord BasicFunct3D::dcordFunct_dksi(ECoord p)
{
	Coord dCoord_dKsi;
	//std::cout << "DcordDKsi" << std::endl;
	for (int i = 0; i < N; i++)
	{
		//std::cout << dCoord_dKsi.x << " " << dCoord_dKsi.y << " " << dCoord_dKsi.z << std::endl;
		//std::cout << "Coord " << coords[i].x << " " << coords[i].y << " " << coords[i].z << "DpiDksi " << dphi_dksi(p, i) << std::endl;
		dCoord_dKsi += coords[i] * dcoordFunct_dksi_i(p, i);
	}
	//std::cout << dCoord_dKsi.x << " " << dCoord_dKsi.y << " " << dCoord_dKsi.z << std::endl;

	return dCoord_dKsi;
}

Coord BasicFunct3D::dcordFunct_dnu(ECoord p)
{
	Coord dCoord_dNu;
	for (int i = 0; i < N; i++)
		dCoord_dNu += coords[i] * dcoordFunct_dnu_i(p, i);

	return dCoord_dNu;
}

Coord BasicFunct3D::dcordFunct_detta(ECoord p)
{
	Coord dCoord_dEtta;
	for (int i = 0; i < N; i++)
		dCoord_dEtta += coords[i] * dcoordFunct_detta_i(p, i);

	return dCoord_dEtta;
}

Coord BasicFunct3D::cordFunct(ECoord p)
{
	Coord coordXYZ;
	for (int i = 0; i < N; i++)
		coordXYZ += coords[i] * coordFunct_i(p, i);

	return coordXYZ;
}

void BasicFunct3D::buildJ(double J[DIM][DIM], ECoord p)
{
	//std::cout << "Coords" << std::endl;
	//for (int i = 0; i < 8; i++)
	///	std::cout << coords[i].x << " " << coords[i].y << " " << coords[i].z << std::endl;

	//std::cout << "Coord P" << std::endl;
	//std::cout << p.ksi << " " << p.nu << " " << p.etta << std::endl;
	Coord dCordXYZ = dcordFunct_dksi(p);
	J[0][0] = dCordXYZ.x;
	J[0][1] = dCordXYZ.y;
	J[0][2] = dCordXYZ.z;

	dCordXYZ = dcordFunct_dnu(p);
	J[1][0] = dCordXYZ.x;
	J[1][1] = dCordXYZ.y;
	J[1][2] = dCordXYZ.z;

	dCordXYZ = dcordFunct_detta(p);
	J[2][0] = dCordXYZ.x;
	J[2][1] = dCordXYZ.y;
	J[2][2] = dCordXYZ.z;
}

void BasicFunct3D::init(Coord coords[N])
{
	for (int i = 0; i < N; i++)
		this->coords[i] = coords[i];
	
}

double BasicFunct3D::coordFunct_x(ECoord p)
{
	double x = 0;
	for (int i = 0; i < N; i++)
		x += coords[i].x * coordFunct_i(p, i);

	return x;
}
double BasicFunct3D::coordFunct_y(ECoord p)
{
	double y = 0;
	for (int i = 0; i < N; i++)
		y += coords[i].y * coordFunct_i(p, i);

	return y;
}
double BasicFunct3D::coordFunct_z(ECoord p)
{
	double z = 0;
	for (int i = 0; i < N; i++)
		z += coords[i].z * coordFunct_i(p, i);

	return z;
}