#include <iostream>
//#include "mesh.h"
#include "matrix.h"
#include "LOS.h"
#include "difequ.h"
#include <string>

using namespace std;

const int POW_GAUSSE = 5;
const int DIM = 3;
const int N2 = 4;
const int N = 8;
const int N1 = 2;
const double R_GAUSSE[POW_GAUSSE] = {-0.90617984593866399279763,-0.5384693101056830910363, 0, 0.5384693101056830910363,0.90617984593866399279763};
const double W_GAUSSE[POW_GAUSSE] = { 0.23692688505618908751426, 0.4786286704993664680412, 0.5688888888888888888889 , 0.4786286704993664680412, 0.23692688505618908751426 };

const double M2_1[N2][N2] = { 4, 2, 2, 1,
							2, 4, 1, 2,
							2, 1, 4, 2,
							1, 2, 2, 4};

const double M2_2[N2][N2] = { 2, 2, 1, 1,
							2, 6, 1, 3,
							1, 1, 2, 2,
							1, 3, 2, 6 };

const double M2_3[N2][N2] = { 2, 1, 2, 1,
							1, 2, 1, 2,
							2, 1, 6, 3,
							1, 2, 3, 6 };

double funct_u1(Coord point) { return point.x; }
double funct_f(Coord point) { return point.x; }
double funct_du_dn(Coord point) { return 1; }
double funct_u_b(Coord point) { return 1 + point.x; }
double integrateGausse3(Coord l, Coord R,  Coord points[N], int i, int j, double (*funct)(Coord, Coord[N], int, int));

//double function(Coord point);
double phi(double ksi, double nu, double etta, int i);
double dphi_dksi(double ksi, double nu, double etta, int i);
double dphi_dnu(double ksi, double nu, double etta, int i);
double dphi_detta(double ksi, double nu, double etta, int i);
Coord dcordFunct_dksi(Coord cordF[N], double ksi, double nu, double etta);
void buildJ(double J[DIM][DIM], Coord F[N], double ksi, double nu, double etta);
int solveSLAU3(double A[DIM][DIM], double BufM[DIM][DIM], double b[DIM], double bufb[DIM], double x[DIM]);
void countGrad_i(double grad_i[DIM], int i, double ksi, double nu, double etta);
double countScalV3(double v1[DIM], double v2[DIM]);
double countDetMatrix3(double J[DIM][DIM]);
double funcIntegrM(Coord transform, Coord F[N], int i, int j);
double funcIntegrG(Coord transform, Coord F[N], int i, int j);
int getGlobeNum(int i, int j, int k, int n_x, int n_y, int n_z) { return k * n_x * n_y + j * n_x + i; };
int u(int i);
int v(int i);
int g(int i);
int buildLocalM(double M[N][N], Coord F[N]);
int buildLocalG(double G[N][N], Coord F[N]);
int buildLocalMFor2Cond(double M2[N2][N2], double cord1[N2], double cord2[N2]);
int multiplyMatrix8(double M[N][N], double v[N], double res[N]); // Умножение матрицы 8*8 на вектор размерности 8
int multiplyMatrix4(double M[N2][N2], double v[N2], double res[N2]); // Умножение матрицы 4*4 на вектор размерности 4
void buildLocalB(double C[N][N], double localF[N], double localB[N], Coord F[N]);
void fillPointsOfElem(Coord loc[N], int i, int j, int k, StoreMeshKnots& XYZ);
int buildGlobalAandB(matrix& A, CalculationArea& mesh, double* b);
int addFirstConditions(matrix& A, double* b, CalculationArea& mesh);
int addSecondConditions(matrix& A, double* b, CalculationArea& mesh);
int addThirdConditions(matrix& A, double* b, CalculationArea& mesh);
void initVector(double* v, int n, double mean);
void printVector(double* q, int n, string name);
void printSolution(double* q, int nX, int nY, int nZ);
void printErr(double* q, double* X, double* Y, double* Z, int nX, int nY, int nZ);
void writeSolInVect(double* q, double* X, double* Y, double* Z, int nX, int nY, int nZ);
void writeComparison(double* x, double* y, int n, string name);

int main()
{
	CalculationArea сalculationArea("Cord.txt", "Sep.txt", "borders.txt", true);
	//mesh.readBorders("borders.txt");
	//mesh.writeMesh("FullCord.txt");
	
	matrix A(сalculationArea.XYZ);
	A.init();
	A.printMatrix();
	double* b = new double[A.n];
	initVector(b, A.n, 0);
	buildGlobalAandB(A, сalculationArea, b);
	double* x = new double[A.n];
	double* r = new double[A.n];
	double* z = new double[A.n];
	double* p = new double[A.n];
	double* buf = new double[A.n];
	//printVector(b, A.n, "B before conditions");
	addSecondConditions(A, b, сalculationArea);
	addThirdConditions(A, b, сalculationArea);
	addFirstConditions(A, b, сalculationArea);
	
	//printVector(b, A.n, "B after conditions");
	//A.printMatrix();
	initVector(x, A.n, 0);
	//writeSolInVect(r, mesh.X, mesh.Y, mesh.Z, mesh.nX, mesh.nY, mesh.nZ);
	mult_mat_vec(A, r, p);
	//printVector(p, A.n, "Supposed b");
	writeComparison(b, p, A.n, "Difference between fact and supposed b");
	count_LOS(A, x, r, z, p, b, buf, 10000, 1e-12);
	//printVector(x, A.n, "Solution");
	cout << "Solution: " << endl;
	printSolution(x, сalculationArea.XYZ.nX, сalculationArea.XYZ.nY, сalculationArea.XYZ.nZ);
	//printErr(x, mesh.X, mesh.Y, mesh.Z, mesh.nX, mesh.nY, mesh.nZ);

	//A.printMatrix();
	double s[DIM];
	
	return 0;
}
/*
double function(Coord p)
{
	return pow(p.x, 1)  + pow(p.y, 0) *p.z;
}
*/
int u(int i) { return i % 2; }
int v(int i) { return (i / 2) % 2; }
int g(int i) { return i / 4; }

double phi(double ksi, double nu, double etta, int i)
{
	double res = 1;
	res *= (u(i) ? ksi : 1 - ksi);
	res *= (v(i) ? nu : 1 - nu);
	res *= (g(i) ? etta : 1 - etta);
	
	return res;
}

double dphi_dksi(double ksi, double nu, double etta, int i)
{
	double res = 1;
	res *= (u(i) ? 1 : - 1);
	res *= (v(i) ? nu : 1 - nu);
	res *= (g(i) ? etta : 1 - etta);

return res;
}

double dphi_dnu(double ksi, double nu, double etta, int i)
{
	double res = 1;
	res *= (u(i) ? ksi : 1 - ksi);
	res *= (v(i) ? 1 : -1);
	res *= (g(i) ? etta : 1 - etta);

	return res;
}

double dphi_detta(double ksi, double nu, double etta, int i)
{
	double res = 1;
	res *= (u(i) ? ksi : 1 - ksi);
	res *= (v(i) ? nu : 1 - nu);
	res *= (g(i) ? 1 : -1);
	return res;
}

Coord dcordFunct_dksi(Coord cordF[N], double ksi, double nu, double etta)
{
	Coord p;
	for (int i = 0; i < N; i++)
		p += cordF[i] * dphi_dksi(ksi, nu, etta, i);

	return p;
}

Coord dcordFunct_dnu(Coord cordF[N], double ksi, double nu, double etta)
{
	Coord p;
	for (int i = 0; i < N; i++)
		p += cordF[i] * dphi_dnu(ksi, nu, etta, i);

	return p;
}

Coord dcordFunct_detta(Coord cordF[N], double ksi, double nu, double etta)
{
	Coord p;
	for (int i = 0; i < N; i++)
		p += cordF[i] * dphi_detta(ksi, nu, etta, i);

	return p;
}

Coord cordFunct(Coord cordF[N], double ksi, double nu, double etta)
{
	Coord p;
	for (int i = 0; i < N; i++)
		p += cordF[i] * phi(ksi, nu, etta, i);

	return p;
}

void buildJ(double J[DIM][DIM], Coord F[N], double ksi, double nu, double etta)
{
	Coord p = dcordFunct_dksi(F, ksi, nu, etta);
	J[0][0] = p.x;
	J[0][1] = p.y;
	J[0][2] = p.z;

	p = dcordFunct_dnu(F, ksi, nu, etta);
	J[1][0] = p.x;
	J[1][1] = p.y;
	J[1][2] = p.z;

	p = dcordFunct_detta(F, ksi, nu, etta);
	J[2][0] = p.x;
	J[2][1] = p.y;
	J[2][2] = p.z;
}

int buildLocalM(double M[N][N], Coord points[N])
{
	int i, j;
	Coord left(0, 0, 0);
	Coord right(1, 1, 1);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
			M[i][j] = integrateGausse3(left, right, points, i, j, funcIntegrM);
		
	}

	return 0;
}

int buildLocalG(double G[N][N], Coord points[N])
{
	int i, j;
	Coord left(0, 0, 0);
	Coord right(1, 1, 1);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
			G[i][j] = integrateGausse3(left, right, points, i, j, funcIntegrG);
	}

	return 0;
}

double countDetMatrix3(double J[DIM][DIM])
{
	return J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] + J[0][2] * J[1][0] * J[2][1] - J[0][2] * J[1][1] * J[2][0] - J[0][1] * J[1][0] * J[2][2] - J[0][0] * J[1][2] * J[2][1];
}

double integrateGausse3(Coord left, Coord right, Coord points[N], int i, int j, double (*funct)(Coord, Coord[N], int, int)) // Метод численного интегрирования Гаусса
{
	Coord avr = (left + right) / 2.0;
	Coord hCoordHalf = (right - left) / 2.0;
	Coord r_gausse;
	double sum = 0;
	//double avrX = (xL + xR) / 2.0, avrY = (yL + yR) / 2.0, avrZ = (zL + zR) / 2.0;
	//double hXhalf = (xR - xL) / 2.0, hYhalf = (yR - yL) / 2.0, hZhalf = (zR - zL) / 2.0;
	int l, m, k;
	for (l = 0; l < POW_GAUSSE; l++)
	{
		for (m = 0; m < POW_GAUSSE; m++)
		{
			for (k = 0; k < POW_GAUSSE; k++)
			{
				r_gausse = { R_GAUSSE[l], R_GAUSSE[m],  R_GAUSSE[k]};
				//sum += W_GAUSSE[l] * W_GAUSSE[m] * W_GAUSSE[k] * funct(avrX + R_GAUSSE[l] * hXhalf, avrY + R_GAUSSE[m] * hYhalf, avrZ + R_GAUSSE[k] * hZhalf, xF, yF, zF, i, j);
				sum += W_GAUSSE[l] * W_GAUSSE[m] * W_GAUSSE[k] * funct(avr + r_gausse*hCoordHalf, points, i, j);
			}
				
		}
	}

	//return hXhalf * hYhalf * hZhalf * sum;
	return hCoordHalf.x * hCoordHalf.y * hCoordHalf.z * sum;
}

double funcIntegrM(Coord transform, Coord F[N], int i, int j)
{
	double ksi = transform.x,
		   nu = transform.y,
		   etta = transform.z;
	double J[DIM][DIM];
	buildJ(J, F, ksi, nu, etta);
	double DetJ = countDetMatrix3(J);
	return phi(ksi, nu, etta, i) * phi(ksi, nu, etta, j) * DetJ;
}

double funcIntegrG(Coord transform, Coord F[N], int i, int j)
{
	double ksi = transform.x,
		   nu = transform.y,
		   etta = transform.z;
	double J[DIM][DIM], BufJ[DIM][DIM];
	double grad_i[DIM], grad_j[DIM];
	double bufv1[DIM], x1[DIM], x2[DIM];
	buildJ(J, F, ksi, nu, etta);
	countGrad_i(grad_i, i, ksi, nu, etta);
	countGrad_i(grad_j, j, ksi, nu, etta);
	double DetJ = countDetMatrix3(J);
	solveSLAU3(J, BufJ, grad_i, bufv1, x1);
	solveSLAU3(J, BufJ, grad_j, bufv1, x2);

	return  countScalV3(x1, x2)* DetJ;
}

void countGrad_i(double grad_i[DIM], int i, double ksi, double nu, double etta)
{
	grad_i[0] = dphi_dksi(ksi, nu, etta, i);
	grad_i[1] = dphi_dnu(ksi, nu, etta, i);
	grad_i[2] = dphi_detta(ksi, nu, etta, i);
}

double countScalV3(double v1[DIM], double v2[DIM])
{
	double sum = 0;
	for (int i = 0; i < DIM; i++)
		sum += v1[i] * v2[i];
	return sum;
}
int solveSLAU3(double A[DIM][DIM], double BufM[DIM][DIM], double b[DIM], double bufb[DIM], double x[DIM])
{
	int i, j, k;
	double buf;
	bool f;
	for (i = 0; i < DIM; i++)
	{
		bufb[i] = b[i];
		for (j = 0; j < DIM; j++)
			BufM[i][j] = A[i][j];
	}
		
	double a;
	for (i = 0; i < DIM; i++)
	{
		a = BufM[i][i];
		if(!a)
		{
			f = true;
			for (j = i + 1; f && j < DIM; j++)
				if (BufM[j][i]) f = false;

			if (f) return -1;
			for (k = i; k < DIM; k++)
			{
				buf = BufM[i][k];
				BufM[i][k] = BufM[j][k];
				BufM[j][k] = buf;
			}
			buf = bufb[i];
			bufb[i] = bufb[j];
			bufb[j] = buf;
		}

		for (j = i; j < DIM; j++)
			BufM[i][j] /= a;

		bufb[i] /= a;

		for (j = i + 1; j < DIM; j++)
		{
			a = BufM[j][i];
			for (k = i; k < DIM; k++)
			{
				BufM[j][k] -= a * BufM[i][k];
			}
			bufb[j] -= a * bufb[i];
		}
	}
	

	for (i = DIM; i > -1; i--)
	{
		for (j = i + 1; j < DIM; j++)
			bufb[i] -= BufM[i][j] * x[j];
		x[i] = bufb[i];
	}

	return 0;
}

int multiplyMatrix8(double M[N][N], double v[N], double res[N]) // Умножение матрицы 8*8 на вектор размерности 8
{
	int i, j;
	for (i = 0; i < N; i++)
	{
		res[i] = 0;
		for (j = 0; j < N; j++)
			res[i] += M[i][j] * v[j];
	}

	return 0;
}

int multiplyMatrix4(double M[N2][N2], double v[N2], double res[N2])
{
	int i, j;
	for (i = 0; i < N2; i++)
	{
		res[i] = 0;
		for (j = 0; j < N2; j++)
			res[i] += M[i][j] * v[j];
	}

	return 0;
}

void buildLocalB(double C[N][N], double localF[N], double localB[N], Coord coords[N])
{
	for (int i = 0; i < N; i++)
		localF[i] = funct_f(coords[i]);

	multiplyMatrix8(C, localF, localB);
}

void fillPointsOfElem(Coord loc[N], int i, int j, int k, StoreMeshKnots &XYZ)
{
	int p;
	int c, d, e, g;
	Coord* points = XYZ.data;
	for(c = 0, g = 0; c < N1; c++)
		for(d = 0; d < N1; d++)
			for (e = 0; e < N1; e++, g++)
			{
				p = XYZ.getKnotIndex(i + e, j + d, k + c);
				loc[g] = points[p];
				/*
				Xloc[g] = X[p];
				Yloc[g] = Y[p];
				Zloc[g] = Z[p];
				*/
			}
}

int buildGlobalAandB(matrix& A, CalculationArea& mesh, double* b)
{
	double G_local[N][N];
	double M_local[N][N];
	//double localC[N][N];
	double localF[N];
	double localB[N];
	//double Xloc[N], Yloc[N], Zloc[N];
	Coord loc[N];

	int L[N]; // Массив для храниения соответсвия глобальных номеров и локальных

	Area areaCur;

	int n = A.n;
	int n_gg = A.n_gg;
	double* di = A.di;
	double* gg = A.gg;
	int* ig = A.ig;
	int* jg = A.jg;
	bool* isUntouch = A.isUntouch;
	int i, j, k/*, ind_knot*/;
	int k_left, k_right;

	//int p; // номер узла (точки)
	int n_X = mesh.XYZ.nX;
	int n_Y = mesh.XYZ.nY;
	int n_Z = mesh.XYZ.nZ;
	//double* X = mesh.X;
	//double* Y = mesh.Y;
	//double* Z = mesh.Z;
	//Coord* XYZarray = mesh.XYZ.data;

	Area* areas = mesh.Areas.data;
	int N_ar = mesh.Areas.count;

	int* IX_w = mesh.XYZ.IXYZW.IXW;
	int* IY_w = mesh.XYZ.IXYZW.IYW;
	int* IZ_w = mesh.XYZ.IXYZW.IZW;

	int ind;
	int n_X1 = n_X - 1;

	int r, l, m, d/*, b_i, b_j*/;
	//int block_ind;
	int i_beg, i_end, j_beg, j_end, k_beg, k_end;
	double lambda_r, /*hi_r, delta_r,*/ gamma_r;
	for (r = 0; r < N_ar; r++)
	{
		areaCur = areas[r];
		i_beg = IX_w[areaCur.lXW];
		i_end = IX_w[areaCur.rXW];
		j_beg = IY_w[areaCur.lYW];
		j_end = IY_w[areaCur.rYW];
		k_beg = IZ_w[areaCur.lZW];
		k_end = IZ_w[areaCur.rZW];
		lambda_r = areaCur.lambda;
		gamma_r = areaCur.gamma;
		for (k = k_beg; k < k_end; k++) // по высоте
		{
			for (j = j_beg; j < j_end; j++) // по вертикали
			{
				for (i = i_beg; i < i_end; i++) // по горизонтали
				{
					
					fillPointsOfElem(loc, i, j, k, mesh.XYZ);
					buildLocalG(G_local, loc);
					buildLocalM(M_local, loc);
					buildLocalB(M_local, localF, localB, loc);
					int c = 0;

					for (d = 0, c = 0; d < N1; d++)
						for (m = 0; m < N1; m++)
							for (l = 0; l < N1; l++, c++)
								L[c] = getGlobeNum(i + l, j + m, k + d, n_X, n_Y, n_Z);

					for (m = 0; m < N; m++)
					{
						di[L[m]] += lambda_r * G_local[m][m] + gamma_r*M_local[m][m];
						isUntouch[L[m]] = false;
						b[L[m]] += localB[m];
					}

					for (m = 1; m < N; m++)
					{
						k_left = ig[L[m]];
						for (l = 0; l < m; l++)
						{
							k_right = ig[L[m] + 1];
							while (jg[k_left] != L[l])
							{
								ind = (k_left + k_right) / 2; // djpvj;yj
								if (jg[ind] <= L[l])
								{
									k_left = ind;
								}
								else
								{
									k_right = ind;
								}
							}

							gg[k_left] += lambda_r * G_local[m][l] + gamma_r * M_local[m][l];
							k_left++;
						}
					}
				}
			}
		}
	}

	for (i = 0; i < n; i++)
		if (isUntouch[i]) di[i] = 1;

	return 0;
}

void initVector(double* v, int n, double mean)
{
	for (int i = 0; i < n; i++)
		v[i] = mean;
}

int addFirstConditions(matrix& A, double* b, CalculationArea& mesh)
{
	int n = A.n;
	int n_gg = A.n_gg;
	double* di = A.di;
	double* gg = A.gg;
	int* ig = A.ig;
	int* jg = A.jg;

	int i, ind, i_ind, ind_beg, ind_end;
	double cond1Mean;

	int i_x1, i_x2, j_y1, j_y2, k_z1, k_z2;
	int i_x, j_y, k_z;
	int ind_current;

	int* IX_w = mesh.XYZ.IXYZW.IXW;
	int* IY_w = mesh.XYZ.IXYZW.IYW;
	int* IZ_w = mesh.XYZ.IXYZW.IZW;

	BorderFirstCond* borderFirstCondAr = mesh.borders.firstCondStor.data;

	int nFirstCond = mesh.borders.firstCondStor.count;
	Coord* points = mesh.XYZ.data;
	int b_ind;
	/*
	double* X = mesh.X;
	double* Y = mesh.Y;
	double* Z = mesh.Z;

	int n_X = mesh.nX;
	int n_Y = mesh.nY;
	int n_Z = mesh.nZ;
	*/
	for (i = 0; i < nFirstCond; i++)
	{
		i_x1 = IX_w[borderFirstCondAr[i].horizL];
		i_x2 = IX_w[borderFirstCondAr[i].horizR];

		j_y1 = IY_w[borderFirstCondAr[i].vertL];
		j_y2 = IY_w[borderFirstCondAr[i].vertR];

		k_z1 = IZ_w[borderFirstCondAr[i].heightL];
		k_z2 = IZ_w[borderFirstCondAr[i].heightR];

		i_x2++;
		j_y2++;
		k_z2++;

		for (k_z = k_z1; k_z < k_z2; k_z++)
		{
			for (j_y = j_y1; j_y < j_y2; j_y++)
			{
				for (i_x = i_x1; i_x < i_x2; i_x++)
				{
					ind = mesh.XYZ.getKnotIndex(i_x, j_y, k_z);
					//ind = getGlobeNum(i_x, j_y, k_z, n_X, n_Y, n_Z);
					cond1Mean = funct_u1(points[ind]);
					ind_beg = ig[ind];
					ind_end = ig[ind + 1];
					for (i_ind = ind_beg; i_ind < ind_end; i_ind++)
					{
						b_ind = jg[i_ind];

						b[b_ind] -= gg[i_ind] * cond1Mean;
						gg[i_ind] = 0;
					}

					ind_current = ind;
					for (ind++; ind < n; ind++)
					{
						ind_beg = ig[ind];
						ind_end = ig[ind + 1];
						for (i_ind = ind_beg; i_ind < ind_end; i_ind++)
						{
							if (jg[i_ind] == ind_current)
							{

								b[ind] -= gg[i_ind] * cond1Mean;
								gg[i_ind] = 0;
							}
						}
					}
					di[ind_current] = 1;
					b[ind_current] = cond1Mean;
				}
			}
		}
	}
	return 0;
}
int addSecondConditions(matrix& A, double* b, CalculationArea& mesh)
{
	int i;
	double M2[N2][N2];

	int horL, horR, verL, verR, heiL, heiR;
	int iHor, jVer, kHei;

	int* IX_w = mesh.XYZ.IXYZW.IXW;
	int* IY_w = mesh.XYZ.IXYZW.IYW;
	int* IZ_w = mesh.XYZ.IXYZW.IZW;

	BorderSecondCond* borderSecondCondAr = mesh.borders.secondCondStor.data;

	int nSecondCond = mesh.borders.secondCondStor.count;
	/*
	double* X = mesh.X;
	double* Y = mesh.Y;
	double* Z = mesh.Z;

	int n_X = mesh.nX;
	int n_Y = mesh.nY;
	int n_Z = mesh.nZ;
	*/
	double cord1[N2], cord2[N2];
	int pIndex[N2];
	double localEtta[N2];
	double addLocalB[N2];

	Coord* points = mesh.XYZ.data;

	int k, l;
	for (i = 0; i < nSecondCond; i++)
	{
		horL = IX_w[borderSecondCondAr[i].horizL];
		horR = IX_w[borderSecondCondAr[i].horizR];

		verL = IY_w[borderSecondCondAr[i].vertL];
		verR = IY_w[borderSecondCondAr[i].vertR];

		heiL = IZ_w[borderSecondCondAr[i].heightL];
		heiR = IZ_w[borderSecondCondAr[i].heightR];

		if (horL == horR)
		{
			horR++;
			verR++;
			heiR++;
			iHor = horL;

			for (kHei = heiL; kHei < heiR; kHei++)
			{
				for (jVer = verL; jVer < verR; jVer++)
				{
					pIndex[0] = mesh.XYZ.getKnotIndex(iHor, jVer, kHei);
					pIndex[1] = mesh.XYZ.getKnotIndex(iHor, jVer + 1, kHei);
					pIndex[2] = mesh.XYZ.getKnotIndex(iHor, jVer, kHei + 1);
					pIndex[3] = mesh.XYZ.getKnotIndex(iHor, jVer + 1, kHei + 1);

					/*
					pIndex[0] = getGlobeNum(iHor, jVer, kHei, n_X, n_Y, n_Z);
					pIndex[1] = getGlobeNum(iHor, jVer + 1, kHei, n_X, n_Y, n_Z);
					pIndex[2] = getGlobeNum(iHor, jVer, kHei + 1, n_X, n_Y, n_Z);
					pIndex[3] = getGlobeNum(iHor, jVer + 1, kHei + 1, n_X, n_Y, n_Z);
					*/

					cord1[0] = points[pIndex[0]].y;
					cord1[1] = points[pIndex[1]].y;
					cord1[2] = points[pIndex[2]].y;
					cord1[3] = points[pIndex[3]].y;

					cord2[0] = points[pIndex[0]].z;
					cord2[1] = points[pIndex[1]].z;
					cord2[2] = points[pIndex[2]].z;
					cord2[3] = points[pIndex[3]].z;

					for (k = 0; k < N2; k++)
						localEtta[k] = funct_du_dn(points[pIndex[k]]);

					buildLocalMFor2Cond(M2, cord1, cord2);
					multiplyMatrix4(M2, localEtta, addLocalB);

					for (k = 0; k < N2; k++)
						b[pIndex[k]] += addLocalB[k];
				}
			}
		}

		if (verL == verR)
		{
			horR++;
			verR++;
			heiR++;
			jVer = verL;

			for (kHei = heiL; kHei < heiR; kHei++)
			{
				for (iHor = horL; iHor < horR; iHor++)
				{
					pIndex[0] = mesh.XYZ.getKnotIndex(iHor, jVer, kHei);
					pIndex[1] = mesh.XYZ.getKnotIndex(iHor + 1, jVer, kHei);
					pIndex[2] = mesh.XYZ.getKnotIndex(iHor, jVer, kHei + 1);
					pIndex[3] = mesh.XYZ.getKnotIndex(iHor + 1, jVer + 1, kHei + 1);
					/*
					pIndex[0] = getGlobeNum(iHor, jVer, kHei, n_X, n_Y, n_Z);
					pIndex[1] = getGlobeNum(iHor + 1, jVer, kHei, n_X, n_Y, n_Z);
					pIndex[2] = getGlobeNum(iHor, jVer, kHei + 1, n_X, n_Y, n_Z);
					pIndex[3] = getGlobeNum(iHor + 1, jVer, kHei + 1, n_X, n_Y, n_Z);
					
					/*
					cord1[0] = X[pIndex[0]];
					cord1[1] = X[pIndex[1]];
					cord1[2] = X[pIndex[2]];
					cord1[3] = X[pIndex[3]];

					cord2[0] = Z[pIndex[0]];
					cord2[1] = Z[pIndex[1]];
					cord2[2] = Z[pIndex[2]];
					cord2[3] = Z[pIndex[3]];
					*/

					cord1[0] = points[pIndex[0]].x;
					cord1[1] = points[pIndex[1]].x;
					cord1[2] = points[pIndex[2]].x;
					cord1[3] = points[pIndex[3]].x;

					cord2[0] = points[pIndex[0]].z;
					cord2[1] = points[pIndex[1]].z;
					cord2[2] = points[pIndex[2]].z;
					cord2[3] = points[pIndex[3]].z;

					for (k = 0; k < N2; k++)
						localEtta[k] = funct_du_dn(points[pIndex[k]]);

					buildLocalMFor2Cond(M2, cord1, cord2);
					multiplyMatrix4(M2, localEtta, addLocalB);

					for (k = 0; k < N2; k++)
						b[pIndex[k]] += addLocalB[k];
				}
			}
		}

		if (heiL == heiR)
		{
			horR++;
			verR++;
			heiR++;
			kHei = heiL;

			for (jVer = verL; jVer < verR; jVer++)
			{
				for (iHor = horL; iHor < horR; iHor++)
				{
					pIndex[0] = mesh.XYZ.getKnotIndex(iHor, jVer, kHei);
					pIndex[1] = mesh.XYZ.getKnotIndex(iHor + 1, jVer, kHei);
					pIndex[2] = mesh.XYZ.getKnotIndex(iHor, jVer + 1, kHei);
					pIndex[3] = mesh.XYZ.getKnotIndex(iHor + 1, jVer + 1, kHei);

					/*
					pIndex[0] = getGlobeNum(iHor, jVer, kHei, n_X, n_Y, n_Z);
					pIndex[1] = getGlobeNum(iHor + 1, jVer, kHei, n_X, n_Y, n_Z);
					pIndex[2] = getGlobeNum(iHor, jVer + 1, kHei, n_X, n_Y, n_Z);
					pIndex[3] = getGlobeNum(iHor + 1, jVer + 1, kHei, n_X, n_Y, n_Z);

					cord1[0] = X[pIndex[0]];
					cord1[1] = X[pIndex[1]];
					cord1[2] = X[pIndex[2]];
					cord1[3] = X[pIndex[3]];

					cord2[0] = Y[pIndex[0]];
					cord2[1] = Y[pIndex[1]];
					cord2[2] = Y[pIndex[2]];
					cord2[3] = Y[pIndex[3]];
					*/

					cord1[0] = points[pIndex[0]].x;
					cord1[1] = points[pIndex[1]].x;
					cord1[2] = points[pIndex[2]].x;
					cord1[3] = points[pIndex[3]].x;

					cord2[0] = points[pIndex[0]].y;
					cord2[1] = points[pIndex[1]].y;
					cord2[2] = points[pIndex[2]].y;
					cord2[3] = points[pIndex[3]].y;

					for (k = 0; k < N2; k++)
						localEtta[k] = funct_du_dn(points[pIndex[k]]);

					buildLocalMFor2Cond(M2, cord1, cord2);
					multiplyMatrix4(M2, localEtta, addLocalB);

					for (k = 0; k < N2; k++)
						b[pIndex[k]] += addLocalB[k];
				}
			}
		}
	}

	return 0;
}

int addThirdConditions(matrix& A, double* b, CalculationArea& mesh)
{
	/*
	int n = A.n;
	int n_gg = A.n_gg;
	double* di = A.di;
	double* gg = A.gg;
	int* ig = A.ig;
	int* jg = A.jg;

	int i, j; 

	double M2[N2][N2];

	int horL, horR, verL, verR, heiL, heiR;
	int iHor, jVer, kHei;
	int ind;

	int* IX_w = mesh.IXW;
	int* IY_w = mesh.IYW;
	int* IZ_w = mesh.IZW;

	BorderThirdCond* borderThirdCondAr = mesh.borders.thirdCondStor.data;

	int nThirdCond = mesh.borders.thirdCondStor.count;

	double* X = mesh.X;
	double* Y = mesh.Y;
	double* Z = mesh.Z;

	int n_X = mesh.nX;
	int n_Y = mesh.nY;
	int n_Z = mesh.nZ;

	double cord1[N2], cord2[N2];
	int pIndex[N2];
	double localUb[N2];
	double addLocalB[N2];
	double betta;

	int m;

	int k_left, k_right;
	int k, l;
	for (i = 0; i < nThirdCond; i++)
	{
			horL = IX_w[borderThirdCondAr[i].horizL];
			horR = IX_w[borderThirdCondAr[i].horizR];

			verL = IY_w[borderThirdCondAr[i].vertL];
			verR = IY_w[borderThirdCondAr[i].vertR];

			heiL = IZ_w[borderThirdCondAr[i].heightL];
			heiR = IZ_w[borderThirdCondAr[i].heightR];

			betta = borderThirdCondAr[i].betta;

			if (horL == horR)
			{
				//horR++;
				//verR++;
				//heiR++;
				iHor = horL;

				for (kHei = heiL; kHei < heiR; kHei++)
				{
					for (jVer = verL; jVer < verR; jVer++)
					{
						pIndex[0] = getGlobeNum(iHor, jVer, kHei, n_X, n_Y, n_Z);
						pIndex[1] = getGlobeNum(iHor, jVer + 1, kHei, n_X, n_Y, n_Z);
						pIndex[2] = getGlobeNum(iHor, jVer, kHei + 1, n_X, n_Y, n_Z);
						pIndex[3] = getGlobeNum(iHor, jVer + 1, kHei + 1, n_X, n_Y, n_Z);

						cord1[0] = Y[pIndex[0]];
						cord1[1] = Y[pIndex[1]];
						cord1[2] = Y[pIndex[2]];
						cord1[3] = Y[pIndex[3]];

						cord2[0] = Z[pIndex[0]];
						cord2[1] = Z[pIndex[1]];
						cord2[2] = Z[pIndex[2]];
						cord2[3] = Z[pIndex[3]];

						for (k = 0; k < N2; k++)
							localUb[k] = funct_u_b(X[pIndex[k]], Y[pIndex[k]], Z[pIndex[k]]);

						buildLocalMFor2Cond(M2, cord1, cord2);
						multiplyMatrix4(M2, localUb, addLocalB);


						for (k = 0; k < N2; k++)
						{
							di[pIndex[k]] += betta * M2[k][k];
							b[pIndex[k]] += addLocalB[k];
						}
						
						for (m = 1; m < N2; m++)
						{
							k_left = ig[pIndex[m]];
							for (l = 0; l < m; l++)
							{
								k_right = ig[pIndex[m] + 1];
								while (jg[k_left] != pIndex[l])
								{
									ind = (k_left + k_right) / 2; // djpvj;yj
									if (jg[ind] <= pIndex[l])
									{
										k_left = ind;
									}
									else
									{
										k_right = ind;
									}
								}

								gg[k_left] += betta * M2[m][l];;
								k_left++;
							}
						}
						
					}
				}
			}

			if (verL == verR)
			{
				//horR++;
				//verR++;
				//heiR++;
				jVer = verL;

				for (kHei = heiL; kHei < heiR; kHei++)
				{
					for (iHor = horL; iHor < horR; iHor++)
					{
						pIndex[0] = getGlobeNum(iHor, jVer, kHei, n_X, n_Y, n_Z);
						pIndex[1] = getGlobeNum(iHor + 1, jVer, kHei, n_X, n_Y, n_Z);
						pIndex[2] = getGlobeNum(iHor, jVer, kHei + 1, n_X, n_Y, n_Z);
						pIndex[3] = getGlobeNum(iHor + 1, jVer, kHei + 1, n_X, n_Y, n_Z);

						cord1[0] = X[pIndex[0]];
						cord1[1] = X[pIndex[1]];
						cord1[2] = X[pIndex[2]];
						cord1[3] = X[pIndex[3]];

						cord2[0] = Z[pIndex[0]];
						cord2[1] = Z[pIndex[1]];
						cord2[2] = Z[pIndex[2]];
						cord2[3] = Z[pIndex[3]];

						for (k = 0; k < N2; k++)
							localUb[k] = funct_u_b(X[pIndex[k]], Y[pIndex[k]], Z[pIndex[k]]);

						buildLocalMFor2Cond(M2, cord1, cord2);
						multiplyMatrix4(M2, localUb, addLocalB);


						for (k = 0; k < N2; k++)
						{
							di[pIndex[k]] += betta * M2[k][k];
							b[pIndex[k]] += addLocalB[k];
						}

						for (k = 1; k < N2; k++)
						{
							k_left = ig[pIndex[k]];
							for (l = 0; l < k; l++)
							{
								k_right = ig[pIndex[k] + 1];
								while (jg[k_left] != pIndex[l])
								{
									ind = (k_left + k_right) / 2; // djpvj;yj
									if (jg[ind] <= pIndex[l])
									{
										k_left = ind;
									}
									else
									{
										k_right = ind;
									}
								}

								gg[k_left] += betta * M2[k][l];
								k_left++;
							}
						}
					}
				}
			}

			if (heiL == heiR)
			{
				//horR++;
				//verR++;
				//heiR++;
				kHei = heiL;

				for (jVer = verL; jVer < verR; jVer++)
				{
					for (iHor = horL; iHor < horR; iHor++)
					{
						pIndex[0] = getGlobeNum(iHor, jVer, kHei, n_X, n_Y, n_Z);
						pIndex[1] = getGlobeNum(iHor + 1, jVer, kHei, n_X, n_Y, n_Z);
						pIndex[2] = getGlobeNum(iHor, jVer + 1, kHei, n_X, n_Y, n_Z);
						pIndex[3] = getGlobeNum(iHor + 1, jVer + 1, kHei, n_X, n_Y, n_Z);

						cord1[0] = X[pIndex[0]];
						cord1[1] = X[pIndex[1]];
						cord1[2] = X[pIndex[2]];
						cord1[3] = X[pIndex[3]];

						cord2[0] = Y[pIndex[0]];
						cord2[1] = Y[pIndex[1]];
						cord2[2] = Y[pIndex[2]];
						cord2[3] = Y[pIndex[3]];

						for (k = 0; k < N2; k++)
							localUb[k] = funct_u_b(X[pIndex[k]], Y[pIndex[k]], Z[pIndex[k]]);

						buildLocalMFor2Cond(M2, cord1, cord2);
						multiplyMatrix4(M2, localUb, addLocalB);


						for (k = 0; k < N2; k++)
						{
							di[pIndex[k]] += betta * M2[k][k];
							b[pIndex[k]] += addLocalB[k];
						}

						for (k = 1; k < N2; k++)
						{
							k_left = ig[pIndex[k]];
							for (l = 0; l < k; l++)
							{
								k_right = ig[pIndex[k] + 1];
								while (jg[k_left] != pIndex[l])
								{
									ind = (k_left + k_right) / 2; // djpvj;yj
									if (jg[ind] <= pIndex[l])
									{
										k_left = ind;
									}
									else
									{
										k_right = ind;
									}
								}

								gg[k_left] += betta * M2[k][l];
								k_left++;
							}
						}
					}
				}
			}

			
		}
*/
	return 0;
}

void printVector(double* q, int n, string name)
{
	cout << name << ": " << endl;
	cout << "\t";
	for (int i = 0; i < n; i++)
		cout << q[i] << " ";
	cout << endl;
}

void printSolution(double* q, int nX, int nY, int nZ)
{
	int i, j, k, c;
	cout << "\tSolution:" << endl;
	for (k = 0, c = 0; k < nZ; k++)
	{
		for (j = 0; j < nY; j++)
		{
			for (i = 0; i < nX; i++, c++)
				cout << setw(8) << q[c] << " ";

			cout << " ";
		}

		cout << endl;
	}
}

void printErr(double* q, Coord* points, int nX, int nY, int nZ)
{
	int i, j, k, c;
	cout << "\tError:" << endl;
	for (k = 0, c = 0; k < nZ; k++)
	{
		for (j = 0; j < nY; j++)
		{
			for (i = 0; i < nX; i++, c++)
				cout << q[c] - funct_u1(points[c]) << " ";

			cout << " ";
		}

		cout << endl;
	}
}

void writeSolInVect(double* q, Coord* points, int nX, int nY, int nZ)
{
	int i, j, k, c;;
	for (k = 0, c = 0; k < nZ; k++)
	{
		for (j = 0; j < nY; j++)
		{
			for (i = 0; i < nX; i++, c++)
				q[c] = funct_u1(points[c]);
		}
	}
}

void writeComparison(double* x, double* y, int n, string name)
{
	cout << "\t" << name << endl;
	for (int i = 0; i < n; i++)
	{
		cout << i << "  " << x[i] << "  " << y[i] << "  " << x[i] - y[i] << endl;
	}
}

int buildLocalMFor2Cond(double M2[N2][N2], double cord1[N2], double cord2[N2])
{
	double a0 = (cord1[1] - cord1[0]) * (cord2[2] - cord2[0]) - (cord2[1] - cord2[0]) * (cord1[2] - cord1[0]);
	double a1 = (cord1[1] - cord1[0]) * (cord2[3] - cord2[2]) - (cord2[1] - cord2[0]) * (cord1[3] - cord1[2]);
	double a2 = (cord1[3] - cord1[1]) * (cord2[2] - cord2[0]) - (cord2[3] - cord2[1]) * (cord1[2] - cord1[0]);

	float sign = a0 / abs(a0);
	int i, j;
	for (i = 0; i < N2; i++)
	{
		for (j = 0; j < N2; j++)
		{
			M2[i][j] = sign * (a0 * M2_1[i][j] + (a1 * M2_2[i][j] + a2 * M2_3[i][j]) / 2) / 36;
		}
	}		

	return 0;
}