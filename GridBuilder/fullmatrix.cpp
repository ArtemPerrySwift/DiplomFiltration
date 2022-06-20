#include "fullmatrix.h"

void multMatrix3x2(double mat1[MAT_DIM_3][MAT_DIM_2], double mat2[MAT_DIM_2][MAT_DIM_3], double res[MAT_DIM_3][MAT_DIM_3])
{
	int i, j, k;
	for (i = 0; i < MAT_DIM_3; i++)
	{
		for (j = 0; j < MAT_DIM_3; j++)
		{
			res[i][j] = 0;
			for (k = 0; k < MAT_DIM_2; k++)
				res[i][j] += mat1[i][k] * mat2[k][j];
		}
	}
}

void multMatrix3x2(double mat1[MAT_DIM_3][MAT_DIM_2], double mat2[MAT_DIM_3][MAT_DIM_2], double res[MAT_DIM_3][MAT_DIM_3])
{
	int i, j, k;
	for (i = 0; i < MAT_DIM_3; i++)
	{
		for (j = 0; j < MAT_DIM_3; j++)
		{
			res[i][j] = 0;
			for (k = 0; k < MAT_DIM_2; k++)
				res[i][j] += mat1[i][k] * mat2[j][k];
		}
	}
}

void multMatrixVect(double mat[MAT_DIM_3][MAT_DIM_2], double vect[MAT_DIM_2], double res[MAT_DIM_3])
{
	int i, k;
	for (i = 0; i < MAT_DIM_3; i++)
	{
		res[i] = 0;
		for (k = 0; k < MAT_DIM_2; k++)
			res[i] += mat[i][k] * vect[k];
	}
}

void multMatrixVect(double mat[MAT_DIM_4][MAT_DIM_4], double vect[MAT_DIM_4], double res[MAT_DIM_4])
{
	int i, k;
	for (i = 0; i < MAT_DIM_4; i++)
	{
		res[i] = 0;
		for (k = 0; k < MAT_DIM_4; k++)
			res[i] += mat[i][k] * vect[k];
	}
}

void multMatrixVect(double mat[MAT_DIM_8][MAT_DIM_8], double vect[MAT_DIM_8], double res[MAT_DIM_8])
{
	int i, k;
	for (i = 0; i < MAT_DIM_8; i++)
	{
		res[i] = 0;
		for (k = 0; k < MAT_DIM_8; k++)
			res[i] += mat[i][k] * vect[k];
	}
}

double countDetMatrix3(double J[3][3])
{
	return J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] + J[0][2] * J[1][0] * J[2][1] - J[0][2] * J[1][1] * J[2][0] - J[0][1] * J[1][0] * J[2][2] - J[0][0] * J[1][2] * J[2][1];
}