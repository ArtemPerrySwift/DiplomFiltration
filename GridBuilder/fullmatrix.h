#pragma once
const int MAT_DIM_2 = 2;
const int MAT_DIM_3 = 3;
const int MAT_DIM_4 = 4;
const int MAT_DIM_8 = 8;

void multMatrix3x2(double mat1[MAT_DIM_3][MAT_DIM_2], double mat2[MAT_DIM_2][MAT_DIM_3], double res[MAT_DIM_3][MAT_DIM_3]);

void multMatrix3x2(double mat1[MAT_DIM_3][MAT_DIM_2], double mat2[MAT_DIM_3][MAT_DIM_2], double res[MAT_DIM_3][MAT_DIM_3]);

void multMatrixVect(double mat[MAT_DIM_8][MAT_DIM_8], double vect[MAT_DIM_8], double res[MAT_DIM_8]);

void multMatrixVect(double mat[MAT_DIM_3][MAT_DIM_2], double vect[MAT_DIM_2], double res[MAT_DIM_3]);

void multMatrixVect(double mat[MAT_DIM_4][MAT_DIM_4], double vect[MAT_DIM_4], double res[MAT_DIM_4]);

double countDetMatrix3(double J[3][3]);
