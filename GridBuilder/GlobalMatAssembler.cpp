#include "GlobalMatAssembler.h"
#include <iostream>
#include <cstdlib>
#include <conio.h>
#include "programlog.h"

void GlobalMatAssembler::assembleGlobalMatrix(SLAE<SparseMatrixSym>& slae, CalculationArea& mesh, double* q)
{
	FinitElement* finalElements = mesh.finitElementStore.finitElements;
	int ktr = mesh.finitElementStore.nFinitElement;
	int materialInd;
	FinalElemBase finalElemBase;
	//Material material;
	slae.A.fillMatrix(0);
	
	char x = -37;
	int del = ktr / 35;
	//std::cout << del;
	//std::cout << setlocale(LC_ALL, NULL) << std::endl;
	//setlocale(LC_ALL, "C");
	//std::cout << setlocale(LC_ALL, NULL) << std::endl;
	//std::cout << std::endl;
	programlog::ProgressBar progressBar;
	progressBar.begin();
	for (int i = 0; i < ktr; i++)
	{
		progressBar.showProgress(double(i + 1) / ktr * 100.0);
		int j, k;
		//materialInd = mesh.finalElementMaterialStorage.finalElemMaterials[i];
		//material = mesh.materialStorage.findMaterialByNum(materialInd);
		finalElemBase.setFinalElement(finalElements[i], mesh.coordsStore);
		//std::cout << "\t Построение локальной матрицы жёсткости" << std::endl;
		finalElemBase.buildLocalG(localG);
		
		/*
		std::cout << "Local G" << std::endl;
		double sum;

		for (j = 0; j < 8; j++)
		{
			sum = 0;
			for (k = 0; k < 8; k++)
			{
				sum += localG[j][k];
				std::cout << localG[j][k] << " ";
			}
			std::cout << sum << std::endl;
		}
		*/
		//std::cout << "\t Построение локального вектора b" << std::endl;
		finalElemBase.buildLocalB(localB);
		finalElemBase.buildLocalM(localM, 0);

		/*
		std::cout << "Local M" << std::endl;
		for (j = 0; j < 8; j++)
		{
			sum = 0;
			for (k = 0; k < 8; k++)
			{
				sum += localM[j][k];
				std::cout << localM[j][k] << " ";
			}
			std::cout << sum << std::endl;
		}

		std::cout << "Local B" << std::endl;
		for (j = 0; j < 8; j++)
		{
			std::cout << localB[j] << " ";
		}
		*/
		addLocalToGlobalAll(slae, localG, localM, localB, finalElemBase.globalInd);
	}
	progressBar.end();
	std::cout << std::endl;
	//getchar();
}

void GlobalMatAssembler::addLocalToGlobalAll(SLAE<SparseMatrixSym>& slae, double G_local[N_KNOTS_3D][N_KNOTS_3D], double M_local[N_KNOTS_3D][N_KNOTS_3D], double localB[N_KNOTS_3D], int L[N_KNOTS_3D])
{
	int m, l;
	int ind;
	//int k_left, k_right;
	double* di = slae.A.di;
	double* gg = slae.A.gg;
	//int* ig = slae.A.ig;
	//int* jg = slae.A.jg;
	double* b = slae.b;
	for (m = 0; m < N_KNOTS_3D; m++)
	{
		di[L[m]] += G_local[m][m] + M_local[m][m];
		b[L[m]] += localB[m];
	}
	int s, d;
	for (m = 1; m < N_KNOTS_3D; m++)
	{
		for (l = 0; l < m; l++)
		{
			//s = L[m] > L[l] ? m : l;
			//d = L[m] > L[l] ? l : m;
			/*
			//if(L[m] > L[l]) s = 
			k_left = ig[L[s]];
			k_right = ig[L[s] + 1];
			while (jg[k_left] != L[d])
			{
				ind = (k_left + k_right) / 2; // djpvj;yj
				if (jg[ind] <= L[d])
				{
					k_left = ind;
				}
				else
				{
					k_right = ind;
				}
			}

			ggl[k_left] += G_local[s][d] + M_local[s][d];
			//ggu[k_left] += G_local[d][s] + M_local[d][s];
			//k_left++;
			*/
			ind = slae.A.getElemIndG(L[m], L[l]);
			gg[ind] += G_local[m][l] + M_local[m][l];
		}
	}
}

/*
void NewtAssembler::assembleLocal(FinalElemBase finalElemBase)
{
	finalElemBase.buildNewt(localG, localB);
	finalElemBase.buildLocalM(localM, 0);
}
*/
void SimpleAssembler::assembleLocal(FinalElemBase finalElemBase)
{
	finalElemBase.buildLocalG(localG);
	finalElemBase.buildLocalM(localM, 0);
	finalElemBase.buildLocalB(localB);
}