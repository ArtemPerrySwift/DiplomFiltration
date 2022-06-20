#include "sparsematrix.h"
#include <utility>
#include <iomanip>
#include <iostream>
#include "array.h"
#include <set>
using namespace arrayspace;

namespace matrix
{
	SparseMatrix::SparseMatrix()
	{
		n = 0;
		ig = NULL;
		jg = NULL;
		di = NULL;
	}

	void SparseMatrix::copyPortrait(SparseMatrix& M)
	{
		n = M.n;
		ig = M.ig;
		jg = M.jg;
	}	

	void SparseMatrix::allocateMemoryForElems()
	{
		di = new double[n];
	}

	void SparseMatrix::freeMemory()
	{
		delete[] di;
	}

	bool SparseMatrix::copyElems(SparseMatrix& M)
	{
		if (n != M.n)
			return false;
		arrayspace::copy(di, M.di, n);
		return true;
	}
	void SparseMatrix::changePortrait(SparseMatrix& M)
	{
		if (n != M.n || ig[n] != M.ig[M.n])
		{
			freeMemory();
			copyPortrait(M);
			allocateMemoryForElems();
		}
	}

	SparseMatrixAsym::SparseMatrixAsym() : SparseMatrix()
	{
		ggl = NULL;
		ggu = NULL;
	}

	void SparseMatrixAsym::allocateMemoryForElems() 
	{	
		SparseMatrix::allocateMemoryForElems();
		ggl = new double[ig[n]];
		ggu = new double[ig[n]];
	}

	bool SparseMatrixAsym::copyElems(SparseMatrixAsym& M)
	{
		if (!SparseMatrix::copyElems(M)) return false;
		if (M.ig[n] != ig[n]) return false;
		arrayspace::copy(ggl, M.ggl, n);
		arrayspace::copy(ggu, M.ggu, n);
		return true;
	}

	bool SparseMatrixAsym::copyElems(SparseMatrixSym& M)
	{
		if (!SparseMatrix::copyElems(M)) return false;
		if (M.ig[n] != ig[n]) return false;
		arrayspace::copy(ggl, M.gg, n);
		arrayspace::copy(ggu, M.gg, n);
		return true;
	}

	void SparseMatrixSym::allocateMemoryForElems()
	{
		SparseMatrix::allocateMemoryForElems();
		gg = new double[ig[n]];
	}

	void SparseMatrixAsym::copy(SparseMatrixAsym& M)
	{
		SparseMatrix::changePortrait(M);
		copyElems(M);
	}

	void SparseMatrixAsym::copy(SparseMatrixSym& M)
	{
		SparseMatrix::changePortrait(M);
		copyElems(M);
	}

	void SparseMatrixAsym::freeMemory()
	{
		SparseMatrix::freeMemory();
		delete[] ggl;
		delete[] ggu;
	}

	void SparseMatrixSym::freeMemory()
	{
		SparseMatrix::freeMemory();
		delete[] gg;
	}
	
	SparseMatrixSym::SparseMatrixSym() : SparseMatrix()
	{
		gg = NULL;
	}

	bool SparseMatrixSym::copyElems(SparseMatrixSym& M)
	{
		if (!SparseMatrix::copyElems(M)) return false;
		if (M.ig[n] != ig[n]) return false;
		arrayspace::copy(gg, M.gg, n);
		return true;
	}

	void SparseMatrixSym::copy(SparseMatrixSym& M)
	{
		SparseMatrix::changePortrait(M);
		copyElems(M);
	}

	void SparseMatrix::mult(double* v, double* res)
	{
		int i, j, ij;
		int i_beg, i_end;

		for (i = 0; i < n; i++)
		{
			res[i] = 0;
			i_beg = ig[i];
			i_end = ig[i + 1];
			for (ij = i_beg; ij < i_end; ij++)
			{
				j = jg[ij];
				multMatElemAndVecElem(i, j, ij, v, res);
			}
			res[i] += di[i] * v[i];
		}
	}
	double* Matrix:: operator * (double* v)
	{
		double* res = new double[n];
		mult(v, res);
		return res;
	}

	void SparseMatrixSym::multMatElemAndVecElem(int i, int j, int ij, double* v, double* res) 
	{
		res[i] += gg[ij] * v[j];
		res[j] += gg[ij] * v[i];
	}

	void SparseMatrixAsym::multMatElemAndVecElem(int i, int j, int ij, double* v, double* res)
	{
		res[i] += ggl[ij] * v[j];
		res[j] += ggu[ij] * v[i];
	}

	int SparseMatrixAsym::decomp_mat_LU(SparseMatrixAsym& LU)
	{
		std::cout << "Разложение Холесского для матрицы начато" << std::endl;
		double* LU_di = LU.di;
		double* LU_ggl = LU.ggl;
		double* LU_ggu = LU.ggu;
		int* LU_ig = LU.ig;
		int* LU_jg = LU.jg;
		
		int i, j;
		int i_beg, i_end, j_beg, j_end;
		int nerazl = 0;
		int ii;
		int k, c, d;
		int ik, kj, jk, ki, ij;
		//bool fl;
		double sumij, sumji, sumdi;
		
		for (i = 0; i < n; i++)
		{
			//double sum = 0;
			sumdi = 0;
			for (j = 0; j < i; j++)
			{
				//fl = true;
				for (ij = ig[i]; ij < ig[i + 1] && jg[ij] < j; ij++);

				if (jg[ij] != j) continue;

				sumij = 0;
				sumji = 0;
				for (k = 0; k < j; k++)
				{
					for (ik = ig[i]; ik < ig[i + 1] && jg[ik] < k; ik++);
					for (jk = ig[j]; jk < ig[j + 1] && jg[jk] < k; jk++);

					for (ki = ig[k]; ki < ig[k + 1] && jg[ki] < i; ki++);
					for (kj = ig[k]; kj < ig[k + 1] && jg[kj] < j; kj++);

					//for (c = ig[i]; c < ig[i + 1] && jg[c] < k; c++);
					//for (d = ig[j]; d < ig[j + 1] && jg[d] < k; d++);
					if(jg[ik] == k && jg[kj] == j)
						sumij += LU.ggl[ik] * LU.ggu[kj];

					if(jg[jk] == k && jg[ki] == i)
						sumji += LU.ggl[jk] * LU.ggu[ki];
				}

				LU.ggl[ij] = (ggl[ij] - sumij);
				LU.ggu[ij] = (ggu[ij] - sumji)/LU.di[j];
				//if (LU.ggl[ij] != 0)
				//	cout << "Not 0 " << LU.ggl[ij] << " " << i << " " << ij << endl;
				sumdi += LU.ggl[ij] * LU.ggu[ij];
			}
			//for (m = 0; m < i; m++)
				//sumdi += LU.di[m];
			//cout << "sumdi[" << i << "]" << sumdi << endl;
			//if (di[i] - sumdi < 0)
			//{
			//	cout << "Di < 0 i = " << i;
			//}
				
			LU.di[i] = di[i] - sumdi;
			//cout << "Lu.di[" << i << "]" << LU.di[i] << endl;
			if (abs(LU.di[i] / sumdi) < 1e-12)
				cout << "Warning: Diagonal element can be equal 0" << endl;
		}
		/*
		for (i = 0; i < n; i++)
		{
			int nail = ig[ii];
			int koil = ig[ii + 1];
			double sumdi = 0;
			while (nail < koil)
			{
				for (i = 0; i < koil; i++)
				{
					double suml = 0;
					double sumu = 0;
					int jj = jg[i];
					if (i != nail)
					{
						int najj = ig[jj];
						int kojj = ig[jj + 1];
					}
				}
			}
		}
		*/
		/*
		for (i = 0; i < n; i++)
		{
			i_beg = ig[i];
			i_end = ig[i + 1];
			LU_di[i] = di[i];
			for (ij = i_beg; ij < i_end; ij++)
			{
				j = jg[ij];
				j_beg = ig[j];
				j_end = ig[j + 1];
				LU_ggl[ij] = ggl[ij];
				LU_ggu[ij] = ggu[ij];
				for (ik = i_beg, kj = j_beg; kj < j_end && ik < ij;)
				{
					if (jg[ik] == jg[kj])
					{
						LU_ggl[ij] -= LU_ggl[ik] * LU_ggu[kj];
						LU_ggu[ij] -= LU_ggl[kj] * LU_ggu[ik];
						ik++;
						kj++;
					}
					else if (jg[ik] > jg[kj]) kj++;
						 else ik++;
				}
				LU_ggu[ij] /= LU_di[j];
				LU_di[i] -= LU_ggl[ij] * LU_ggu[ij];
			}
		}
		*/

		std::cout << "Разложение Холесского для матрицы успешно построено" << std::endl;
		return 0;
		
	}

	void SparseMatrixAsym::fillMatrix(double mean)
	{
		fill_vec(di, n, mean);
		fill_vec(ggl, ig[n], mean);
		fill_vec(ggu, ig[n], mean);
	}

	void SparseMatrixSym::fillMatrix(double mean)
	{
		fill_vec(di, n, mean);
		fill_vec(gg, ig[n], mean);
	}

	void SparseMatrixAsym::printFullMatrix()
	{
		int COUT_WIDTH = 10;
		cout << setw(COUT_WIDTH);
		bool isIndexFound = false; // Был ли найден соответствующий индекс в массиве jg 
		int c;
		for (int i = 0; i < n; i++)
		{
			int i_beg = ig[i];
			int i_end = ig[i + 1];
			for (c = 0; c < jg[i_beg]; c++)
				cout << setw(COUT_WIDTH) << 0.0 << " ";

			for (int j = i_beg; j < i_end; j++)
			{
				for (; c < jg[j]; c++)
					cout << setw(COUT_WIDTH) << 0.0 << " ";
				cout << setw(COUT_WIDTH) << ggu[j] << " ";
				c++;
			}
			for (; c < i; c++)
				cout << setw(COUT_WIDTH) << 0.0 << " ";

			cout << setw(COUT_WIDTH) << di[i] << " ";

			for (int j = i + 1; j < n; j++)
			{
				int j_beg = ig[j];
				int j_end = ig[j + 1];
				isIndexFound = false;
				for (int k = j_beg; k < j_end && !isIndexFound; k++)
				{
					isIndexFound = (jg[k] == i);
					if (isIndexFound)
					{
						cout << setw(COUT_WIDTH) << ggl[k] << " ";
					}
				}
				if (!isIndexFound)
					cout << setw(COUT_WIDTH) << 0.0 << " ";
			}
			cout << setw(COUT_WIDTH) << endl;
		}
	}

	void SparseMatrix::buildPortrait(CoordStorage coordStorage, FinitElementStore finitElementStore)
	{
		set<int>* map;
		map = new set<int>[coordStorage.count];
		int kuzlov = coordStorage.count;
		int ktr = finitElementStore.nFinitElement;
		FinitElement* finalElements = finitElementStore.finitElements;
		n = kuzlov;
		int indexes[VER_NUM];
		
		int i, j, k, l;
		for (k = 0; k < ktr; k++)
		{
			for (l = 0; l < VER_NUM; l++)
				indexes[l] = finalElements[k].ver[l];

			/*
			indexes[0] = finalElements[k].ver1;
			indexes[1] = finalElements[k].ver2;
			indexes[2] = finalElements[k].ver3;
			indexes[3] = finalElements[k].ver4;
			*/

			for (i = 0; i < VER_NUM; i++)
			{
				for (j = 0; j < VER_NUM; j++)
				{
					if (indexes[i] > indexes[j])
						map[indexes[i]].insert(indexes[j]);
				}
			}
		}
		ig = new int[kuzlov + 1];
		ig[0] = 0;

		for (i = 0; i < kuzlov; i++)
		{
			ig[i + 1] = ig[i] + map[i].size();
		}
		jg = new int[ig[kuzlov]];

		int ijCount;
		//int* elem;

		for (i = 0; i < kuzlov; i++)
		{
			j = ig[i];
				
			for (set<int>::iterator elem = map[i].begin(); elem != map[i].end(); elem++, j++)
				jg[j] = *elem;

				//for (auto& item : map[i])
				//{
				//	jg[j] = item;
				//	j++;
				//}
		}
		
		allocateMemoryForElems();
	}


	/**
	SparseMatrixPortrait::SparseMatrixPortrait(SparseMatrixPortrait& portarait)
	{
		int* ig = portarait.ig;
		int* jg = portarait.jg;
		int n = portarait.n;

		this->di = new double[n];
		this->n = n;
		this->ig = new int[n + 1];
		arrayspace::copy(this->ig, ig, n + 1);

		int n_jg = ig[n];
		this->jg = new int[n_jg];
		arrayspace::copy(this->jg, jg, n_jg);

		//allmemory();
		//allocateMemoryForElems();
	}
	*/

	void DecompSparseMatrixLDLt<SparseMatrixSym>::decompMatrix(SparseMatrixSym& M)
	{
		copy(M);
		std::cout << "Разложение Холесского для матрицы начато" << std::endl;
		int i, j;
		int i_beg, i_end, j_beg, j_end;
		int nerazl = 0;
		int ii;
		int k, c, d;
		int ik, kj, jk, ki, ij;
		//bool fl;
		double sumij, sumji, sumdi;

		for (i = 0; i < n; i++)
		{
			//double sum = 0;
			sumdi = 0;
			for (j = 0; j < i; j++)
			{
				//fl = true;
				for (ij = ig[i]; ij < ig[i + 1] && jg[ij] < j; ij++);

				if (jg[ij] != j) continue;

				sumij = 0;
				sumji = 0;
				for (k = 0; k < j; k++)
				{
					for (ik = ig[i]; ik < ig[i + 1] && jg[ik] < k; ik++);
					for (jk = ig[j]; jk < ig[j + 1] && jg[jk] < k; jk++);

					//for (ki = ig[k]; ki < ig[k + 1] && jg[ki] < i; ki++);
					//for (kj = ig[k]; kj < ig[k + 1] && jg[kj] < j; kj++);
					//
					//for (c = ig[i]; c < ig[i + 1] && jg[c] < k; c++);
					//for (d = ig[j]; d < ig[j + 1] && jg[d] < k; d++);
					if (jg[ik] == k && jg[kj] == j)
						sumij += gg[ik] *di[k]* gg[kj];

					//if (jg[jk] == k && jg[ki] == i)
					//	sumji += gg[jk] * gg[ki];
				}

				gg[ij] = (gg[ij] - sumij)/di[j];
				//gg[ij] = (ggu[ij] - sumji) / LU.di[j];
				//if (LU.ggl[ij] != 0)
				//	cout << "Not 0 " << LU.ggl[ij] << " " << i << " " << ij << endl;
				sumdi += gg[ij] * gg[ij] * di[j];
			}
			//for (m = 0; m < i; m++)
				//sumdi += LU.di[m];
			//cout << "sumdi[" << i << "]" << sumdi << endl;
			//if (di[i] - sumdi < 0)
			//{
			//	cout << "Di < 0 i = " << i;
			//}

			di[i] = di[i] - sumdi;
			//cout << "Lu.di[" << i << "]" << LU.di[i] << endl;
			if (abs(di[i] / sumdi) < 1e-12)
				cout << "Warning: Diagonal element can be equal 0" << endl;
		}
	}

	void DecompSparseMatrixLDLt<SparseMatrixSym>::solveLy(double* f, double* y)
	{
		int i, j;
		int i_beg, i_end;
		int ij;
		for (i = 0; i < n; i++)
		{
			y[i] = f[i];
			i_beg = ig[i];
			i_end = ig[i + 1];
			for (ij = i_beg; ij < i_end; ij++)
			{
				j = jg[ij];
				y[i] -= gg[ij] * y[j];
			}
		}
		solveDz(y, y);
	}

	void DecompSparseMatrixLDLt<SparseMatrixSym>::solveUx(double* y, double* x)
	{
		int i, j;
		int i_beg, i_end;
		int ij;
		for (i = 0; i < n; i++)
			x[i] = y[i];

		for (i = n - 1; i > -1; i--)
		{
			i_beg = ig[i];
			i_end = ig[i + 1];
			for (ij = i_beg; ij < i_end; ij++)
			{
				j = jg[ij];
				x[j] -= gg[ij] * x[i];
			}
		}
	}

	void DecompSparseMatrixLDLt<SparseMatrixSym>::solveDz(double* f, double* z)
	{
		for (int i = 0; i < n; i++)
			z[i] = f[i] / di[i];
	}
}



