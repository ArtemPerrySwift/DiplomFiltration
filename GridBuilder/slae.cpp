#include <cmath>
#include <iostream>
#include <iomanip>
#include "programlog.h"
#include "slae.h"
#include "array.h"

using namespace arrayspace;
namespace slae
{

	//void SLAE<MatrixType>::init(MatrixType A)
	//{
	//	b = new double[A.n];
	//	init(A, b);
	//}
	void SLAE<SparseMatrixAsym>::setOneVariableSolve(int iVar, double varMean)
	{
		int* ig = A.ig;
		int* jg = A.jg;
		double* ggl = A.ggl;
		double* ggu = A.ggu;
		double* di = A.di;
		int n = A.n;

		int ind_beg = ig[iVar];
		int ind_end = ig[iVar + 1];
		int i_ind, b_ind;

		int ind_current;
		for (i_ind = ind_beg; i_ind < ind_end; i_ind++)
		{
			b_ind = jg[i_ind];

			b[b_ind] -= ggu[i_ind] * varMean;

			ggu[i_ind] = 0;
			ggl[i_ind] = 0;
		}

		ind_current = iVar;
		for (iVar++; iVar < n; iVar++)
		{
			ind_beg = ig[iVar];
			ind_end = ig[iVar + 1];
			for (i_ind = ind_beg; i_ind < ind_end; i_ind++)
			{
				if (jg[i_ind] == ind_current)
				{

					b[iVar] -= ggl[i_ind] * varMean;
					ggu[i_ind] = 0;
					ggl[i_ind] = 0;
				}
			}
		}
		di[ind_current] = 1;
		b[ind_current] = varMean;
	}

	void SLAE<SparseMatrixSym>::setOneVariableSolve(int iVar, double varMean)
	{
		int* ig = A.ig;
		int* jg = A.jg;
		double* gg = A.gg;
		double* di = A.di;
		int n = A.n;

		int ind_beg = ig[iVar];
		int ind_end = ig[iVar + 1];
		int i_ind, b_ind;

		int ind_current;
		for (i_ind = ind_beg; i_ind < ind_end; i_ind++)
		{
			b_ind = jg[i_ind];
			b[b_ind] -= gg[i_ind] * varMean;
			gg[i_ind] = 0;
		}

		ind_current = iVar;
		for (iVar++; iVar < n; iVar++)
		{
			ind_beg = ig[iVar];
			ind_end = ig[iVar + 1];
			for (i_ind = ind_beg; i_ind < ind_end; i_ind++)
			{
				if (jg[i_ind] == ind_current)
				{
					b[iVar] -= gg[i_ind] * varMean;
					gg[i_ind] = 0;
				}
			}
		}
		di[ind_current] = 1;
		b[ind_current] = varMean;
	}

	/*
	void SLAE<MatrixType>::init(SparseMatrixAsym A, double* b)
	{
		r = new double[A.n];
		z = new double[A.n];
		p = new double[A.n];
		//f = new double[A.n];
		buf_v = new double[A.n];
		buf_v1 = new double[A.n];
		this->A = A;
		this->b = b;
		LU.copy(A);
	}
	*/
	/*
	int SLAE::count_LOS_Simple(real* x, int maxiter, real eps)
	{
		int n = A.n;
		int i, j, k;
		f = b;
		//cout << endl;
		real err;
		real begErr;
		real norm_ff = sqrt(scal(f, f, n));
		bool fl = true;
		i = 0;
		//programlog::ProgressBar progressBar;
		//progressBar.begin();
		while (fl)
		{
			//cout << endl << "Begin again" << endl;
			A.mult(x, r);
			for (j = 0; j < n; j++)
			{
				r[j] = f[j] - r[j];
				z[j] = r[j];
			}
			A.mult(z, p);

			real a_k, b_k;
			real skal_pp;

			err = sqrt(scal(r, r, n)) / norm_ff;
			begErr = abs(err - eps);
			fl = i < maxiter&& err > eps;
			for (k = 0; fl && k < maxiter; i++, k++)
			{
				skal_pp = scal(p, p, n);
				a_k = scal(p, r, n) / skal_pp;
				for (j = 0; j < n; j++)
				{
					x[j] += a_k * z[j];
					r[j] -= a_k * p[j];
				}
				A.mult(r, buf_v);
				b_k = -scal(p, buf_v, n) / skal_pp;
				for (j = 0; j < n; j++)
				{
					z[j] = r[j] + b_k * z[j];
					p[j] = buf_v[j] + b_k * p[j];
				}
				err = sqrt(scal(r, r, n)) / norm_ff;
				fl = i < maxiter&& err > eps;
				//progressBar.showProgress(((abs(err - eps)/ begErr* abs(err - eps) / begErr)*(1 - eps/err) + eps / err* eps / err) * 100.0);
				//cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
			}
			//itarations = i;
			//nev = err;
		}
		//progressBar.end();
		std::cout << std::endl;
		cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
		if (i == maxiter) return 1;
		return 0;
	}
	
	*/
	/*
	int SLAE::count_LOS(real* x, int maxiter, real eps)
	{
		A.decomp_mat_LU(LU);
		bool is_L_diag_1 = false;
		bool is_U_diag_1 = false;
		f = b;
		int n = A.n;
		int i, j, k;
		real norm_ff = sqrt(scal(f, f, n));
		real err;
		*/
		/*
		cout << "b" << endl;
		for (int i = 0; i < n; i++)
		{
			if(f[i] != 0.0)
			cout << i << " " << f[i] << " ";
		}
			
		cout << endl;
		*/
	/*
		cout << "normff " << norm_ff << endl;
		bool fl = true;
		i = 0;
		while (fl)
		{
			A.mult(x, buf_v);
			for (j = 0; j < n; j++)
				buf_v[j] = f[j] - buf_v[j];

			cout << "scal(buf_v, buf_v) = " << scal(buf_v, buf_v, n) << "   ";

			calc_Lx(LU, r, buf_v, is_L_diag_1);
			cout << "scal(r, r) = " << scal(r, r, n) << "   ";

			calc_Ux(LU, z, r, is_U_diag_1);

			A.mult(z, buf_v);
			calc_Lx(LU, p, buf_v, is_L_diag_1);

			real a_k, b_k;
			real skal_pp;

			err = sqrt(scal(r, r, n)) / norm_ff;
			fl = i < maxiter&& err > eps;
			for (k = 0; fl && k < 1000; i++, k++)
			{

				skal_pp = scal(p, p, n);
				a_k = scal(p, r, n) / skal_pp;

				for (j = 0; j < n; j++)
				{
					x[j] += a_k * z[j];
					r[j] -= a_k * p[j];
				}
				//cout << "scal(x, x) = " << scal(x, x, n) << "   ";
				//cout << "scal(r, x) = " << scal(x, x, n) << "   ";
				calc_Ux(LU, buf_v, r, is_U_diag_1);

				A.mult(buf_v, buf_v1);
				calc_Lx(LU, buf_v, buf_v1, is_L_diag_1);

				b_k = -scal(p, buf_v, n) / skal_pp;
				calc_Ux(LU, buf_v1, r, is_U_diag_1);
				for (j = 0; j < n; j++)
				{
					z[j] = buf_v1[j] + b_k * z[j];
					p[j] = buf_v[j] + b_k * p[j];
				}
				err = sqrt(scal(r, r, n)) / norm_ff;
				fl = i < maxiter&& err > eps;
				cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
				//itarations = i;
				//nev = err;
			}
			cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
		}
		return 0;
		
	}
	*/
/*
	int SLAE::calc_Ux(SparseMatrixAsym& LU, real* x, real* f, bool is_diag_1)
	{
		int n = LU.n;

		real* U_di = LU.di;
		real* U_ggu = LU.ggu;
		int* U_ig = LU.ig;
		int* U_jg = LU.jg;

		int i, j;
		int i_beg, i_end;
		int ij;
		for (i = 0; i < n; i++)
			x[i] = f[i];

		if (is_diag_1)
			for (i = n - 1; i > -1; i--)
			{
				i_beg = U_ig[i];
				i_end = U_ig[i + 1];
				for (ij = i_beg; ij < i_end; ij++)
				{
					j = U_jg[ij];
					x[j] -= U_ggu[ij] * x[i];
				}
			}
		else
			for (i = n - 1; i > -1; i--)
			{
				i_beg = U_ig[i];
				i_end = U_ig[i + 1];
				x[i] /= U_di[i];
				for (ij = i_beg; ij < i_end; ij++)
				{
					j = U_jg[ij];
					x[j] -= U_ggu[ij] * x[i];
				}
			}

		return 0;
	}

	int SLAE::calc_Lx(SparseMatrixAsym& LU, real* x, real* f, bool is_diag_1)
	{
		int n = LU.n;

		real* L_di = LU.di;
		real* L_ggl = LU.ggl;
		int* L_ig = LU.ig;
		int* L_jg = LU.jg;

		int i, j;
		int i_beg, i_end;
		int ij;
		if (is_diag_1)
			for (i = 0; i < n; i++)
			{
				x[i] = f[i];
				i_beg = L_ig[i];
				i_end = L_ig[i + 1];
				for (ij = i_beg; ij < i_end; ij++)
				{
					j = L_jg[ij];
					x[i] -= L_ggl[ij] * x[j];
				}
			}
		else
			for (i = 0; i < n; i++)
			{
				x[i] = f[i];
				i_beg = L_ig[i];
				i_end = L_ig[i + 1];
				for (ij = i_beg; ij < i_end; ij++)
				{
					j = L_jg[ij];
					x[i] -= L_ggl[ij] * x[j];
				}
				if (abs(L_di[i]) < 1e-14)
					cout << "di[" << i << "] == 0" << endl;
				x[i] /= L_di[i];
				//cout << L_di[i] << endl;
			}
		return 0;
	}

	void SLAE::setOneVariableSolve(int iVar, double varMean)
	{
		int* ig = A.ig;
		int* jg = A.jg;
		double* ggl = A.ggl;
		double* ggu = A.ggu;
		double* di = A.di;
		int n = A.n;

		int ind_beg = ig[iVar];
		int ind_end = ig[iVar + 1];
		int i_ind, b_ind;

		int ind_current;
		for (i_ind = ind_beg; i_ind < ind_end; i_ind++)
		{
			b_ind = jg[i_ind];

			b[b_ind] -= ggu[i_ind] * varMean;

			ggu[i_ind] = 0;
			ggl[i_ind] = 0;
		}

		ind_current = iVar;
		for (iVar++; iVar < n; iVar++)
		{
			ind_beg = ig[iVar];
			ind_end = ig[iVar + 1];
			for (i_ind = ind_beg; i_ind < ind_end; i_ind++)
			{
				if (jg[i_ind] == ind_current)
				{

					b[iVar] -= ggl[i_ind] * varMean;
					ggu[i_ind] = 0;
					ggl[i_ind] = 0;
				}
			}
		}
		di[ind_current] = 1;
		b[ind_current] = varMean;
	}
*/
	int solveSLAU3(double A[3][3], double b[3], double x[3])
	{
		double BufM[3][3];
		double bufb[3];
		int i, j, k;
		double buf;
		bool f;
		for (i = 0; i < 3; i++)
		{
			bufb[i] = b[i];
			for (j = 0; j < 3; j++)
				BufM[i][j] = A[i][j];
		}

		double a;
		double aModule;
		int jMax;
		for (i = 0; i < 3; i++)
		{
			a = BufM[i][i];
			aModule = abs(a);
			jMax = i;
			for (j = i + 1; j < 3; j++)
			{
				if (abs(BufM[j][i]) > aModule) jMax = j;
			}

			if (!BufM[jMax][i]) return -1;

			if (jMax != i)
			{
				for (k = i; k < 3; k++)
				{
					buf = BufM[i][k];
					BufM[i][k] = BufM[jMax][k];
					BufM[jMax][k] = buf;
				}
				buf = bufb[i];
				bufb[i] = bufb[jMax];
				bufb[jMax] = buf;
			}

			a = BufM[i][i];
			for (j = i; j < 3; j++)
				BufM[i][j] /= a;

			bufb[i] /= a;

			for (j = i + 1; j < 3; j++)
			{
				a = BufM[j][i];
				for (k = i; k < 3; k++)
				{
					BufM[j][k] -= a * BufM[i][k];
				}
				bufb[j] -= a * bufb[i];
			}
		}

		for (i = 2; i > -1; i--)
		{
			for (j = i + 1; j < 3; j++)
				bufb[i] -= BufM[i][j] * x[j];
			x[i] = bufb[i];
		}

		return 0;
	}

}
