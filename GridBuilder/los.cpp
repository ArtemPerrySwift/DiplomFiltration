#include "los.h"
#include "array.h";
LOS::LOS()
{
	n = 0;
	r = z = p = buf_v = f = NULL;
}

bool LOS::allocateMemory(unsigned int n)
{
	this->n = n;
	if (n == 0) return false;
	r = new double[n];
	z = new double[n];
	p = new double[n];
	f = new double[n];
	buf_v = new double[n];
	return true;
}

bool LOS::changeMemory(unsigned int n)
{
	if (n == this->n) return true;
	deleteMemory();
	allocateMemory(n);
	return true;
}

void LOS::deleteMemory()
{
	delete[] r;
	delete[] z;
	delete[] p;
	delete[] buf_v;
	r = z = p = buf_v = NULL;
}

int LOS::solve(Matrix &A, double* b, double* x, int maxiter, double eps)
{
	int n = A.n;
	changeMemory(n);
	int i, j, k;
	arrayspace::copy(f, b, n);
	f = b;

	//cout << endl;
	double err;
	double begErr;
	double norm_ff = sqrt(arrayspace::scal(f, f, n));
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

		double a_k, b_k;
		double skal_pp;

		err = sqrt(arrayspace::scal(r, r, n)) / norm_ff;
		begErr = abs(err - eps);
		fl = i < maxiter&& err > eps;
		for (k = 0; fl && k < maxiter; i++, k++)
		{
			skal_pp = arrayspace::scal(p, p, n);
			a_k = arrayspace::scal(p, r, n) / skal_pp;
			for (j = 0; j < n; j++)
			{
				x[j] += a_k * z[j];
				r[j] -= a_k * p[j];
			}
			A.mult(r, buf_v);
			b_k = -arrayspace::scal(p, buf_v, n) / skal_pp;
			for (j = 0; j < n; j++)
			{
				z[j] = r[j] + b_k * z[j];
				p[j] = buf_v[j] + b_k * p[j];
			}
			err = sqrt(arrayspace::scal(r, r, n)) / norm_ff;
			fl = i < maxiter&& err > eps;
			//progressBar.showProgress(((abs(err - eps)/ begErr* abs(err - eps) / begErr)*(1 - eps/err) + eps / err* eps / err) * 100.0);
			//cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
		}
		//itarations = i;
		//nev = err;
	}
	//progressBar.end();
	//std::cout << std::endl;
	//cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
	if (i == maxiter) return 1;
	return 0;
}

int LOS_precond::solve(SparseMatrixSym& A, double* b, double* x, int maxiter, double eps)
{
	changeMemory(A.n);
	LLt.decompMatrix(A);
	arrayspace::copy(f, b, n);
	f = b;
	int n = A.n;
	int i, j, k;
	double norm_ff = sqrt(arrayspace::scal(f, f, n));
	double err;	
	cout << "normff " << norm_ff << endl;
	bool fl = true;
	i = 0;
	while (fl)
	{
		A.mult(x, buf_v);
		for (j = 0; j < n; j++)
			buf_v[j] = f[j] - buf_v[j];

		cout << "arrayspace::scal(buf_v, buf_v) = " << arrayspace::scal(buf_v, buf_v, n) << "   ";

		LLt.solveLy(buf_v, r);
		cout << "arrayspace::scal(r, r) = " << arrayspace::scal(r, r, n) << "   ";

		LLt.solveUx(r, z);

		A.mult(z, buf_v);
		LLt.solveLy(buf_v, p);

		double a_k, b_k;
		double skal_pp;

		err = sqrt(arrayspace::scal(r, r, n)) / norm_ff;
		fl = i < maxiter&& err > eps;
		for (k = 0; fl && k < 1000; i++, k++)
		{

			skal_pp = arrayspace::scal(p, p, n);
			a_k = arrayspace::scal(p, r, n) / skal_pp;

			for (j = 0; j < n; j++)
			{
				x[j] += a_k * z[j];
				r[j] -= a_k * p[j];
			}
			//cout << "arrayspace::scal(x, x) = " << arrayspace::scal(x, x, n) << "   ";
			//cout << "arrayspace::scal(r, x) = " << arrayspace::scal(x, x, n) << "   ";
			LLt.solveUx(r, buf_v);
			//calc_Ux(LU, buf_v, r, is_U_diag_1);

			A.mult(buf_v, buf_v1);
			LLt.solveLy(buf_v1, buf_v);
			//calc_Lx(LU, buf_v, buf_v1, is_L_diag_1);

			b_k = -arrayspace::scal(p, buf_v, n) / skal_pp;
			LLt.solveUx(r, buf_v1);
			//calc_Ux(LU, buf_v1, r, is_U_diag_1);
			for (j = 0; j < n; j++)
			{
				z[j] = buf_v1[j] + b_k * z[j];
				p[j] = buf_v[j] + b_k * p[j];
			}
			err = sqrt(arrayspace::scal(r, r, n)) / norm_ff;
			fl = i < maxiter&& err > eps;
			//cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
			//itarations = i;
			//nev = err;
		}
		cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
	}
	return 0;
}

void LOS_precond::deleteMemory()
{
	LOS::deleteMemory();
	delete[] buf_v1;
}

bool LOS_precond::allocateMemory(unsigned int n)
{
	if (!LOS::allocateMemory(n)) return false;
	buf_v1 = new double[n];
}