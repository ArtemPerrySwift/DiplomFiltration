#pragma once
#include <fstream>

using namespace std;
namespace arrayspace
{
	template<typename type>
	void copy(type* arReciver, type* arSource, unsigned int n);

	void copy(double* arReciver, double* arSource, unsigned int n);

	void copy(int* arReciver, int* arSource, unsigned int n);

	void copy(signed char* arReciver, signed char* arSource, unsigned int n);

	template<typename type>
	int read_vec(type* vec, int n, ifstream& in);

	int read_vec(double* vec, int n, ifstream& in);

	int read_vec(int* vec, int n, ifstream& in);

	template<typename type>
	int fill_vec(type* vec, int n, type elem);

	int fill_vec(double* vec, int n, double elem);

	int fill_vec(int* vec, int n, int elem);

	double scal(double* v1, double* v2, int n);

	void minus(double* v1, double* v2, double* res, int n);

	void mult(double* v, double a, double* res, int n);

	bool sort(int* mas, int iBeg, int n);

	int serchInSorted(int* mas, int n, int elem);

	int serchInSorted(int* mas, int iBeg, int n, int elem);

	bool isSameWithoutOrdUniqe(const int* mas1, const int* mas2, const int n);

	bool copyElemsByInd(int* masInd, double* masRes, double* masSor, int nResElems);

	void plus(double* v, double* res, double a, int n);

	void plus(double* v1, double* v2, double* res, int n);
}
