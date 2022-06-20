#include "array.h"
#include <cstring>
#include <fstream>

using namespace std;
namespace arrayspace
{
	template<typename type>
	void copy(type* arReciver, type* arSource, unsigned int n)
	{
		for (int i = 0; i < n; i++)
			arReciver[i] = arSource[i];
		//size_t generalArrraySize = n * sizeof(type); // Размер массивов в байтах
		//memcpy(arReciver, arSource, generalArrraySize); // Копируем сразу весь массив целиком
	}

	void copy(double* arReciver, double* arSource, unsigned int n)
	{
		for (int i = 0; i < n; i++)
			arReciver[i] = arSource[i];
		//size_t generalArrraySize = n * sizeof(double); // Размер массивов в байтах
		//memcpy(arReciver, arSource, generalArrraySize); // Копируем сразу весь массив целиком
	}

	void copy(int* arReciver, int* arSource, unsigned int n)
	{
		for (int i = 0; i < n; i++)
			arReciver[i] = arSource[i];
		//size_t generalArrraySize = n * sizeof(int); // Размер массивов в байтах
		//memcpy(arReciver, arSource, generalArrraySize); // Копируем сразу весь массив целиком
	}

	template<typename type>
	int read_vec(type* vec, int n, ifstream& in)
	{
		for (int i = 0; i < n; i++)
			in >> vec[i];
		return 0;
	}

	int read_vec(double* vec, int n, ifstream& in)
	{
		for (int i = 0; i < n; i++)
			in >> vec[i];
		return 0;
	}

	int read_vec(int* vec, int n, ifstream& in)
	{
		for (int i = 0; i < n; i++)
			in >> vec[i];
		return 0;
	}

	template<typename type>
	int fill_vec(type* vec, int n, type elem)
	{
		for (int i = 0; i < n; i++)
			vec[i] = elem;
		return 0;
	}

	int fill_vec(double* vec, int n, double elem)
	{
		for (int i = 0; i < n; i++)
			vec[i] = elem;
		return 0;
	}

	int fill_vec(int* vec, int n, int elem)
	{
		for (int i = 0; i < n; i++)
			vec[i] = elem;
		return 0;
	}

	double scal(double* v1, double* v2, int n)
	{
		double sum = 0;
		for (int i = 0; i < n; i++)
		{
			sum += v1[i] * v2[i];
		}
			
		return sum;
	}

	void minus(double* v1, double* v2, double* res, int n)
	{
		for (int i = 0; i < n; i++)
			res[i] = v1[i] - v2[i];

	}

	void mult(double* v, double a, double* res, int n)
	{
		for (int i = 0; i < n; i++)
			res[i] = a*v[i];
	}
	
	bool sort(int* mas, int iBeg, int n)
	{
		int tmp = 0;
		int iEnd = iBeg + n;
		for (int i = iBeg; i < iEnd; i++) {
			for (int j = iEnd - 1; j >= (i + 1); j--) {
				if (mas[j] < mas[j - 1]) {
					tmp = mas[j];
					mas[j] = mas[j - 1];
					mas[j - 1] = tmp;
				}
			}
		}
		return true;
	}

	int serchInSorted(int* mas, int n, int elem)
	{
		int i;
		for (i = 0; i < n && elem < mas[i]; i++);

		if (i == n) return -1;
		if (elem != mas[i]) return -2;

		return i;
	}

	int serchInSorted(int* mas, int iBeg, int n, int elem)
	{
		int i;
		int iEnd = iBeg + n;
		for (i = iBeg; i < iEnd && elem > mas[i]; i++);

		if (i == iEnd) return -1;
		if (elem != mas[i]) return -2;

		return i;
	}

	bool isSameWithoutOrdUniqe(int* mas1, int* mas2, int n)
	{
		int i, j;
		int elem;
		bool fl = false;
		for (i = 0; i < n && !fl; i++)
		{
			elem = mas1[i];
			fl = true;
			for (j = 0; j < n && fl; j++)
			{
				//if (i == j) continue;
				fl = elem != mas2[i];
			}
		}

		return !fl;
	}

	void copy(signed char* arReciver, signed char* arSource, unsigned int n)
	{
		for (int i = 0; i < n; i++)
			arReciver[i] = arSource[i];
	}

	bool copyElemsByInd(int* masInd, double* masRes, double* masSor, int nResElems)
	{
		for (int i = 0; i < nResElems; i++)
			masRes[i] = masSor[masInd[i]];

		return true;
	}

	void plus(double* v, double* res, double a, int n)
	{
		for (int i = 0; i < n; i++)
			res[i] = v[i] + a;
	}
}