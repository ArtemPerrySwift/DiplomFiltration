#include "balance.h"
#include "array.h"
#include "programlog.h"
#include <set>

Balance::Balance()
{
	nFace = 0;
	betta = dQ = diBuf = bNew = NULL;
	epsBalance = 1e-10;
	alpha = 1;
}
void Balance::init(CalculationArea calculationArea)
{
	faceStore = calculationArea.faceStore;
	flowStore = calculationArea.flowStore;
	finitElementStore = calculationArea.finitElementStore;
	changeMemory(faceStore.count);
	buildPortrait();
}

bool Balance::allocateMemory(unsigned int n)
{
	nFace = n;
	if (n == 0) return false;
	betta = new double[n];
	dQ = new double[n];
	diBuf = new double[n];
	bNew = new double[n];
	return true;
}

bool Balance::changeMemory(unsigned int n)
{
	if (n == nFace) return true;
	deleteMemory();
	allocateMemory(n);
	return true;
}

void Balance::deleteMemory()
{
	delete[] betta;
	delete[] dQ;
	delete[] diBuf;
	delete[] bNew;
	betta = dQ = diBuf = bNew = NULL;
}

void Balance::buildPortrait()
{
	std::set<int>* map;
	int nFaces = faceStore.count;
	map = new std::set<int>[nFaces];
	int ktr = finitElementStore.nFinitElement;
	FinitElement* finalElements = finitElementStore.finitElements;
	slae.A.n = nFaces;
	int indexes[FACES_NUM];

	int i, j, k, l;
	for (k = 0; k < ktr; k++)
	{
		for (l = 0; l < FACES_NUM; l++)
			indexes[l] = finalElements[k].faces[l];

		for (i = 0; i < FACES_NUM; i++)
		{
			for (j = 0; j < FACES_NUM; j++)
			{
				if (indexes[i] > indexes[j])
					map[indexes[i]].insert(indexes[j]);
			}
		}
	}
	slae.A.ig = new int[nFaces + 1];
	int* ig = slae.A.ig;
	ig[0] = 0;

	for (i = 0; i < nFaces; i++)
	{
		ig[i + 1] = ig[i] + map[i].size();
	}
	
	slae.A.jg = new int[ig[nFaces]];
	int* jg = slae.A.jg;
	int ijCount;

	for (i = 0; i < nFaces; i++)
	{
		j = ig[i];

		for (std::set<int>::iterator elem = map[i].begin(); elem != map[i].end(); elem++, j++)
			jg[j] = *elem;
	}

	slae.A.allocateMemoryForElems();
	slae.b = new double[slae.A.n];
}

void Balance::fillSlae()
{
	//reculcBetta();
	int nElems = finitElementStore.nFinitElement;
	FinitElement* finitElements = finitElementStore.finitElements;
	double* flows = flowStore.flows;
	int iElem;
	//FinitElement iFinitElement;
	int* iFinElemFaces;
	int i, j;
	int iFace, jFace;

	//int* ig = slae.A.ig;
	//int* jg = slae.A.jg;
	double* gg = slae.A.gg;
	double* di = slae.A.di;
	int ind;
	signed char* signs;
	double sumQ;
	double* b = slae.b;
	arrayspace::fill_vec(b, nFace, 0);
	for (iElem = 0; iElem < nElems; iElem++)
	{
		iFinElemFaces = finitElements[iElem].faces;
		signs = finitElements[iElem].flowSign;
		sumQ = 0;
		for (i = 0; i < FACES_NUM; i++)
		{
			iFace = iFinElemFaces[i];
			for (j = i + 1; j < FACES_NUM; j++)
			{
				jFace = iFinElemFaces[j];
				//ind = slae.A.getElemIndG(iFace, jFace);
				slae.A.setElem(iFace, jFace, betta[i] * signs[i] * signs[j]);
			}
			sumQ += flows[iFace];
		}

		for (i = 0; i < FACES_NUM; i++)
		{
			iFace = iFinElemFaces[i];
			di[iFace] = 2;
			b[iFace] -= signs[i] * sumQ;
		}
	}
}

double Balance::reculcBetta()
{
	int nElems = finitElementStore.nFinitElement;
	FinitElement* finitElements = finitElementStore.finitElements;
	double* flows = flowStore.flows;
	int iElem;
	int* iFinElemFaces;
	signed char* signs;
	double sumQ;
	int i;
	int iFace;
	double epsCur = 0;
	for (iElem = 0; iElem < nElems; iElem++)
	{
		iFinElemFaces = finitElements[iElem].faces;
		signs = finitElements[iElem].flowSign;
		sumQ = 0;
		for (i = 0; i < FACES_NUM; i++)
		{
			iFace = iFinElemFaces[i];
			sumQ += signs[i]*(flows[iFace] + dQ[iFace]);
		}
		epsCur += betta[iFace] * abs(sumQ);
		betta[iFace] = 1 / (abs(sumQ));
	}
	return epsCur;
}

double Balance::functMin(double mean)
{
	double* errB = bNew;
	arrayspace::plus(diBuf, slae.A.di, 1.0/mean, slae.A.n);
	los_prec.solve(slae.A, slae.b, dQ, 10000, 1e-14);
	slae.A.mult(dQ, bNew);
	arrayspace::minus(slae.b, bNew, errB, slae.A.n);
	return abs(arrayspace::scal(errB, errB, slae.A.n) - epsBalance);
}

double Balance::findAlpha() { return 1.0/calcMin(alpha); }

bool Balance::setEpsBalanse(double eps)
{
	if (eps <= 0)
	{
		programlog::writeErr("Balanse epsilon is less or equal 0");
		return false;
	}
		
	epsBalance = eps;
	return true;
}

void Balance::balanceFlows()
{
	arrayspace::fill_vec(dQ, nFace, 0);
	arrayspace::fill_vec(betta, nFace, 1);
	while (reculcBetta() > epsBalance)
	{
		fillSlae();
		arrayspace::copy(diBuf, slae.A.di, nFace);
		alpha = findAlpha();
		arrayspace::plus(diBuf, slae.A.di, alpha, nFace);
		los_prec.solve(slae.A, slae.b, dQ, 10000, 1e-14);
	}
	int* iFinElemFaces = finitElementStore.finitElements[0].faces;
	signed char* signs = finitElementStore.finitElements[0].flowSign;
	double sumQ = 0;

	arrayspace::plus(flowStore.flows, dQ, flowStore.flows, nFace);

	for (int i = 0; i < FACES_NUM; i++)
	{
		int iFace = iFinElemFaces[i];
		sumQ += signs[i] * (flowStore.flows[iFace] + dQ[iFace]);
	}

}