#include "balance.h"
#include "array.h"
#include "programlog.h"
#include <set>

Balance::Balance()
{
	nFaces = 0;
	betta = dQ = diBuf = bNew = NULL;
	epsBalance = 1e-10;
	alpha = 1;
}
void Balance::init(CalculationArea calculationArea)
{
	faceStore = calculationArea.faceStore;
	flowStore = calculationArea.flowStore;
	finitElementStore = calculationArea.finitElementStore;
	changeMemory(faceStore.count, finitElementStore.nFinitElement);
	buildPortrait();
}

bool Balance::allocateMemory(unsigned int n)
{
	nFaces = n;
	if (n == 0) return false;
	betta = new double[n];
	dQ = new double[n];
	diBuf = new double[n];
	bNew = new double[n];
	return true;
}

bool Balance::changeMemory(unsigned int nFaces, unsigned int nElems)
{
	if (this->nElems != nElems)
	{
		delete[] betta;
		betta = nElems == 0 ? NULL : new double[nElems];
		this->nElems = nElems;
	}
	if (this->nFaces != nFaces)
	{
		delete[] dQ;
		delete[] diBuf;
		delete[] bNew;
		if (nFaces == 0)
			dQ = diBuf = bNew = NULL;
		else
		{
			dQ = new double[nFaces];
			diBuf = new double[nFaces];
			bNew = new double[nFaces];
		}
		this->nFaces = nFaces;
	}
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
	double bettaElem;
	arrayspace::fill_vec(b, nFaces, 0);
	arrayspace::fill_vec(di, nFaces, 0);
	for (iElem = 0; iElem < nElems; iElem++)
	{
		iFinElemFaces = finitElements[iElem].faces;
		signs = finitElements[iElem].flowSign;
		
		sumQ = 0;
		bettaElem = betta[iElem];
		for (i = 0; i < FACES_NUM; i++)
		{
			iFace = iFinElemFaces[i];
			for (j = i + 1; j < FACES_NUM; j++)
			{
				jFace = iFinElemFaces[j];
				//ind = slae.A.getElemIndG(iFace, jFace);
				slae.A.setElem(iFace, jFace, bettaElem * signs[i] * signs[j]);
			}
			sumQ += signs[i] * flows[iFace];
			di[iFace] += betta[iElem];
		}

		for (i = 0; i < FACES_NUM; i++)
		{
			iFace = iFinElemFaces[i];
			b[iFace] -= bettaElem*signs[i] * sumQ;
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
	double maxFlow;
	double epsCur = 0;
	for (iElem = 0; iElem < nElems; iElem++)
	{
		iFinElemFaces = finitElements[iElem].faces;
		signs = finitElements[iElem].flowSign;
		sumQ = 0;
		maxFlow = 0;
		for (i = 0; i < FACES_NUM; i++)
		{
			iFace = iFinElemFaces[i];
			sumQ += signs[i]*(flows[iFace] + dQ[iFace]);
			maxFlow = maxFlow < abs(flows[iFace] + dQ[iFace]) ? abs(flows[iFace] + dQ[iFace]) : maxFlow;
		}
		betta[iElem] = 1 / (abs(maxFlow));
		epsCur += betta[iElem] * abs(sumQ);
		
	}
	return epsCur;
}

double Balance::calcDisbalanse()
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
	//double maxFlow;
	double epsCur = 0;
	for (iElem = 0; iElem < nElems; iElem++)
	{
		iFinElemFaces = finitElements[iElem].faces;
		signs = finitElements[iElem].flowSign;
		sumQ = 0;
		//maxFlow = 0;
		for (i = 0; i < FACES_NUM; i++)
		{
			iFace = iFinElemFaces[i];
			sumQ += signs[i] * (flows[iFace] + dQ[iFace]);
			//maxFlow = maxFlow < abs(flows[iFace] + dQ[iFace]) ? abs(flows[iFace] + dQ[iFace]) : maxFlow;
		}
		//betta[iElem] = 1 / (abs(maxFlow));
		epsCur += betta[iElem] * abs(sumQ);

	}
	return epsCur;
}

double Balance::functMin(double mean)
{
	double* errB = bNew;
	arrayspace::plus(diBuf, slae.A.di, 1.0/mean, slae.A.n);
	add2Cond();
	los.solve(slae.A, slae.b, dQ, 10000, 1e-14);

	//slae.A.mult(dQ, bNew);
	//arrayspace::minus(slae.b, bNew, errB, slae.A.n);
	return abs(reculcBetta());
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
	arrayspace::fill_vec(dQ, nFaces, 0);
	arrayspace::fill_vec(betta, nElems, 1);
	alpha = 1;
	while (reculcBetta() > epsBalance)
	{
		//arrayspace::fill_vec(betta, nFaces, 1);
		fillSlae();
		//slae.A.printFullMatrix();
		//for (int i = 0; i < nFaces; i++)
		//	std::cout << slae.b[i] << std::endl;
		arrayspace::copy(diBuf, slae.A.di, nFaces);
		alpha = findAlpha();
		arrayspace::plus(diBuf, slae.A.di, alpha, nFaces);
		add2Cond();
		los.solve(slae.A, slae.b, dQ, 10000, 1e-14);
		//std::cout << "Balance Solut" << std::endl;
		//for (int i = 0; i < nFaces; i++)
		//	std::cout << dQ[i] << std::endl;
	}

	int* iFinElemFaces = finitElementStore.finitElements[0].faces;
	signed char* signs = finitElementStore.finitElements[0].flowSign;
	double sumQ = 0;

	arrayspace::plus(flowStore.flows, dQ, flowStore.flows, nFaces);

	//for (int i = 0; i < FACES_NUM; i++)
	//{
	//	int iFace = iFinElemFaces[i];
	//	sumQ += signs[i] * (flowStore.flows[iFace] + dQ[iFace]);
	//}
	refreshFlowSigns();
}

void Balance::add2Cond()
{
	int n2Faces = faces2CondStore.nFaces;
	int* iFaces = faces2CondStore.iFaces;

	for (int i = 0; i < n2Faces; i++)
		slae.setOneVariableSolve(iFaces[i], 0);
}

void Balance::refreshFlowSigns()
{
	int iElem;
	int iLocalFace;
	int iGlobalFace;
	int iNegbFinElem;
	int* currVerFace;
	int	NeighbVerFace[VER_NUM_FACE];
	int k;
	int nElems = finitElementStore.nFinitElement;
	FinitElement iFinitElement;
	FinitElement* finitElements = finitElementStore.finitElements;
	double* flows = flowStore.flows;
	Face* faces = faceStore.faces;
	for (iElem = 0; iElem < nElems; iElem++)
	{
		iFinitElement = finitElements[iElem];
		for (iLocalFace = 0; iLocalFace < FACES_NUM; iLocalFace++)
		{
			iGlobalFace = iFinitElement.faces[iLocalFace];
			if (flows[iGlobalFace] < 0)
			{
				finitElements[iElem].flowSign[iLocalFace] *= -1;
				flows[iGlobalFace] = abs(flows[iGlobalFace]);
				iNegbFinElem = faceStore.findNeighboringFinitElem(iElem, iFinitElement, iLocalFace);
				if (iNegbFinElem > -1)
				{
					currVerFace = faces[iGlobalFace].knots;
					for (k = 0; k < FACES_NUM; k++)
					{
						finitElements[iNegbFinElem].getFaceGlobalNum(k, NeighbVerFace);
						if (arrayspace::isSameWithoutOrdUniqe(NeighbVerFace, currVerFace, VER_NUM_FACE))
							finitElements[iNegbFinElem].flowSign[k] *= -1;
					}
				}
			}
		}
	}
}