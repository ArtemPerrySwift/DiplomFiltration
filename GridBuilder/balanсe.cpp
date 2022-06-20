#include "balance.h"
#include <set>

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
	int* jg = slae.A.jg;
	slae.A.jg = new int[ig[nFaces]];

	int ijCount;

	for (i = 0; i < nFaces; i++)
	{
		j = ig[i];

		for (std::set<int>::iterator elem = map[i].begin(); elem != map[i].end(); elem++, j++)
			jg[j] = *elem;
	}

	slae.A.allocateMemoryForElems();
}

void Balance::fillSlae()
{
	reculcBetta();
	int nElems = finitElementStore.nFinitElement;
	FinitElement* finitElements = finitElementStore.finitElements;
	int iElem;
	//FinitElement iFinitElement;
	int* iFinElemFaces;
	int i, j;
	int iFace, jFace;

	int* ig = slae.A.ig;
	int* jg = slae.A.jg;
	double* gg = slae.A.gg;
	int iBeg, iEnd;
	for (iElem = 0; iElem < nElems; iElem++)
	{
		iFinElemFaces = finitElements[iElem].faces;
		for (i = 0; i < FACES_NUM; i++)
		{
			iFace = iFinElemFaces[i];
			for (j = 0; j < FACES_NUM; j++)
			{
				jFace = iFinElemFaces[j];

			}
		}
	}
}