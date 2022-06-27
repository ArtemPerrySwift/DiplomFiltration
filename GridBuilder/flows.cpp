#include "flows.h"
#include "array.h"
#include "facebasic.h"
#include "DifferentEquParams.h"
#include "flowelemculcer.h"
#include <iostream>
#include <iomanip>

FlowCulcer::FlowCulcer()
{
	qN = 0;
	nFlows = 0;
	q = new double[qN];
	flows = NULL;
}

void FlowCulcer::init(CalculationArea calculationArea, double* q)
{
	this->calculationArea = calculationArea;
	int n = calculationArea.coordsStore.count;
	nFlows = calculationArea.flowStore.count;
	flows = calculationArea.flowStore.flows;
	outPhaseStorage = calculationArea.borderFacesStore.phaseStorage;
	/*
	if (qN != n)
	{
		delete[] this->q;
		this->q = new double[n];
		qN = n;
	}

	arrayspace::copy(this->q, q, n);
	*/
	this->q = q;

}

bool FlowCulcer::isInnerN(int iLocalFace)
{
	switch (iLocalFace)
	{
	case 0:
		return true;
	case 1:
		return false;
	case 2:
		return true;
	case 3:
		return true;
	case 4:
		return false;
	case 5:
		return false;
	}
}

void FlowCulcer::calcFlows()
{
	arrayspace::fill_vec(flows, calculationArea.flowStore.count, 0);
	FinitElement* finitElements = calculationArea.finitElementStore.finitElements;
	int nElems = calculationArea.finitElementStore.nFinitElement;
	Face* faces = calculationArea.faceStore.faces;
	int i, j, k;
	int iFace;
	FaceBasicFlow faceBasicFlow;
	FinitElement finitElement;
	double flow;

	//double flowSigns[FACES_NUM];
	double flowLoc[FACES_NUM];
	FlowElemCulcer flowElemCulcer;
	int iNeibElem;
	int faceVer[VER_NUM_FACE];
	int faceVerCur[VER_NUM_FACE];
	bool fl;
	for (i = 0; i < nElems; i++)
	{
		finitElement = finitElements[i];
		flowElemCulcer.init(finitElement, calculationArea.coordsStore, q);
		flowElemCulcer.getFlows(flowLoc, finitElements[i].flowSign);
		/*
		for (j = 0; j < FACES_NUM; j++)
		{
			iFace = finitElement.faces[j];
			if (flows[iFace] * flowLoc[j] < 0)
			{
				flows[iFace] -= flowLoc[j] / 2.0;
			}
			else
			{
				iNeibElem = calculationArea.faceStore.findNeighboringFinitElem(i, finitElements[i], j);
				for (k = 0; k < FACES_NUM; k++)
				{
					finitElements[iNeibElem].getFaceGlobalNum(k, faceVer);
					if (arrayspace::isSameWithoutOrdUniqe(finitElement.faces, faceVer, VER_NUM_FACE))
					{
						if (abs(flows[iFace]) > abs(flowLoc[j] / 2))
						{
							finitElements[iNeibElem].flowSign[k] = finitElements[i].flowSign[j];
							finitElements[i].flowSign[j] *= -1;
							flows[iFace] -= flowLoc[j] / 2;

						}
						else
						{
							finitElements[iNeibElem].flowSign[k] = -finitElements[i].flowSign[j];
							flows[iFace] = flowLoc[j] / 2 - flows[iFace];
						}
					}
						
				}				
			}
		}
		*/
		for (j = 0; j < FACES_NUM; j++)
		{
			iFace = finitElement.faces[j];

			iNeibElem = calculationArea.faceStore.findNeighboringFinitElem(i, finitElements[i], j);
			if (iNeibElem != -1)
			{
				if (finitElements[i].flowSign[j] > 0)
				{
					flowLoc[j] *= DifferentEquParams::lambda(finitElement.K, finitElement.phaseStorage.phases, finitElement.phaseStorage.count);
				}
				else
				{
					flowLoc[j] *= DifferentEquParams::lambda(finitElements[iNeibElem].K, finitElements[iNeibElem].phaseStorage.phases, finitElements[iNeibElem].phaseStorage.count);
				}

				if (iNeibElem > i)
				{
					flows[iFace] = flowLoc[j] / 2;
				}
				else
				{
					fl = true;
					finitElement.getFaceGlobalNum(j, faceVerCur);
					for (k = 0; k < FACES_NUM && fl; k++)
					{
						finitElements[iNeibElem].getFaceGlobalNum(k, faceVer);
						if (arrayspace::isSameWithoutOrdUniqe(faceVerCur, faceVer, VER_NUM_FACE))
						{
							finitElements[iNeibElem].flowSign[k] = -finitElements[i].flowSign[j];
							fl = false;
						}
					}
					k--;
					if (finitElements[iNeibElem].flowSign[k] == finitElements[i].flowSign[j])
					{
						if (flows[iFace] > flowLoc[j] / 2)
						{
							finitElements[i].flowSign[j] = -finitElements[iNeibElem].flowSign[k];
							flows[iFace] -= flowLoc[j] / 2;
						}
						else
						{
							finitElements[iNeibElem].flowSign[k] = -finitElements[i].flowSign[j];
							flows[iFace] = flowLoc[j] / 2 - flows[iFace];
						}
					}
					else
					{
						flows[iFace] += flowLoc[j] / 2;
					}
				}
			}
			else
			{
				if (finitElements[i].flowSign[j] > 0)
				{
					flowLoc[j] *= DifferentEquParams::lambda(finitElement.K, finitElement.phaseStorage.phases, finitElement.phaseStorage.count);
				}
				else
				{
					flowLoc[j] *= DifferentEquParams::lambda(finitElement.K, outPhaseStorage.phases, outPhaseStorage.count);
				}
				flows[iFace] = flowLoc[j];
			}
				
			/*
			if (flows[iFace] * flowLoc[j] < 0)
			{
				flows[iFace] -= flowLoc[j] / 2.0;
			}
			else
			{
				if (abs(flows[iFace]) > abs(flowLoc[j] / 2))
				{
					finitElements[i].flowSign[j] *= -1;
					flows[iFace] -= flowLoc[j] / 2;
				}
				else
				{
					iNeibElem = calculationArea.faceStore.findNeighboringFinitElem(i, finitElements[i], j);
					if(iNeibElem != -1)
						for (k = 0; k < FACES_NUM; k++)
						{
							finitElements[iNeibElem].getFaceGlobalNum(k, faceVer);
							if (arrayspace::isSameWithoutOrdUniqe(finitElement.faces, faceVer, VER_NUM_FACE))
								finitElements[iNeibElem].flowSign[k] = -finitElements[i].flowSign[j];
						}
					flows[iFace] = flowLoc[j] / 2 - flows[iFace];
				}

			}
			*/
		}

		/*if (!i)
		{
			std::cout << "Flow of 0 element" << std::endl;
			for (j = 0; j < FACES_NUM; j++)
			{
				std::cout << flowLoc[j] << " ";
			}
			std::cout << std::endl;
		}*/
		/*
		for (j = 0; j < FACES_NUM; j++)
		{
			iFace = finitElement.faces[j];
			if (flows[iFace] * flowLoc[j] > 0)
			{
				flows[iFace] -= flowLoc[j];
				flows[iFace] /= 2;
				if (flows[iFace] < 0)
				{
					iNeibElem = calculationArea.faceStore.findNeighboringFinitElem(i, finitElements[i], j);
					for (k = 0; k < FACES_NUM; k++)
					{
						finitElements[iNeibElem].getFaceGlobalNum(k, faceVer);
						if (arrayspace::isSameWithoutOrdUniqe(finitElement.faces, faceVer, VER_NUM_FACE))
							finitElements[iNeibElem].flowSign[k] *= -1;
					}
				}
				else
				{
					finitElements[i].flowSign[j] *= -1;
				}
			}
			else
			{
				flows[iFace] = (flows[iFace] == 0.0) ? flowLoc[j] / 2.0 : -flowLoc[j] / 2.0;
			}
			
			//faceBasicFlow.init(faces[iFace], calculationArea.coordsStore, q, isInnerN(i), DifferentEquParams::lambda(finitElement.K, finitElement.phases.phases, finitElement.phases.count));
			//low = faceBasicFlow.getFlow();
			//
			//flowSigns[j] = flow < 0 ? -1: 1;
			
		}
		*/
		//arrayspace::copy(finitElements[i].flowSign, flowSigns, FACES_NUM);
		
	}

	for (i = 0; i < nFlows; i++)
		flows[i] = abs(flows[i]);
	
}

double* FlowCulcer::getflows() 
{
	double* flowsRes;
	flowsRes = new double[calculationArea.faceStore.count];
	return flowsRes;
}

void FlowCulcer::printFlows(std::ofstream &out)
{
	out << "____________________________________" << std::endl;
	out << "|----N of Face----|------Flow-------|" << std::endl;
	out << "____________________________________" << std::endl;
	for (int i = 0; i < nFlows; i++)
		out << "|" << setw(17) << i << "|" << setw(17) << flows[i] << "|" << std::endl;
	out << "___________________________________" << std::endl;
}