#include "gridbuilder.h"
#include <iostream>
#include <iomanip>
#include "programlog.h"
#include "array.h"
#include "qadregu.h"
#include "fullmatrix.h"
//#include "FinalElemetsBasic.h"

//---------------Border-------------//
bool Border::readData(std::ifstream& in)
{
	in >> fun_num >> horizL >> horizR >> vertL >> vertR >> heightL >> heightR;
	/*
	in.read((char*)&(fun_num), sizeof(int)); if (in.gcount() != sizeof(int)) return false;
	in.read((char*)&(horizL), sizeof(int)); if (in.gcount() != sizeof(int)) return false;
	in.read((char*)&(horizR), sizeof(int)); if (in.gcount() != sizeof(int)) return false;
	in.read((char*)&(vertL), sizeof(int)); if (in.gcount() != sizeof(int)) return false;
	in.read((char*)&(vertR), sizeof(int)); if (in.gcount() != sizeof(int)) return false;
	in.read((char*)&(heightL), sizeof(int)); if (in.gcount() != sizeof(int)) return false;
	in.read((char*)&(heightR), sizeof(int)); if (in.gcount() != sizeof(int)) return false;
	*/
	return true;
}

//---------------Border3C-------------//
bool BorderEmptyC::readData(std::ifstream& in)
{
	Border::readData(in);
	in >> sat1;
	in >> sat2;
	//in.read((char*)&(betta), sizeof(double)); if (in.gcount() != sizeof(double)) return false;
	return true;
}

//---------------Border3C-------------//
bool Border3C::readData(std::ifstream& in)
{
	Border::readData(in);
	in >> betta;
	//in.read((char*)&(betta), sizeof(double)); if (in.gcount() != sizeof(double)) return false;
	return true;
}

//---------------BorderStorage-------------//
bool BorderStorage::allocateMemory()
{
	borders = new Border[count];
	return true;
}

bool BorderStorage::readData(std::ifstream& in, int i)
{
	Border borderBuf;
	if (!borders[i].readData(in)) return false;
	return true;
}

//---------------Border1Storage-------------//
bool Border1Storage::allocateMemory()
{
	borders = new Border1C[count];
	return true;
}

bool Border1Storage::readData(std::ifstream& in, int i)
{
	Border borderBuf;
	if (!borders[i].readData(in)) return false;
	return true;
}

//---------------Border2Storage-------------//
bool Border2Storage::allocateMemory()
{
	borders = new Border2C[count];
	return true;
}

bool Border2Storage::readData(std::ifstream& in, int i)
{
	Border borderBuf;
	if (!borders[i].readData(in)) return false;
	return true;
}

//---------------Border3Storage-------------//
bool Border3Storage::allocateMemory()
{
	borders = new Border3C[count];
	return true;
}

bool Border3Storage::readData(std::ifstream& in, int i)
{
	Border borderBuf;
	if (!borders[i].readData(in)) return false;
	return true;
}

//---------------BorderEmptyStorage-------------//
bool BorderEmptyStorage::allocateMemory()
{
	borders = new BorderEmptyC[count];
	return true;
}

bool BorderEmptyStorage::readData(std::ifstream& in, int i)
{
	Border borderBuf;
	if (!borders[i].readData(in)) return false;
	return true;
}

//---------------ContainerBorders-------------//
bool ContainerBorders::readBordersCond1(std::ifstream& in, std::string fileName) { return firstCondStor.read(in, fileName) && firstCondStorWell.read(in, fileName); }

bool ContainerBorders::readBordersCond2(std::ifstream& in, std::string fileName) { return secondCondStor.read(in, fileName) && secondCondStorWell.read(in, fileName); }

bool ContainerBorders::readBordersCond3(std::ifstream& in, std::string fileName) { return thirdCondStor.read(in, fileName) && thirdCondStorWell.read(in, fileName); }

bool ContainerBorders::readBordersCondEmpty(std::ifstream& in, std::string fileName) { return emptyCondStor.read(in, fileName); }

bool ContainerBorders::readBorders(std::ifstream& in, std::string fileName) { return  readBordersCond1(in, fileName) && readBordersCond2(in, fileName) && readBordersCond3(in, fileName) && readBordersCondEmpty(in, fileName); }

bool ContainerBorders::readBorders(std::string fileName)
{
	std::ifstream in;
	in.open(fileName);
	bool res = readBorders(in, fileName);
	in.close();
	return res;
}


//---------------AreaStorage-------------//

bool AreaStorage::read(std::ifstream& in, std::string path, int count)
{
	std::string read_file_err = "Error of reading file " + path;

	if (!in.is_open())
	{
		programlog::writeErr(read_file_err + " - Cannot open the file");
		return false;
	}

	if (count == -1)
	{
		in >> count >> nPhases;
		if (in.fail())
		{
			programlog::writeErr(read_file_err + " - Error of reading elements count");
			return false;
		}
	}
	else
		if (nPhases == -1)
		{
			programlog::writeErr(read_file_err + " - there is no information about phases count");
			return false;
		}

	if (count < 0)
	{
		programlog::writeErr(read_file_err + " - area count cannot be < 0");
		return false;
	}

	if (nPhases < 0)
	{
		programlog::writeErr(read_file_err + " - phases count cannot be < 0");
		return false;
	}

	this->count = count;
	return Storage::read(in, path, count);
}

bool AreaStorage::allocateMemory()
{
	areas = new Area[count];
	return true;
}

bool AreaStorage::readData(std::ifstream& in, int i)
{
	Area areaBuf;

	in >> areaBuf.lambda >> areaBuf.gamma >> areaBuf.lXW >> areaBuf.rXW >> areaBuf.lYW >> areaBuf.rYW >> areaBuf.lZW >> areaBuf.rZW;
	areaBuf.phaseStorage.count = nPhases;
	areaBuf.phaseStorage.phases = new Phase[nPhases];

	Phase phaseBuf;
	for (int j = 0; j < nPhases; j++)
	{
		in >> phaseBuf.dynamicVisc >> phaseBuf.saturation;
		areaBuf.phaseStorage.phases[j] = phaseBuf;
	}

	areas[i] = areaBuf;
	return true;
}

AreaStorage::AreaStorage() { nPhases = -1; }


//---------------SepParamStorage-------------//
bool SepParamStorage::allocateMemory()
{
	sepParams = new SepParam[count];
	return true;
}

bool SepParamStorage::readData(std::ifstream& in, int i)
{
	SepParam sepParamBuf;
	in >> sepParamBuf.n >> sepParamBuf.q;
	if (sepParamBuf.q < 0)
		sepParamBuf.q = 1 / abs(sepParamBuf.q);
	sepParams[i] = sepParamBuf;
	return true;
}

void SepParamStorage::multiplySeParams(int koeff)
{
	for (int i = 0; i < count; i++)
	{
		sepParams[i].n *= koeff;
		sepParams[i].q = pow(sepParams[i].q, 1.0 / koeff);
	}
}

//---------------CoordStorage-------------//

bool CoordStorage::copyElemsByInd(int* masInd, Coord* coordRes, int nResCoord)
{
	for (int i = 0; i < nResCoord; i++)
		coordRes[i] = coords[masInd[i]];

	return true;
}

bool CoordStorage::allocateMemory()
{
	coords = new Coord[count];
	return true;
}

bool CoordStorage::readData(std::ifstream& in, int i)
{
	in >> coords[i].x >> coords[i].y >> coords[i].z;
	return true;
}

void CoordStorage::write(std::ofstream &out)
{
	out << "_________________________________________________________________" << std::endl;
	out << "|---N of Coord--|-------x-------|-------y-------|-------z-------|" << std::endl;
	out << "_________________________________________________________________" << std::endl;
	//out.width(15);
	for (int i = 0; i < count; i++)
		out << "|" << std::setw(15) << i << "|" << std::setw(15) << coords[i].x << "|" << std::setw(15) << coords[i].y << "|" << std::setw(15) << coords[i].z << "|" << std::endl;
	out << "_________________________________________________________________" << std::endl;
}
//---------------CoordStorageSpecif-------------//
bool CoordStorageSpecif::allocateMemory()
{
	coords = new Coord[count];
	return true;
}

bool CoordStorageSpecif::readData(std::ifstream& in, int i)
{
	Coord coordBuf;
	in >> coordBuf.x >> coordBuf.y >> coordBuf.z;
	coords[i] = coordBuf;

	return true;
}

bool CoordStorageSpecif::read(std::string path, int count)
{
	std::string read_file_err = "Error of reading file " + path;
	std::ifstream in;
	in.open(path);

	return read(in, path);
}

//---------------CoordStorageXYZ-------------//
bool CoordStorageXYZ::init(int nX, int nY, int nZ)
{
	this->nX = nX;
	this->nY = nY;
	this->nZ = nZ;

	this->count = nX * nY * nZ;
	this->knotsInLine = nX;
	this->knotsInPlane = nY * knotsInLine;
	this->nKnots = nZ * knotsInPlane;
	this->nKnotsRes = this->nKnots;
	return this->allocateMemory();
}

int CoordStorageXYZ::getKnotIndex(int ix, int jy, int kz)
{
	return kz * knotsInPlane + jy * knotsInLine + ix;
}

bool CoordStorageXYZ::write(std::string fileName)
{
	int i, j, k, p;
	std::ofstream out;
	out.open(fileName);
	out << nX << " " << nY << " " << nZ << std::endl;
	for (k = 0, p = 0; k < nZ; k++)
	{
		for (j = 0; j < nY; j++)
		{
			for (i = 0; i < nX; i++, p++)
				out << coords[p].x << " " << coords[p].y << " " << coords[p].z << "  ";
			out << std::endl;
		}
		out << std::endl;
	}
	return 0;
}

bool CoordStorageXYZ::fillFinitElems(FinitElementStore& finitElementStore)
{
	int nX1 = nX - 1, nY1 = nY - 1, nZ1 = nZ - 1;
	int nFE = finitElemNum();
	finitElementStore.nFinitElement = nFE;
	finitElementStore.finitElements = new FinitElement[nFE];

	int i, j, k;
	FinitElement FinitElementBuf;
	for (k = 0; k < nZ1; k++)
	{
		for (j = 0; j < nY1; j++)
		{
			for (i = 0; i < nX1; i++)
			{
				FinitElementBuf.ver[0] = getKnotIndex(i, j, k);
				FinitElementBuf.ver[1] = getKnotIndex(i + 1, j, k);
				FinitElementBuf.ver[2] = getKnotIndex(i, j + 1, k);
				FinitElementBuf.ver[3] = getKnotIndex(i + 1, j + 1, k);
				FinitElementBuf.ver[4] = getKnotIndex(i, j, k + 1);
				FinitElementBuf.ver[5] = getKnotIndex(i + 1, j, k + 1);
				FinitElementBuf.ver[6] = getKnotIndex(i, j + 1, k + 1);
				FinitElementBuf.ver[7] = getKnotIndex(i + 1, j + 1, k + 1);
			}
		}
	}
	return true;
}

bool CoordStorageXYZ::read(std::ifstream& in, std::string path, int count)
{
	std::string read_file_err = "Error of reading file " + path;

	if (!in.is_open())
	{
		programlog::writeErr(read_file_err + " - Cannot open the file");
		return false;
	}

	if (count == -1)
	{
		in >> nX >> nY >> nZ;
		if (in.fail())
		{
			programlog::writeErr(read_file_err + " - Error of reading elements count");
			return false;
		}
	}

	if (nX < 0)
	{
		programlog::writeErr(read_file_err + " - coordinates count on X cannot be < 0");
		return false;
	}

	if (nY < 0)
	{
		programlog::writeErr(read_file_err + " - coordinates count on Y cannot be < 0");
		return false;
	}

	if (nZ < 0)
	{
		programlog::writeErr(read_file_err + " - coordinates count on Z cannot be < 0");
		return false;
	}

	this->count = count = nX * nY * nZ;
	knotsInLine = nX;
	knotsInPlane = nY * knotsInLine;
	nKnots = nZ * knotsInPlane;

	return Storage::read(in, path, count);

}

void CoordStorageXYZ::getAreaFinElemsNums(Area area, int* &indFinElem, int& nIndFinElems)
{
	int lx = area.lXW;
	int ly = area.lYW;
	int lz = area.lZW;
	int rx = area.rXW;
	int ry = area.rYW;
	int rz = area.rZW;

	int i, j, k;
	int indKnot;
	nIndFinElems = (rx - lx) * (ry - ly) * (rz - lz);
	indFinElem = new int[nIndFinElems];
	int nMissKnotsElemsInSloy = nY + nX - 1;
	int iEl;
	for (i = lx, iEl = 0; i < rx; i++)
	{
		for (j = ly; j < ry; j++)
		{
			for (k = lz; k < rz; k++, iEl++)
			{
				indKnot = getKnotIndex(i, j, k);
				indKnot -= k * (nMissKnotsElemsInSloy) + j;
				indFinElem[iEl] = indKnot;
			}
		}
	}
}

void CoordStorageXYZ::getAreaKnotsNums(Area area, int*& indKnot, int& nIndKnots)
{
	int lx = area.lXW;
	int ly = area.lYW;
	int lz = area.lZW;
	int rx = area.rXW + 1;
	int ry = area.rYW + 1;
	int rz = area.rZW + 1;
	

	int i, j, k;
	nIndKnots = (rx - lx) * (ry - ly) * (rz - lz);
	indKnot = new int[nIndKnots];
	int nMissKnotsElemsInSloy = nY + nX - 1;
	int iKnot;
	for (k = lz, iKnot = 0; k < rz; k++)
		for (j = ly; j < ry; j++)
			for (i = lx; i < rx; i++, iKnot++)
				indKnot[iKnot] = getKnotIndex(i, j, k);	
}

int CoordStorageXYZ::finitElemNum()
{
	return (nX - 1) * (nY - 1) * (nZ - 1);
}

//---------------ContainerSepParam-------------//
ContainerSepParam::ContainerSepParam(int nSepX, int nSepY, int nSepZ, int nSepW, std::string fileName) { readSepParams(nSepX, nSepY, nSepZ, nSepW, fileName); }

ContainerSepParam::ContainerSepParam() {};

bool ContainerSepParam::readSepParams(int nSepX, int nSepY, int nSepZ, int nSepW, std::string fileName)
{
	std::string read_file_err = "Error of reading file " + fileName;
	std::ifstream in;
	in.open(fileName);

	if (!in.is_open())
	{
		programlog::writeErr(read_file_err + " - Cannot open the file");
		return false;
	}

	this->nSepX = nSepX;
	sepX.read(in, fileName, nSepX);

	this->nSepY = nSepY;
	sepY.read(in, fileName, nSepY);

	this->nSepZ = nSepZ;
	sepZ.read(in, fileName, nSepZ);

	this->nSepW = nSepW;
	sepW.read(in, fileName, nSepW);

	int coeff;
	in >> coeff;
	if (coeff != 1)
	{
		sepX.multiplySeParams(coeff);
		sepY.multiplySeParams(coeff);
		sepZ.multiplySeParams(coeff);
		sepW.multiplySeParams(coeff);
	}
	return true;
}


//---------------IndCoordStore-------------//
bool IndCoordStore::calcCoordInd(SepParamStorage sepParamStorage)
{
	SepParam* sepParamC = sepParamStorage.sepParams;

	int nSepC = sepParamStorage.count;
	coordInds = new int[nSepC + 1];
	int nC = 0;
	for (int i = 0; i < nSepC; i++)
	{
		coordInds[i] = nC;
		nC += sepParamC[i].n;
	}
	coordInds[nSepC] = nC;
	indCount = nSepC + 1;

	return true;
}

bool IndCoordStore::write(char seperator)
{
	for (int i = 0; i < indCount; i++)
		std::cout << coordInds[i] << seperator;
	std::cout << std::endl << std::endl;

	return true;
}


//---------------ContainerIXYZW-------------//
ContainerIXYZW::ContainerIXYZW(ContainerSepParam containerSepParam) { create(containerSepParam); }
ContainerIXYZW::ContainerIXYZW() {};

void ContainerIXYZW::create(ContainerSepParam containerSepParam)
{
	IXW.calcCoordInd(containerSepParam.sepX);
	IYW.calcCoordInd(containerSepParam.sepY);
	IZW.calcCoordInd(containerSepParam.sepZ);
}

//---------------CalculationArea-------------//
CalculationArea::CalculationArea(std::string fileNameCord, std::string fileNameSep, std::string fileNameBorders, bool writeInCordFile)
{
	setlocale(LC_ALL, "Russian");
	std::ifstream in(fileNameCord);

	/*Считывание основных ломанных*/
	std::cout << "Считывание опорных точек" << std::endl;
	CoordStorageXYZ XYZW;
	XYZW.read(in, fileNameCord);

	/*Считывание информации о подобластях*/
	std::cout << "Считывание информации о подобластях" << std::endl;
	AreaStorage areas;
	areas.read(in, fileNameCord);
	nPhases = areas.nPhases;
	
	//areasOut.read(in, fileNameCord);
	std::cout << "Считывание информации о скважинах" << std::endl;
	WellStorage wellStorage;
	wellStorage.read(in, fileNameCord);

	PhaseStorage phaseStorage;
	phaseStorage.read(in, fileNameCord, areas.nPhases);

	in >> endT;
	in.close();
	/*Считывание параметров разбиения*/
	//in.open(fileNameSep);

	ContainerSepParam containerSepParam;
	containerSepParam.readSepParams(XYZW.nX - 1, XYZW.nY - 1, XYZW.nZ - 1, wellStorage.count, fileNameSep);

	std::cout << "Построение сетки" << std::endl;
	CrushedMeshCoordStorage XYZ;
	XYZ.initDividesStore(XYZW, containerSepParam);

	XYZ.write("FullCords.txt");
	std::cout << "Встройка скважин" << std::endl;
	wellStorage.insertWells(XYZ, containerSepParam.sepW.sepParams, areas, XYZ.IXYZW);

	fillGeneralCoords(XYZ, wellStorage);
	std::ofstream out;
	out.open("OutCoords.txt");
	coordsStore.write(out);
	out.close();
	std::cout << "Чтение границ с краевыми условиями" << std::endl;
	ContainerBorders borders;
	borders.readBorders(fileNameBorders);


	std::cout << "Иницилизация конечных элементов" << std::endl;
	//XYZ.fillFinitElems(finitElementStore, areas);
	fillGeneralFinElems(XYZ, wellStorage, areas);
	out.open("OutFinitElems.txt");
	finitElementStore.write(out);
	out.close();
	finitElementStore.calcSquares(coordsStore.coords);

	std::cout << "Определение граней" << std::endl;
	faceStore.init(finitElementStore, coordsStore.count);

	std::cout << "Определение граней" << std::endl;
	faceStore.createFaceStore(finitElementStore, borderFacesStore);
	borderFacesStore.phaseStorage = phaseStorage;

	std::cout << "Выделение поамяти для хранения значений потоков" << std::endl;
	flowStore.init(faceStore);

	std::cout << "Заполнение информации об узлах с первыми краевыми условиями" << std::endl;
	fill1CondKnots(XYZ, borders, wellStorage);

	std::cout << "Заполнение информации о гранях со вторыми краевыми условиями" << std::endl;
	fillCondFaces( XYZ, borders, wellStorage);



	std::cout << "Выделение памяти под объёмы фаз" << std::endl;
	containerPhaseVol.init(faceStore, areas.nPhases);

	/*Генерация сетки*/
}

bool CalculationArea::fill1CondKnots(CrushedMeshCoordStorage XYZ, ContainerBorders borders, WellStorage wellStorage)
{
	int i_x1, i_x2, j_y1, j_y2, k_z1, k_z2;
	int i_x, j_y, k_z;
	int ind;

	int* IX_w = XYZ.IXYZW.IXW.coordInds;
	int* IY_w = XYZ.IXYZW.IYW.coordInds;
	int* IZ_w = XYZ.IXYZW.IZW.coordInds;

	Border1C* borderFirstCondAr = borders.firstCondStor.borders;
	int nFirstCond = borders.firstCondStor.count;
	Coord* points = XYZ.coords;
	//int b_ind;

	int i;
	int nx, ny, nz;
	int n = 0;
	for (i = 0; i < nFirstCond; i++)
	{
		i_x1 = IX_w[borderFirstCondAr[i].horizL];
		i_x2 = IX_w[borderFirstCondAr[i].horizR];

		j_y1 = IY_w[borderFirstCondAr[i].vertL];
		j_y2 = IY_w[borderFirstCondAr[i].vertR];

		k_z1 = IZ_w[borderFirstCondAr[i].heightL];
		k_z2 = IZ_w[borderFirstCondAr[i].heightR];

		nx = i_x2 - i_x1 + 1;
		ny = j_y2 - j_y1 + 1;
		nz = k_z2 - k_z1 + 1;
		n += nx * ny * nz;
	}
	BorderStorageWell borderStorageWell = borders.firstCondStorWell;
	int nCondWells = borderStorageWell.count;
	int* iCondWells = borderStorageWell.iWells;

	for (i = 0; i < nCondWells; i++)
		n += wellStorage.wells[iCondWells[i]].getNumWellKnots();


	knots1Cond.IKnots = new int[n];
	knots1Cond.kt1 = n;
	int iknot = 0;
	int* Iknots = knots1Cond.IKnots;
	for (i = 0; i < nFirstCond; i++)
	{
		i_x1 = IX_w[borderFirstCondAr[i].horizL];
		i_x2 = IX_w[borderFirstCondAr[i].horizR];

		j_y1 = IY_w[borderFirstCondAr[i].vertL];
		j_y2 = IY_w[borderFirstCondAr[i].vertR];

		k_z1 = IZ_w[borderFirstCondAr[i].heightL];
		k_z2 = IZ_w[borderFirstCondAr[i].heightR];

		i_x2++;
		j_y2++;
		k_z2++;

		for (k_z = k_z1; k_z < k_z2; k_z++)
		{
			for (j_y = j_y1; j_y < j_y2; j_y++)
			{
				for (i_x = i_x1; i_x < i_x2; i_x++, iknot++)
				{
					ind = XYZ.getKnotIndex(i_x, j_y, k_z);
					Iknots[iknot] = ind;
				}
			}
		}
	}

	Well wellBuf;
	int iCoord, jCircle, kSirf;
	int coordBeg = 0, coordEnd;
	int sirfBeg = 0, sirfEnd;
	for (i = 0; i < nCondWells; i++)
	{
		wellBuf = wellStorage.wells[iCondWells[i]];
		coordEnd = wellBuf.nCircleCoords;
		sirfEnd = wellBuf.nSirfs;
		jCircle = wellBuf.nCircle - 1;
		int iCoord, kSirf;
		for (kSirf = sirfBeg; kSirf < sirfEnd; kSirf++)
		{
			for (iCoord = coordBeg; iCoord < coordEnd; iCoord++, iknot++)
			{
				Iknots[iknot] = wellBuf.globI[wellBuf.getKnotIndex(iCoord, jCircle, kSirf)];
			}
		}

	}
	return 0;
}


CalculationArea::CalculationArea() {};

bool CalculationArea::fillCondFaces(CrushedMeshCoordStorage XYZ, ContainerBorders borders, WellStorage wellStorage)
{
	FaceStore Faces2CondStore; // Грани со вторыми краевыми
	int i_x1, i_x2, j_y1, j_y2, k_z1, k_z2;
	int i_x, j_y, k_z;
	int ind;

	int horL, horR, verL, verR, heiL, heiR;
	int iHor, jVer, kHei;

	int* IX_w = XYZ.IXYZW.IXW.coordInds;
	int* IY_w = XYZ.IXYZW.IYW.coordInds;
	int* IZ_w = XYZ.IXYZW.IZW.coordInds;

	Border* borderSecondCondAr = borders.secondCondStor.borders;

	int nSecondCond = borders.secondCondStor.count;
	Coord* points = XYZ.coords;
	int b_ind;

	int i;
	int nx, ny, nz;
	int n = 0;
	for (i = 0; i < nSecondCond; i++)
	{
		i_x1 = IX_w[borderSecondCondAr[i].horizL];
		i_x2 = IX_w[borderSecondCondAr[i].horizR];

		j_y1 = IY_w[borderSecondCondAr[i].vertL];
		j_y2 = IY_w[borderSecondCondAr[i].vertR];

		k_z1 = IZ_w[borderSecondCondAr[i].heightL];
		k_z2 = IZ_w[borderSecondCondAr[i].heightR];

		nx = i_x2 - i_x1;
		ny = j_y2 - j_y1;
		nz = k_z2 - k_z1;

		nx = (nx == 0) ? 1 : nx;
		ny = (ny == 0) ? 1 : ny;
		nz = (nz == 0) ? 1 : nz;

		n += nx * ny * nz;
	}
	BorderStorageWell borderStorageWell = borders.secondCondStorWell;
	int nCondWells = borderStorageWell.count;
	int* iCondWells = borderStorageWell.iWells;
	for (i = 0; i < nCondWells; i++)
		n += wellStorage.wells[iCondWells[i]].getNumWellFaces();

	Faces2CondStore.faces = new Face[n];
	Faces2CondStore.count = n;
	Face* cond2Faces = Faces2CondStore.faces;
	int iFace = 0;
	int* localIKnots;

	for (i = 0; i < nSecondCond; i++)
	{
		horL = IX_w[borderSecondCondAr[i].horizL];
		horR = IX_w[borderSecondCondAr[i].horizR];

		verL = IY_w[borderSecondCondAr[i].vertL];
		verR = IY_w[borderSecondCondAr[i].vertR];

		heiL = IZ_w[borderSecondCondAr[i].heightL];
		heiR = IZ_w[borderSecondCondAr[i].heightR];

		if (horL == horR)
		{
			//horR++;
			//verR++;
			//heiR++;
			iHor = horL;

			for (kHei = heiL; kHei < heiR; kHei++)
			{
				for (jVer = verL; jVer < verR; jVer++, iFace++)
				{
					localIKnots = cond2Faces[iFace].knots;
					localIKnots[0] = XYZ.getKnotIndex(iHor, jVer, kHei);
					localIKnots[1] = XYZ.getKnotIndex(iHor, jVer + 1, kHei);
					localIKnots[2] = XYZ.getKnotIndex(iHor, jVer, kHei + 1);
					localIKnots[3] = XYZ.getKnotIndex(iHor, jVer + 1, kHei + 1);

				}
			}
		}

		if (verL == verR)
		{
			//horR++;
			//verR++;
			//heiR++;
			jVer = verL;

			for (kHei = heiL; kHei < heiR; kHei++)
			{
				for (iHor = horL; iHor < horR; iHor++, iFace++)
				{
					localIKnots = cond2Faces[iFace].knots;
					localIKnots[0] = XYZ.getKnotIndex(iHor, jVer, kHei);
					localIKnots[1] = XYZ.getKnotIndex(iHor + 1, jVer, kHei);
					localIKnots[2] = XYZ.getKnotIndex(iHor, jVer, kHei + 1);
					localIKnots[3] = XYZ.getKnotIndex(iHor + 1, jVer + 1, kHei + 1);
				}
			}
		}

		if (heiL == heiR)
		{
			//horR++;
			//verR++;
			//heiR++;
			kHei = heiL;

			for (jVer = verL; jVer < verR; jVer++)
			{
				for (iHor = horL; iHor < horR; iHor++, iFace++)
				{
					localIKnots = cond2Faces[iFace].knots;
					localIKnots[0] = XYZ.getKnotIndex(iHor, jVer, kHei);
					localIKnots[1] = XYZ.getKnotIndex(iHor + 1, jVer, kHei);
					localIKnots[2] = XYZ.getKnotIndex(iHor, jVer + 1, kHei);
					localIKnots[3] = XYZ.getKnotIndex(iHor + 1, jVer + 1, kHei);
				}
			}
		}
		
	}

	Well wellBuf;
	int iCoord, jCircle, kSirf;
	int coordBeg = 0, coordEnd;
	int sirfBeg = 0, sirfEnd;
	int* globI;
	for (i = 0; i < nCondWells; i++)
	{
		wellBuf = wellStorage.wells[iCondWells[i]];
		coordEnd = wellBuf.nCircleCoords - 1;
		sirfEnd = wellBuf.nSirfs - 1;
		jCircle = wellBuf.nCircle - 1;
		int iCoord, kSirf;
		globI = wellBuf.globI;
		for (kSirf = sirfBeg; kSirf < sirfEnd; kSirf++)
		{
			for (iCoord = coordBeg; iCoord < coordEnd; iCoord++, iFace++)
			{
				localIKnots = cond2Faces[iFace].knots;
				localIKnots[0] = globI[wellBuf.getKnotIndex(iCoord, jCircle, kSirf)];
				localIKnots[1] = globI[wellBuf.getKnotIndex(iCoord + 1, jCircle, kSirf)];
				localIKnots[2] = globI[wellBuf.getKnotIndex(iCoord, jCircle, kSirf + 1)];
				localIKnots[3] = globI[wellBuf.getKnotIndex(iCoord + 1, jCircle, kSirf + 1)];
			}

			localIKnots = cond2Faces[iFace].knots;
			localIKnots[0] = globI[wellBuf.getKnotIndex(iCoord, jCircle, kSirf)];
			localIKnots[1] = globI[wellBuf.getKnotIndex(coordBeg, jCircle, kSirf)];
			localIKnots[2] = globI[wellBuf.getKnotIndex(iCoord, jCircle, kSirf + 1)];
			localIKnots[3] = globI[wellBuf.getKnotIndex(coordBeg, jCircle, kSirf + 1)];
			iFace++;
		}

	}

	faces2CondStore.nFaces = n;
	faces2CondStore.iFaces = new int[n];
	int* iFaces = faces2CondStore.iFaces;
	for (i = 0; i < n; i++)
		iFaces[i] = faceStore.findFaceIndex(cond2Faces[i]);

	return 0;
}

bool CalculationArea::fillGeneralFinElems(CrushedMeshCoordStorage XYZ, WellStorage wellStorage, AreaStorage areas)
{
	FinitElementStore finElStoreXYZ;
	FinitElementStore finElStoreWell;
	FinitElementStore finElStoreBuf;

	XYZ.fillFinitElems(finElStoreXYZ, areas);
	
	int nWells = wellStorage.count;
	Area area;
	int nBufElems;
	int i, j, k;

	for (i = 0, nBufElems = 0; i < nWells; i++)
		nBufElems += wellStorage.wells[i].finitElemNum();

	nBufElems += finElStoreXYZ.nFinitElement;

	finElStoreBuf.nFinitElement = nBufElems;
	finElStoreBuf.finitElements = new FinitElement[nBufElems];
	FinitElement* finElemsXYZ = finElStoreXYZ.finitElements;
	FinitElement* finElemsBuf = finElStoreBuf.finitElements;

	for (j = 0; j < finElStoreXYZ.nFinitElement; j++)
		finElemsBuf[j] = finElemsXYZ[j];

	int nElems = finElStoreXYZ.nFinitElement;
	int* indDelAreas;
	int nIndDelAreas;
	FinitElement* finElemsWell;
	int nFinElemsWell;
	for (i = 0; i < nWells; i++)
	{
		wellStorage.wells[i].fillFinitElems(finElStoreWell);
		finElemsWell = finElStoreWell.finitElements;
		nFinElemsWell = finElStoreWell.nFinitElement;
		area = wellStorage.wells[i].getXYZMeshInd();
		XYZ.getAreaFinElemsNums(area, indDelAreas, nIndDelAreas);

		for (j = 0; j < nIndDelAreas; j++)
			finElemsBuf[indDelAreas[j]] = finElemsWell[j];

		for (; j < nFinElemsWell; j++, nElems++)
		{
			finElemsBuf[nElems] = finElemsWell[j];
		}
		finElStoreWell.clear();
		delete[] indDelAreas;
	}
	finitElementStore.finitElements = new FinitElement[nElems];
	finitElementStore.nFinitElement = nElems;
	FinitElement* finitElemGen = finitElementStore.finitElements;

	for (i = 0; i < nElems; i++)
		finitElemGen[i] = finElemsBuf[i];
	
	return true;
}

bool CalculationArea::fillGeneralCoords(CoordStorage XYZ, WellStorage wellStorage)
{
	Coord* coordsXYZ;
	Coord* coordWell;
	Coord* coord;
	int* iCoordWell;
	int i, j;

	//CoordStorage coordStorageBuf;
	int nCoordStorage = wellStorage.iCoordMax;
	int nWell = wellStorage.count;
	//for (i = 0; i < nWell; i++)
	//	nCoordStorageBuf += wellStorage.wells[i].count;
	
	coordsStore.count = nCoordStorage;
	coordsStore.coords = new Coord[nCoordStorage];
	coord = coordsStore.coords;
	coordsXYZ = XYZ.coords;
	for (i = 0; i < XYZ.count; i++)
		coord[i] = coordsXYZ[i];
	
	int nWellCoord;
	for (i = 0; i < nWell; i++)
	{
		coordWell = wellStorage.wells[i].coords;
		iCoordWell = wellStorage.wells[i].globI;
		nWellCoord = wellStorage.wells[i].count;

		for (j = 0; j < nWellCoord; j++)
			coord[iCoordWell[j]] = coordWell[j];
		
	}
	return true;
}
//---------------CrushedMeshCoordStorage-------------//
CrushedMeshCoordStorage::CrushedMeshCoordStorage() {};
CrushedMeshCoordStorage::CrushedMeshCoordStorage(CoordStorageXYZ XYZW, ContainerSepParam containerSepParam) { initDividesStore(XYZW, containerSepParam); }

void CrushedMeshCoordStorage::initDividesStore(CoordStorageXYZ XYZW, ContainerSepParam containerSepParam)
{
	IXYZW.create(containerSepParam);

	nX = IXYZW.IXW.coordInds[IXYZW.IXW.indCount - 1] + 1;
	nY = IXYZW.IYW.coordInds[IXYZW.IYW.indCount - 1] + 1;
	nZ = IXYZW.IZW.coordInds[IXYZW.IZW.indCount - 1] + 1;

	this->init(nX, nY, nZ);

	knotsInLine = nX; // Количество точек в горизонтальной линии
	knotsInPlane = nY * knotsInLine; // Количество точек в вертикальной линии
	nKnots = nZ * knotsInPlane;

	for (int i = 0; i < nKnots; i++)
	{
		coords[i].x = coords[i].y = coords[i].z = -1;
	}
		
	int knotsWInLine = XYZW.nX; // Количество основных точек в горизонтальной линии
	int knotsWInPlane = XYZW.nY * knotsWInLine; // Количество основных точек в вертикальной линии

	SepParam* sepX = containerSepParam.sepX.sepParams;
	SepParam* sepY = containerSepParam.sepY.sepParams;
	SepParam* sepZ = containerSepParam.sepZ.sepParams;

	SepParam sepParamX, sepParamY, sepParamZ;


	int pCurr, // текущая точка (общие координаты)
		pNextHeight, // следующая точка по высоте (общие координаты)
		pNextHoriz, // следующая горизонтальная точка (общие координаты)
		pNextVert; // следующая вертикальная точка (общие координаты)


	Coord dCoord;
	Coord coord;

	int iX = 0, iY = 0, iZ = 0;

	int nXW = XYZW.nX,
		nYW = XYZW.nY,
		nZW = XYZW.nZ;

	int* IXW = IXYZW.IXW.coordInds;
	int* IYW = IXYZW.IYW.coordInds;
	int* IZW = IXYZW.IZW.coordInds;

	int iXW, iYW, iZW;
	double qHoriz, qVert, qHeight;
	double nHoriz, nVert, nHeight;

	int pWCurr, // текущая точка (основные координаты)
		pWNextHeight, // следующая точка по высоте (основные координаты)
		pWNextHoriz, // следующая горизонтальная точка (основные координаты)
		pWNextVert; // следующая вертикальная точка (основные координаты)

	double dx, dy, dz;
	double x, y, z;
	int i, j, k, c;
	double koef;
	int n_x, n_y, n_z;

	/*
	std::ifstream in1;
	in1.open("FullCord.txt");
	in1 >> n_x;
	in1 >> n_y;
	in1 >> n_z;
	for (k = 0, c = 0; k < n_z; k++)
	{
		for (j = 0; j < n_y; j++)
		{
			for (i = 0; i < n_x; i++, c++)
			{

				in1 >> x >> y >> z;
				coord.x = x;
				coord.y = y;
				coord.z = z;

				this->data[c] = coord;
			}
		}
	}
	*/


	copyBasePoints(XYZW);

	for (iZW = 0; iZW < nZW; iZW++)
	{
		iZ = IZW[iZW];
		for (iYW = 0; iYW < nYW; iYW++)
		{
			iY = IYW[iYW];
			for (iXW = 0; iXW < nXW - 1; iXW++)
			{
				sepParamX = sepX[iXW];
				iX = IXW[iXW];

				pCurr = getKnotIndex(iX, iY, iZ);
				pNextHoriz = pCurr + sepX[iXW].n;

				copyLineOfBasePoints(sepParamX, XYZW, pNextHoriz, pCurr);
			}
		}
	}

	for (iZW = 0; iZW < nZW; iZW++)
	{
		iZ = IZW[iZW];
		for (iYW = 0; iYW < nYW - 1; iYW++)
		{
			iY = IYW[iYW];
			sepParamY = sepY[iYW];
			for (iXW = 0; iXW < nXW - 1; iXW++)
			{
				for (iX = IXW[iXW]; iX < IXW[iXW + 1] + 1; iX++)
				{
					pCurr = getKnotIndex(iX, iY, iZ);
					//pNextVert = getKnotIndex(iX, iY + 1, iZ);
					pNextVert = pCurr + sepParamY.n * knotsInLine;
					copyLineOfBasePoints(sepParamY, XYZW, pNextVert, pCurr, knotsInLine);
				}
			}
		}
	}

	for (iZW = 0; iZW < nZW - 1; iZW++)
	{
		iZ = IZW[iZW];
		sepParamZ = sepZ[iZW];
		for (iYW = 0; iYW < nYW - 1; iYW++)
		{
			for (iY = IYW[iYW]; iY < IYW[iYW + 1] + 1; iY++)
			{
				for (iXW = 0; iXW < nXW - 1; iXW++)
				{
					for (iX = IXW[iXW]; iX < IXW[iXW + 1] + 1; iX++)
					{

						pCurr = getKnotIndex(iX, iY, iZ);
						pNextHeight = pCurr + sepParamZ.n * knotsInPlane;
						copyLineOfBasePoints(sepParamZ, XYZW, pNextHeight, pCurr, knotsInPlane);

					}
				}
			}
		}
	}
}

bool CrushedMeshCoordStorage::fillFinitElems(FinitElementStore& finitElementStore, AreaStorage& AreaStorage)
{
	//int n_X = mesh.XYZ.nX;
	//int n_Y = mesh.XYZ.nY;
	//int n_Z = mesh.XYZ.nZ;
	//double* X = mesh.X;
	//double* Y = mesh.Y;
	//double* Z = mesh.Z;
	//Coord* XYZarray = mesh.XYZ.data;

	Area* areas = AreaStorage.areas;
	Area areaCur;
	int N_ar = AreaStorage.count;

	int* IX_w = IXYZW.IXW.coordInds;
	int* IY_w = IXYZW.IYW.coordInds;
	int* IZ_w = IXYZW.IZW.coordInds;

	int ind;
	
	finitElementStore.nFinitElement = (nX - 1) * (nY - 1) * (nZ - 1);
	finitElementStore.finitElements = new FinitElement[finitElementStore.nFinitElement];
	FinitElement* finitElements = finitElementStore.finitElements;
	//int n_X1 = n_X - 1;

	int r, l, m, d/*, b_i, b_j*/;
	int i, j, k;
	//int block_ind;
	int i_beg, i_end, j_beg, j_end, k_beg, k_end;
	double lambda_r, /*hi_r, delta_r,*/ gamma_r;

	FinitElement FinitElementBuf;
	//FinalElemBase finalElemBase;
	int nPhases = AreaStorage.nPhases;
	int iAr = 0;
	for (r = 0; r < N_ar; r++)
	{
		areaCur = areas[r];
		i_beg = IX_w[areaCur.lXW];
		i_end = IX_w[areaCur.rXW];
		j_beg = IY_w[areaCur.lYW];
		j_end = IY_w[areaCur.rYW];
		k_beg = IZ_w[areaCur.lZW];
		k_end = IZ_w[areaCur.rZW];
		lambda_r = areaCur.lambda;
		gamma_r = areaCur.gamma;
		for (k = k_beg; k < k_end; k++) // по высоте
		{
			for (j = j_beg; j < j_end; j++) // по вертикали
			{
				for (i = i_beg; i < i_end; i++, iAr++) // по горизонтали
				{
					FinitElementBuf.ver[0] = getKnotIndex(i, j, k);
					FinitElementBuf.ver[1] = getKnotIndex(i + 1, j, k);
					FinitElementBuf.ver[2] = getKnotIndex(i, j + 1, k);
					FinitElementBuf.ver[3] = getKnotIndex(i + 1, j + 1, k);
					FinitElementBuf.ver[4] = getKnotIndex(i, j, k + 1);
					FinitElementBuf.ver[5] = getKnotIndex(i + 1, j, k + 1);
					FinitElementBuf.ver[6] = getKnotIndex(i, j + 1, k + 1);
					FinitElementBuf.ver[7] = getKnotIndex(i + 1, j + 1, k + 1);

					FinitElementBuf.phaseStorage.initStorage(areaCur.phaseStorage.phases, nPhases);
					FinitElementBuf.K = lambda_r;
					FinitElementBuf.FI = gamma_r;

					finitElements[iAr] = FinitElementBuf;
					//finalElemBase.setFinalElement(FinitElementBuf, *this);
					//FinitElementBuf.sqare = finalElemBase.calcSqare();
				}
			}
		}
	}

	return true;
}

void CrushedMeshCoordStorage::copyBasePoints(CoordStorageXYZ& XYZW)
{
	Coord* XYZ = this->coords;
	Coord* XYZWc = XYZW.coords;
	int iXW, iYW, iZW;
	int* IXW = IXYZW.IXW.coordInds;
	int* IYW = IXYZW.IYW.coordInds;
	int* IZW = IXYZW.IZW.coordInds;
	int nXW = XYZW.nX,
		nYW = XYZW.nY,
		nZW = XYZW.nZ;

	int pCurr, pWCurr;
	for (iZW = 0; iZW < nZW; iZW++)
	{
		for (iYW = 0; iYW < nYW; iYW++)
		{
			for (iXW = 0; iXW < nXW; iXW++)
			{
				pWCurr = XYZW.getKnotIndex(iXW, iYW, iZW);
				pCurr = getKnotIndex(IXW[iXW], IYW[iYW], IZW[iZW]);
				XYZ[pCurr] = XYZWc[pWCurr];
			}
		}
	}
}

void CrushedMeshCoordStorage::copyLineOfBasePoints(SepParam sepParam, CoordStorageSpecif& XYZW, int pNext, int pCurr, int pStep)
{
	Coord* XYZ = this->coords;
	Coord coord, dCoord;
	double koef;
	double q = sepParam.q;
	int n = sepParam.n;

	if (q == 1)
	{
		dCoord = (XYZ[pNext] - XYZ[pCurr]) / n;

	}
	else
	{
		koef = (1 - q) / (1 - pow(q, n));
		dCoord = (XYZ[pNext] - XYZ[pCurr]) * koef;
	}

	coord = XYZ[pCurr] + dCoord;
	dCoord *= q;

	for (int p = pCurr + pStep; p < pNext; p += pStep)
	{
		XYZ[p] = coord;
		coord += dCoord;
		dCoord *= q;
	}
}

//---------------AreaPhaseStorage-------------//
bool AreaPhaseStorage::allocateMemory()
{
	areaPhases = new AreaPhase[count];
	return true;
}

bool AreaPhaseStorage::readData(std::ifstream& in, int i)
{
	AreaPhase AreaPhaseBuf;
	in >> AreaPhaseBuf.iArea >> AreaPhaseBuf.iPhases >> AreaPhaseBuf.PhaseSaturation;
	return true;
}

FaceStore::FaceStore()
{
	nBufLookingCenter = 0;
	nBufLookingArea = 10;
	ig = NULL;
	jg = NULL;
	n = 0;
	count = 0;
	faces = NULL;
}

//---------------FacesStore-------------//
bool FaceStore::init(FinitElementStore& finitElementStore, int nCoord)
{
	FinitElement* finitElems = finitElementStore.finitElements;
	int nFinitElems = finitElementStore.nFinitElement;
	n = nCoord;
	ig = new int[n + 1];
	int* jgHelp;
	int nBuf = n;
	jgHelp = new int[nBuf];
	arrayspace::fill_vec(ig, n + 1, 0);
	int* verBuf;
	int i, j;

	for (i = 0; i < nFinitElems; i++)
	{
		verBuf = finitElems[i].ver;
		for (j = 0; j < VER_NUM; j++)
		{
			ig[verBuf[j]]++;
		}
	}
	int njg = 0;
	int curInd;
	for (i = 0; i < n; i++)
	{
		curInd = njg;
		njg += ig[i];
		ig[i] = curInd;
	}
	ig[i] = njg;
	jg = new int[njg];
	arrayspace::fill_vec(jgHelp, n, 0);

	for (i = 0; i < nFinitElems; i++)
	{
		verBuf = finitElems[i].ver;
		for (j = 0; j < VER_NUM; j++)
		{
			int iVer = verBuf[j];
			int ijg = ig[iVer] + jgHelp[iVer];
			jg[ijg] = i;
			jgHelp[iVer]++;
		}
	}

	int iBeg;
	for (i = 0; i < n; i++)
	{
		iBeg = ig[i];
		arrayspace::sort(jg, iBeg, ig[i + 1] - iBeg);
	}

	delete[] jgHelp;

	return true;
}

bool FaceStore::createFaceStore(FinitElementStore& finitElementStore)
{
	FinitElement* finitElems = finitElementStore.finitElements;
	int nFinitElems = finitElementStore.nFinitElement;

	int i, j, k;
	int iNeighbFinitElem;
	FaceStore facesStoreBuf;
	Face FaceBuf;
	Face FaceBufNeigb;
	facesStoreBuf.faces = new Face[FACES_NUM * finitElementStore.nFinitElement];
	facesStoreBuf.n = FACES_NUM * finitElementStore.nFinitElement;
	int nFaces = 0;
	int* NeighbVerFace = new int[VER_NUM_FACE];
	int* currVerFace;

	for (i = 0; i < nFinitElems; i++)
	{
		for (j = 0; j < FACES_NUM; j++)
		{
			iNeighbFinitElem = findNeighboringFinitElem(i, finitElems[i], j);
			if (iNeighbFinitElem > i || iNeighbFinitElem < 0) // Если смежный по грани элемент имеет больший глобальный номер или по этой грани у элемента нет соседа 
			{
				finitElems[i].getFaceGlobalNum(j, FaceBuf.knots);
				facesStoreBuf.faces[nFaces] = FaceBuf;
				finitElems[i].faces[j] = nFaces;
				if (iNeighbFinitElem > i)
				{
					currVerFace = FaceBuf.knots;
					for (k = 0; k < FACES_NUM; k++)
					{
						finitElems[iNeighbFinitElem].getFaceGlobalNum(k, NeighbVerFace);
						if (arrayspace::isSameWithoutOrdUniqe(NeighbVerFace, currVerFace, VER_NUM_FACE))
							finitElems[i].faces[k] = nFaces;
					}
				}
				nFaces++;			
			}
			/*
			else
			{
				finitElems[i].getFaceGlobalNum(j, currVerFace);
				for (k = 0; k < VER_NUM_FACE; k++)
				{
					finitElems[iNeighbFinitElem].getFaceGlobalNum(k, NeighbVerFace);
					if (arrayspace::isSameWithoutOrd(NeighbVerFace, currVerFace, VER_NUM_FACE))
						finitElems[i].faces[j] = finitElems[iNeighbFinitElem].faces[k];
				}
			}
			*/
		}
	}
	copyStore(facesStoreBuf, nFaces);
	delete[] facesStoreBuf.faces;
	facesStoreBuf.count = 0;
	return true;
}

bool FaceStore::createFaceStore(FinitElementStore& finitElementStore, BorderFacesStore& borderFacesStore)
{
	FinitElement* finitElems = finitElementStore.finitElements;
	int nFinitElems = finitElementStore.nFinitElement;

	int i, j, k;
	int iNeighbFinitElem;
	FaceStore facesStoreBuf;
	Face FaceBuf;
	Face FaceBufNeigb;

	facesStoreBuf.n = FACES_NUM * finitElementStore.nFinitElement;
	facesStoreBuf.faces = new Face[facesStoreBuf.n];
	Face* bufFaces = facesStoreBuf.faces;

	BorderFacesStore borderFacesStoreBuf;
	borderFacesStoreBuf.iFaces = new int[facesStoreBuf.n];
	borderFacesStoreBuf.nFaces = facesStoreBuf.n;
	int* iBordFaces = borderFacesStoreBuf.iFaces;
	int nBordFace = 0;
	int nFaces = 0;

	int NeighbVerFace[VER_NUM_FACE];
	int* currVerFace;

	for (i = 0; i < nFinitElems; i++)
	{
		for (j = 0; j < FACES_NUM; j++)
		{
			finitElems[i].getFaceGlobalNum(j, FaceBuf.knots);
			iNeighbFinitElem = findNeighboringFinitElem(i, finitElems[i], j);
			if (iNeighbFinitElem > i || iNeighbFinitElem == -1)
			{
				finitElems[i].faces[j] = nFaces;
				if (iNeighbFinitElem == -1)
				{
					iBordFaces[nBordFace] = nFaces;
					nBordFace++;
				}
				arrayspace::copy(facesStoreBuf.faces[nFaces].knots, FaceBuf.knots, VER_NUM_FACE);
				nFaces++;
			}
			else
			{		
				//finitElems[i].faces[j] = nFaces;
				currVerFace = FaceBuf.knots;
				for (k = 0; k < FACES_NUM; k++)
				{
					finitElems[iNeighbFinitElem].getFaceGlobalNum(k, NeighbVerFace);
					if (arrayspace::isSameWithoutOrdUniqe(NeighbVerFace, currVerFace, VER_NUM_FACE))
						finitElems[i].faces[j] = finitElems[iNeighbFinitElem].faces[k];
				}

			}
			/*
			if (iNeighbFinitElem < i || iNeighbFinitElem > -1) // Если смежный по грани элемент имеет больший глобальный номер или по этой грани у элемента нет соседа 
			{
				finitElems[i].getFaceGlobalNum(j, FaceBuf.knots);
				facesStoreBuf.faces[nFaces] = FaceBuf;
				finitElems[i].faces[j] = nFaces;
				if (iNeighbFinitElem != -1)
				{
					currVerFace = FaceBuf.knots;
					for (k = 0; k < FACES_NUM; k++)
					{
						finitElems[iNeighbFinitElem].getFaceGlobalNum(k, NeighbVerFace);
						if (arrayspace::isSameWithoutOrdUniqe(NeighbVerFace, currVerFace, VER_NUM_FACE))
							finitElems[i].faces[k] = nFaces;
					}
				}
				else
				{
					iBordFaces[nBordFace] = nFaces;
					nBordFace++;
				}
				nFaces++;
			}
			*/
			/*
			else
			{
				finitElems[i].getFaceGlobalNum(j, currVerFace);
				for (k = 0; k < VER_NUM_FACE; k++)
				{
					finitElems[iNeighbFinitElem].getFaceGlobalNum(k, NeighbVerFace);
					if (arrayspace::isSameWithoutOrd(NeighbVerFace, currVerFace, VER_NUM_FACE))
						finitElems[i].faces[j] = finitElems[iNeighbFinitElem].faces[k];
				}
			}
			*/
		}
	}
	borderFacesStore.copyBorderFacesStore(borderFacesStoreBuf, nBordFace);
	copyStore(facesStoreBuf, nFaces);
	delete[] facesStoreBuf.faces;
	facesStoreBuf.count = 0;
	return true;
}

void FaceStore::copyStore(FaceStore& facesStore, int nFaces)
{
	if (nFaces > facesStore.n)
		programlog::writeErr("Error of trying copy more elements from the store than there are");

	n = count = nFaces;
	faces = new Face[nFaces];
	Face* otherFaces = facesStore.faces;
	for (int i = 0; i < nFaces; i++)
		faces[i] = otherFaces[i];
}

int FaceStore::findFaceIndex(const Face& face)
{
	int i;
	int nFaces = count;
	int lBufLookindBorder = nBufLookingCenter > nBufLookingArea ? nBufLookingCenter - nBufLookingArea : 0;
	int rBufLookindBorder = nBufLookingCenter + nBufLookingArea < nFaces ? nBufLookingCenter + nBufLookingArea : nFaces;
	//FinitElement* finitElements = calculationArea.finitElementStore.finitElements;
	bool isElemFound;
	double ans;
	for (i = nBufLookingCenter; i < rBufLookindBorder; i++)
	{
		//finitElemCulcer.init(finitElements[i], calculationArea.coordsStore, q);
		//ans = finitElemCulcer.countFunct(p, isElemFound);
		if (arrayspace::isSameWithoutOrdUniqe(faces[i].knots, face.knots, VER_NUM_FACE))
		{
			nBufLookingCenter = i;
			return i;
		}
	}
	for (i = lBufLookindBorder; i < nBufLookingCenter; i++)
	{
		if (arrayspace::isSameWithoutOrdUniqe(faces[i].knots, face.knots, VER_NUM_FACE))
		{
			nBufLookingCenter = i;
			return i;
		}
	}

	for (i = 0; i < nFaces; i++)
	{
		if (arrayspace::isSameWithoutOrdUniqe(faces[i].knots, face.knots, VER_NUM_FACE))
		{
			nBufLookingCenter = i;
			return i;
		}
	}
	programlog::writeErr("The point in which solution is searching is out of Calculation Area");
	return -1;
}

int FaceStore::findNeighboringFinitElem(int indFinElem, FinitElement& finitElement, int iLocalFace)
{
	int faceVer[VER_NUM_FACE];
	finitElement.getFaceGlobalNum(iLocalFace, faceVer);

	int ver0 = faceVer[0];
	int finElemsWithVer1 = ig[ver0 + 1] - ig[ver0];
	int ver0Beg = ig[ver0], ver0End = ig[ver0 + 1];
	int i, j;
	int ifinElem;
	int enoghVerNumFace = 3;
	bool fl = false;
	for (i = ver0Beg; i < ver0End && !fl; i++)
	{
		ifinElem = jg[i];
		if (ifinElem == indFinElem) continue;
		fl = true;
		for (j = 1; j < enoghVerNumFace && fl; j++)
		{
			int verj = faceVer[j];
			fl = arrayspace::serchInSorted(jg, ig[verj], ig[verj + 1] - ig[verj], ifinElem) > -1;
		}
	}

	if (!fl) return -1;

	return ifinElem;
}


//---------------FinitElement-------------//

bool FinitElement::getFaceGlobalNum(int iLocalFace, int faceVer[VER_NUM_FACE])
{
	getFaceLocalNum(iLocalFace, faceVer);

	int i, iVerLocal;

	for (i = 0; i < VER_NUM_FACE; i++)
	{
		iVerLocal = faceVer[i];
		faceVer[i] = ver[iVerLocal];
	}

	return true;
}

bool FinitElement::getFaceLocalNum(int iLocalFace, int faceVer[VER_NUM_FACE])
{
	switch (iLocalFace)
	{
	case 0: 
		faceVer[0] = 0;
		faceVer[1] = 1;
		faceVer[2] = 2;
		faceVer[3] = 3;
		return true;
	case 1:
		faceVer[0] = 0;
		faceVer[1] = 1;
		faceVer[2] = 4;
		faceVer[3] = 5;
		return true;
	case 2:
		faceVer[0] = 2;
		faceVer[1] = 3;
		faceVer[2] = 6;
		faceVer[3] = 7;
		return true;
	case 3:
		faceVer[0] = 0;
		faceVer[1] = 2;
		faceVer[2] = 4;
		faceVer[3] = 6;
		return true;
	case 4:
		faceVer[0] = 1;
		faceVer[1] = 3;
		faceVer[2] = 5;
		faceVer[3] = 7;
		return true;
	case 5:
		faceVer[0] = 4;
		faceVer[1] = 5;
		faceVer[2] = 6;
		faceVer[3] = 7;
		return true;
	default:
		programlog::writeErr("Error: wrong local index of face");
		return false;
	}

	return true;
}

bool FinitElement::getOpposFaceLocalNum(int iLocalFace, int faceVer[VER_NUM_FACE])
{
	switch (iLocalFace)
	{
	case 0:
		getFaceLocalNum(5, faceVer);
		return true;
	case 1:
		getFaceLocalNum(2, faceVer);
		return true;
	case 2:
		getFaceLocalNum(1, faceVer);
		return true;
	case 3:
		getFaceLocalNum(4, faceVer);
		return true;
	case 4:
		getFaceLocalNum(3, faceVer);
		return true;
	case 5:
		getFaceLocalNum(0, faceVer);
		return true;
	default:
		programlog::writeErr("Error: wrong local index of face");
		return false;
	}

	return true;
}


//---------------FlowStore-------------//
void FlowStore::init(FaceStore faceStore)
{
	count = faceStore.count;
	flows = new double[count];
}

//---------------PhaseVolStore-------------//
void PhaseVolStore::init(FaceStore faceStore)
{
	count = faceStore.count;
	PhaseVol = new double[count];
}

//---------------ContainerPhaseVol-------------//
void ContainerPhaseVol::init(FaceStore faceStore, int nPhases)
{
	count = nPhases;
	phaseVolStore = new PhaseVolStore[count];

	for (int i = 0; i < nPhases; i++)
		phaseVolStore[i].init(faceStore);

	phaseVolStoreMix.init(faceStore);
}

//---------------Well-------------//
bool Well::findBegApprCenterPlace(CoordStorageXYZ XYZ, int cordEmptyInd[4])
{
	int ix1, ix2, jy1, jy2;
	double xL = center.x - R, 
		xR = center.x + R,
		yL = center.y - R,
		yR = center.y + R;
	Coord* coords = XYZ.coords;

	int i, j;
	for (ix1 = 0; ix1 < XYZ.nX && coords[XYZ.getKnotIndex(ix1, 0, 0)].x < xL; ix1++);

	//ix1--;
	if (ix1 == XYZ.nX)
		return false;

	ix1--;

	for (ix2 = ix1 + 1; ix2 < XYZ.nX && coords[XYZ.getKnotIndex(ix2, 0, 0)].x < xR; ix2++);

	if (ix2 == XYZ.nX)
		return false;
		


	for (jy1 = 0; jy1 < XYZ.nY && coords[XYZ.getKnotIndex(0, jy1, 0)].y < yL; jy1++);

	//ix1--;
	if (jy1 == XYZ.nY)
		return false;

	jy1--;

	for (jy2 = jy1 + 1; jy2 < XYZ.nY && coords[XYZ.getKnotIndex(0, jy2, 0)].y < yR; jy2++);

	if (jy2 == XYZ.nY)
		return false;

	cordEmptyInd[0] = ix1;
	cordEmptyInd[1] = ix2;
	cordEmptyInd[2] = jy1;
	cordEmptyInd[3] = jy2;
}

bool Well::findCenterPlace(CoordStorageXYZ XYZ, int cordEmptyInd[4], int kz)
{
	int ix1 = cordEmptyInd[0],
		ix2 = cordEmptyInd[1],
		jy1 = cordEmptyInd[2],
		jy2 = cordEmptyInd[3];

	double xL = center.x - R,
		xR = center.x + R,
		yL = center.y - R,
		yR = center.y + R;

	bool IsAllIndNotFound = true;

	Coord* coords = XYZ.coords;
	int nX = XYZ.nX;
	int nY = XYZ.nY;
	int i, j;
	while (IsAllIndNotFound)
	{
		for (; ix1 > -1 && !is_x_greater(xL, XYZ, ix1, jy1, jy2, kz); ix1--);
		if (ix1 == -1)
		{
			return false;
		}
		IsAllIndNotFound = cordEmptyInd[0] != ix1;

		for (; ix2 < nX && !is_x_less(xR, XYZ, ix2, jy1, jy2, kz); ix2++);
		if (ix2 == nX)
		{
			return false;
		}
			
		IsAllIndNotFound = cordEmptyInd[1] != ix2;
		
		for (; jy1 > -1 && !is_y_greater(yL, XYZ, jy1, ix1, ix2, kz); jy1++);
		if (jy1 == -1)
		{
			return false;
		}
		IsAllIndNotFound = cordEmptyInd[2] != jy1;

		for (; jy2 < nY && !is_y_less(yR, XYZ, jy2, ix1, ix2, kz); jy2++);
		if (jy2 == nY)
		{
			return false;
		}
		IsAllIndNotFound = cordEmptyInd[3] != jy2;
	}
	return true;
}

bool Well::is_x_less(double x, CoordStorageXYZ XYZ, int ix, int jy1, int jy2, int kz)
{
	Coord* coords = XYZ.coords;
	jy2++;
	for (int jy = jy1; jy < jy2; jy++)
		if (coords[XYZ.getKnotIndex(ix, jy, kz)].x > x)
			return false;

	return true;

}
bool Well::is_y_less(double y, CoordStorageXYZ XYZ, int jy, int ix1, int ix2, int kz)
{
	Coord* coords = XYZ.coords;
	ix2++;
	for (int ix = ix1; ix < ix2; ix++)
		if (coords[XYZ.getKnotIndex(ix, jy, kz)].y > y)
			return false;

	return true;
}
bool Well::is_x_greater(double x, CoordStorageXYZ XYZ, int ix, int jy1, int jy2, int kz)
{
	Coord* coords = XYZ.coords;
	jy2++;
	for (int jy = jy1; jy < jy2; jy++)
		if (coords[XYZ.getKnotIndex(ix, jy, kz)].x < x)
			return false;

	return true;
}

bool Well::is_y_greater(double y, CoordStorageXYZ XYZ, int jy, int ix1, int ix2, int kz)
{
	Coord* coords = XYZ.coords;
	ix2++;
	for (int ix = ix1; ix < ix2; ix++)
		if (coords[XYZ.getKnotIndex(ix, jy, kz)].y < y)
			return false;

	return true;
}

void Well::findWellInds(CoordStorageXYZ XYZ, int cordEmptyInd[4])
{
	findBegApprCenterPlace(XYZ, cordEmptyInd);
	int nZ = XYZ.nZ;
	for (int k = 0; k < nZ; k++)
		findCenterPlace(XYZ, cordEmptyInd, k);
}

bool Well::write(std::string fileName)
{
	return true;
}

int Well::getKnotIndex(int iCircleCoord, int jCircle, int kSurf)
{
	return kSurf * nSirfCoords + nCircleCoords * jCircle + iCircleCoord;
}

int Well::buildWellGrid(CoordStorageXYZ XYZ, SepParam sepParam, AreaStorage& areaStorage, ContainerIXYZW IXYZW, int coordResStorage)
{
	findWellInds(XYZ, cordEmptyInd);

	int ix1 = cordEmptyInd[0],
		ix2 = cordEmptyInd[1],
		jy1 = cordEmptyInd[2],
		jy2 = cordEmptyInd[3];

	int dix = (ix2 - ix1) / 2,
		djy = (jy2 - jy1) / 2;

	sepParam.n = dix > sepParam.n ? dix : sepParam.n;
	sepParam.n = djy > sepParam.n ? djy : sepParam.n;

	nCircle = sepParam.n + 1;
	nSirfs = XYZ.nZ;
	nCircleCoords = 2 * (ix2 - ix1 + jy2 - jy1);
	nSirfCoords = nCircleCoords * nCircle;
	count = nCoords = nSirfCoords * nSirfs;

	wellAreaStore.init(cordEmptyInd, nSirfs, areaStorage, IXYZW);

	coords = new Coord[nCoords];
	globI = new int[nCoords];

	int i, j, k;
	int ix, jy;

	Coord* coordsXYZ = XYZ.coords;
	Coord coordBuf, dCoord;

	double qKoeff;
	int nSep = sepParam.n;
	double koef;

	int indLoc, indGlob;
	double z;
	int iLastCircle = nCircle - 1;
	int iGlobAdd = coordResStorage;

	for (k = 0; k < nSirfs; k++)
	{
		i = 0;
		//Определяем координату z для k плоскости скважины
		z = 0;
		z += coordsXYZ[XYZ.getKnotIndex(ix1, jy1, k)].z;
		z += coordsXYZ[XYZ.getKnotIndex(ix2, jy1, k)].z;
		z += coordsXYZ[XYZ.getKnotIndex(ix1, jy2, k)].z;
		z += coordsXYZ[XYZ.getKnotIndex(ix2, jy2, k)].z;
		z /= 4;

		for (ix = ix1; ix < ix2; ix++, i++)
			iGlobAdd = fillLineFromPointToRadius(XYZ, i, k, ix, jy1, z, iGlobAdd, sepParam);
		

		for(jy = jy1; jy < jy2; jy++, i++)
			iGlobAdd = fillLineFromPointToRadius(XYZ, i, k, ix2, jy, z, iGlobAdd, sepParam);

		for (ix = ix2; ix > ix1; ix--, i++)
			iGlobAdd = fillLineFromPointToRadius(XYZ, i, k, ix, jy2, z, iGlobAdd, sepParam);

		for (jy = jy2; jy > jy1; jy--, i++)
			iGlobAdd = fillLineFromPointToRadius(XYZ, i, k, ix1, jy, z, iGlobAdd, sepParam);
		/*
		for (i = 0; i < nCircleCoords; i++)
		{
			globI[getKnotIndex(i, 0, k)] = XYZ.
		}
		*/
	}
	Area area;
	area.lXW = ix1 + 1;
	area.lYW = jy1 + 1;
	area.lZW = 0;
	area.rXW = ix2 - 1;
	area.rYW = jy2 - 1;
	area.rZW = nSirfs - 1;
	int* indKnots;
	int nIndKnots;
	
	XYZ.getAreaKnotsNums(area, indKnots, nIndKnots);
	int iCoord, jCircle, kSurf;
	int nChanged;
	iGlobAdd = coordResStorage;
	//Можно улучшить, будет время можно будет посмотреть
	for (kSurf = 0, nChanged = 0; kSurf < nSirfs; kSurf++)
	{
		for (jCircle = 1; jCircle < nCircle; jCircle++)
		{
			for (iCoord = 0; iCoord < nCircleCoords; iCoord++)
			{
				if (nChanged < nIndKnots)
				{
					globI[getKnotIndex(iCoord, jCircle, kSurf)] = indKnots[nChanged];
					nChanged++;
				}
				else
				{
					globI[getKnotIndex(iCoord, jCircle, kSurf)] = iGlobAdd;
					iGlobAdd++;
				}
									 
			}
		}
	}
	
	/*
	std::ofstream out;
	out.open("CheckFile.txt");
	int max = 0;
	for (i = 0; i < nCoords; i++)
	{
		if (max < globI[i]) max = globI[i];
		out << globI[i] << std::endl;
	}
	out << "Max = " << max << " Count = " << iGlobAdd - nChanged << std::endl;
	*/
	return iGlobAdd;
}

int Well::fillLineFromPointToRadius(CoordStorageXYZ XYZ, int i, int k, int iGlob, int jGlobe, double z, int iGlobAdd, SepParam sepParam)
{
	int indLoc, indGlob;
	int iLastCircle = nCircle - 1;
	indLoc = getKnotIndex(i, 0, k);
	indGlob = XYZ.getKnotIndex(iGlob, jGlobe, k);
	globI[indLoc] = indGlob;
	Coord coordOuter = XYZ.coords[indGlob];
	coords[indLoc] = coordOuter;

	Coord coordBuf = coordOuter;

	coordBuf.z = z;
	Coord coordInner = findCrossPcircle(coordBuf);
	indLoc = getKnotIndex(i, iLastCircle, k);
	coords[indLoc] = coordInner;
	globI[indLoc] = iGlobAdd;
	iGlobAdd++;

	iGlobAdd = fillCoordLine(coordOuter, coordInner, sepParam, i, k, iGlobAdd);
	return iGlobAdd;
}

int Well::fillCoordLine(Coord coordOuter, Coord coordInner, SepParam sepParam, int iCoord, int iSurf, int iGlobAdd)
{
	double qKoeff = sepParam.q;
	int nSep = sepParam.n;
	double koef;
	int iLastCircle = nCircle - 1;
	Coord dCoord;
	int indLoc;
	if (qKoeff == 1)
	{
		dCoord = (coordInner - coordOuter) / nSep;

	}
	else
	{
		koef = (1 - qKoeff) / (1 - pow(qKoeff, nSep));
		dCoord = (coordInner - coordOuter) * koef;
	}

	Coord coordBuf = coordOuter;
	for (int jCircle = 1; jCircle < iLastCircle; jCircle++)
	{
		indLoc = getKnotIndex(iCoord, jCircle, iSurf);
		coordBuf += dCoord;

		coords[indLoc] = coordBuf;
		globI[indLoc] = iGlobAdd;
		iGlobAdd++;

		dCoord *= qKoeff;
	}

	return iGlobAdd;
}

Coord Well::findCrossPcircle(Coord coordBuf)
{
	double x1, x2;
	double y1, y2;
	Coord coordAns(0, 0, coordBuf.z);
	double err = center.x != 0 ? abs((coordBuf.x - center.x) / center.x) : abs(coordBuf.x - center.x);
	if (err < 1e-14)
	{
		y1 = center.y - r;
		y2 = center.y + r;
		coordAns.x = center.x;
		coordAns.y = abs(coordBuf.y - y1) > abs(coordBuf.y - y2) ? y2 : y1;
		return coordAns;
	}

	double a = (coordBuf.y - center.y) / (coordBuf.x - center.x);
	double b = center.y - center.x * a;
	
	auto nx = QuadrEqu::solve(1.0 + a * a, 2 * a * b, b * b - r * r, x1, x2);
	
	coordAns.x = x1;
	if (nx == 2)
		coordAns.x = abs(x1 - coordBuf.x) < abs(x2 - coordBuf.x) ? x1 : x2;
	

	coordAns.y = a * coordAns.x + b;

	return coordAns;
}

bool Well::fillFinitElems(FinitElementStore& finitElementStore)
{
	int iLestCircleCoords = nCircleCoords  - 1,
		jLastCircle = nCircle - 1,
		kLastSurf = nSirfs - 1;
	finitElementStore.nFinitElement = finitElemNum();
	finitElementStore.finitElements = new FinitElement[finitElementStore.nFinitElement];
	FinitElement* finitElements = finitElementStore.finitElements;
	int kSurf, jCircle, iCoord;

	FinitElement finitElement;
	
	WellArea* wellAreas = wellAreaStore.wellAreas;
	WellArea rWellArea;
	int nWellAreas = wellAreaStore.nWellAreas;
	int r, i, j, k, s;
	//int iCoord, jCircle, kSurf;
	int iCoordBeg, iCoordEnd, 
		jCircleBeg, jCircleEnd, 
		kSurfBeg, kSurfEnd;
	jCircleBeg = 0;
	jCircleEnd = nCircle - 1;
	bool fl;
	for (r = 0, s = 0; r < nWellAreas; r++)
	{
		fl = false;
		rWellArea = wellAreas[r];
		iCoordBeg = rWellArea.iCoordBeg;
		iCoordEnd = rWellArea.iCoordEnd;
		kSurfBeg = rWellArea.kSurfBeg;
		kSurfEnd = rWellArea.kSurfEnd;
		if (iCoordEnd == nCircleCoords)
		{
			iCoordEnd--;
			fl = true;
		}
		for (k = kSurfBeg; k < kSurfEnd; k++)
		{
			for (j = jCircleBeg; j < jCircleEnd; j++)
			{
				for (i = iCoordBeg; i < iCoordEnd; i++, s++)
				{
					finitElement.K = rWellArea.K;
					finitElement.FI = rWellArea.phi;
					finitElement.ver[0] = globI[getKnotIndex(i, j, k)];
					finitElement.ver[1] = globI[getKnotIndex(i + 1, j, k)];
					finitElement.ver[2] = globI[getKnotIndex(i, j + 1, k)];
					finitElement.ver[3] = globI[getKnotIndex(i + 1, j + 1, k)];
					finitElement.ver[4] = globI[getKnotIndex(i, j, k + 1)];
					finitElement.ver[5] = globI[getKnotIndex(i + 1, j, k + 1)];
					finitElement.ver[6] = globI[getKnotIndex(i, j + 1, k + 1)];
					finitElement.ver[7] = globI[getKnotIndex(i + 1, j + 1, k + 1)];
					finitElement.phaseStorage.initStorage(rWellArea.phaseStorage.phases, rWellArea.phaseStorage.count);
					finitElements[s] = finitElement;
				}
				if (fl)
				{
					finitElement.K = rWellArea.K;
					finitElement.FI = rWellArea.phi;
					finitElement.ver[0] = globI[getKnotIndex(i, j, k)];
					finitElement.ver[1] = globI[getKnotIndex(0, j, k)];
					finitElement.ver[2] = globI[getKnotIndex(i, j + 1, k)];
					finitElement.ver[3] = globI[getKnotIndex(0, j + 1, k)];
					finitElement.ver[4] = globI[getKnotIndex(i, j, k + 1)];
					finitElement.ver[5] = globI[getKnotIndex(0, j, k + 1)];
					finitElement.ver[6] = globI[getKnotIndex(i, j + 1, k + 1)];
					finitElement.ver[7] = globI[getKnotIndex(0, j + 1, k + 1)];
					finitElement.phaseStorage.initStorage(rWellArea.phaseStorage.phases, rWellArea.phaseStorage.count);
					finitElements[s] = finitElement;
					s++;
				}
			}
		}
	}
	return true;
}

int Well::finitElemNum()
{
	return (nCircleCoords) * (nCircle - 1) * (nSirfs - 1);
}

bool Well::readData(std::ifstream& in)
{
	in >> center.x >> center.y >> center.z >> r >> R;
	return true;
}

Area Well::getXYZMeshInd()
{
	Area area;
	area.lXW = cordEmptyInd[0];
	area.rXW = cordEmptyInd[1];
	area.lYW = cordEmptyInd[2];
	area.rYW = cordEmptyInd[3];
	area.lZW = 0;
	area.rZW = nSirfs - 1;
	return area;
}

int Well::getNumWellFaces()
{
	return (nCircleCoords) * (nSirfs - 1);
}

int Well::getNumWellKnots()
{
	return (nCircleCoords) * (nSirfs);
}

bool Well::read(std::ifstream& in, std::string path, int count)
{
	std::string read_file_err = "Error of reading file " + path;

	if (!in.is_open())
	{
		programlog::writeErr(read_file_err + " - Cannot open the file");
		return false;
	}

	if (count == -1)
	{
		in >> nCircleCoords >> nCircle >> nSirfs;
		if (in.fail())
		{
			programlog::writeErr(read_file_err + " - Error of reading elements count");
			return false;
		}
	}

	if (nCircleCoords < 0)
	{
		programlog::writeErr(read_file_err + " - coordinates count on circle coordinates cannot be < 0");
		return false;
	}

	if (nCircle < 0)
	{
		programlog::writeErr(read_file_err + " - coordinates count on circles cannot be < 0");
		return false;
	}

	if (nSirfs < 0)
	{
		programlog::writeErr(read_file_err + " - coordinates count on sirfs cannot be < 0");
		return false;
	}

	this->count = count = nCoords = nCircleCoords * nCircle * nSirfs;
	nSirfCoords = nCircleCoords * nCircle;
	//nKnots = nZ * knotsInPlane;

	return Storage::read(in, path, count);
	return true;
}
//---------------WellAreaStore-------------//
void WellAreaStore::init(int cordEmptyInd[4], int kz, AreaStorage& AreaStorage, ContainerIXYZW IXYZW)
{
	Area* areas = AreaStorage.areas;
	int nAreas = AreaStorage.count;

	int ix1 = cordEmptyInd[0],
		ix2 = cordEmptyInd[1],
		jy1 = cordEmptyInd[2],
		jy2 = cordEmptyInd[3];

	int* IX_w = IXYZW.IXW.coordInds;
	int* IY_w = IXYZW.IYW.coordInds;
	int* IZ_w = IXYZW.IZW.coordInds;

	int nCoordy = jy2 - jy1;
	int nCoordx = ix2 - ix1;

	int nWellAreasBuf = 2 * (nCoordx + nCoordy) * (kz - 1);

	WellArea* wellAreasBuf = new WellArea[nWellAreasBuf];
	AreaPhaseStorage areaPhaseStorage;
	int niCord = 2*(nCoordx + nCoordy);
	nWellAreas = 0;

	Area area;
	WellArea wellArea1;
	WellArea wellArea2;
	int ix_beg, ix_end, jy_beg, jy_end, kz_beg, kz_end;
	bool fl;
	for (int r = 0; r < nAreas; r++)
	{
		area = areas[r];
		ix_beg = IX_w[area.lXW];
		ix_end = IX_w[area.rXW];
		jy_beg = IY_w[area.lYW];
		jy_end = IY_w[area.rYW];
		kz_beg = IZ_w[area.lZW];
		kz_end = IZ_w[area.rZW];

		fl = ix_beg > ix2 || ix_end < ix1 || jy_beg > jy2 || jy_end < jy1;
		if (fl) continue;

		ix_beg = (ix_beg < ix1) ? ix1 : ix_beg;
		ix_end = (ix_end > ix2) ? ix2 : ix_end;
		jy_beg = (jy_beg < jy1) ? jy1 : jy_beg;
		jy_end = (jy_end > jy2) ? jy2 : jy_end;
		wellArea1.kSurfBeg = kz_beg;
		wellArea1.kSurfEnd = kz_end;
		wellArea1.K = area.lambda;
		wellArea1.phi = area.gamma;

		if (ix_beg == ix1)
		{
			wellArea1.iCoordBeg = niCord - (jy_end - jy1);
			wellArea1.iCoordEnd = niCord - (jy_beg - jy1);
			wellAreasBuf[nWellAreas] = wellArea1;
			wellAreasBuf[nWellAreas].phaseStorage.initStorage(area.phaseStorage.phases, area.phaseStorage.count);
			nWellAreas++;
		}

		if (ix_end == ix2)
		{
			wellArea1.iCoordBeg = nCoordx + (jy_beg - jy1);
			wellArea1.iCoordEnd = nCoordx + (jy_end - jy1);
			wellAreasBuf[nWellAreas] = wellArea1;
			wellAreasBuf[nWellAreas].phaseStorage.initStorage(area.phaseStorage.phases, area.phaseStorage.count);
			nWellAreas++;
		}

		if (jy_beg == jy1)
		{
			wellArea1.iCoordBeg = ix_beg - ix1;
			wellArea1.iCoordEnd = ix_end - ix1;
			wellAreasBuf[nWellAreas] = wellArea1;
			wellAreasBuf[nWellAreas].phaseStorage.initStorage(area.phaseStorage.phases, area.phaseStorage.count);
			nWellAreas++;
		}

		if (jy_end == jy2)
		{
			wellArea1.iCoordBeg = niCord - nCoordy - (ix_end - ix1);
			wellArea1.iCoordEnd = niCord - nCoordy - (ix_beg - ix1);
			wellAreasBuf[nWellAreas] = wellArea1;
			wellAreasBuf[nWellAreas].phaseStorage.initStorage(area.phaseStorage.phases, area.phaseStorage.count);
			nWellAreas++;
		}
	}

	wellAreas = new WellArea[nWellAreas];
	for (int i = 0; i < nWellAreas; i++)
		wellAreas[i] = wellAreasBuf[i];
	
	delete[] wellAreasBuf;
}

bool WellStorage::allocateMemory()
{
	wells = new Well[count];
	return true;
}

bool WellStorage::readData(std::ifstream& in, int i)
{
	wells[i].readData(in);
	return true;
}

bool WellStorage::insertWells(CoordStorageXYZ XYZ, SepParam* sepParams, AreaStorage& areaStorage, ContainerIXYZW IXYZW)
{
	int resCoords;
	resCoords = XYZ.count;
	int niLast = count;
	for (int i = 0; i < niLast; i++)
		resCoords = wells[i].buildWellGrid(XYZ, sepParams[i], areaStorage, IXYZW, resCoords);

	//iCoordMax = wells[niLast].buildWellGrid(XYZ, sepParams[niLast], areaStorage, IXYZW, resCoords);
	iCoordMax = resCoords;
	return true;
}

//---------------FinitElementStore-------------//
void FinitElementStore::clear()
{
	if (nFinitElement > 0)
		delete[] finitElements;
	nFinitElement = 0;
}

double FinitElementStore::integrFunct(Coord eCoord)
{
	ECoord p;
	p = eCoord;

	buildJ(J, p);
	double DetJ = countDetMatrix3(J);

	return DetJ;
}

void FinitElementStore::calcSquares(Coord* coords)
{
	Coord left(0, 0, 0);
	Coord right(1, 1, 1);
	int i, j;
	for (i = 0; i < nFinitElement; i++)
	{
		for (j = 0; j < N; j++)
			this->coords[j] = coords[finitElements[i].ver[j]];
		finitElements[i].sqare = integrateGausse(left, right);
	}
}

void FinitElementStore::write(std::ofstream &out)
{
	out << "____________________________________________________________________________________________________________________________________________________________________________" << std::endl;
	out << "|N of Finit Element|-------Ver 1------|-------Ver 2------|-------Ver 3------|-------Ver 4------|-------Ver 5------|-------Ver 6------|-------Ver 7------|-------Ver 8------|" << std::endl;
	out << "____________________________________________________________________________________________________________________________________________________________________________" << std::endl;
	//out.width(18);
	int i, j;
	int* verFinElem;
	for (i = 0; i < nFinitElement; i++)
	{
		verFinElem = finitElements[i].ver;
		out << "|" << std::setw(18)  << i << "|";
		for (j = 0; j < VER_NUM; j++)
			out << std::setw(18) << verFinElem[j] << "|";
		out << std::endl;
	}
	out << "____________________________________________________________________________________________________________________________________________________________________________" << std::endl;
}

//---------------BorderStorageWell-------------//
bool BorderStorageWell::allocateMemory()
{
	iWells = new int[count];
	nWells = count;
	return true;
}

bool BorderStorageWell::readData(std::ifstream& in, int i)
{
	in >> iWells[i];
	return true;
}

void BorderFacesStore::copyBorderFacesStore(BorderFacesStore borderFacesStore, int nFaces)
{
	if(nFaces > borderFacesStore.nFaces)
		programlog::writeErr("Error of trying copy more elements from the store than there are");

	iFaces = new int[nFaces];
	int* iFaceOther = borderFacesStore.iFaces;
	this->nFaces = nFaces;
	for (int i = 0; i < nFaces; i++)
		iFaces[i] = iFaceOther[i];
}