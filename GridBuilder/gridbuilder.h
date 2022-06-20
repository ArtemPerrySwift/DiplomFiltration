#pragma once
#include <fstream>
#include "storage.h"
#include "coord.h"
#include "phase.h"
#include "basicfunct3d.h"
#include "gaussintegr.h"

//#include "mesh.h"

//using namespace std;

const int FACES_NUM = 6;
const int VER_NUM = 8;
const int VER_NUM_FACE = 4;


struct Border
{
	int fun_num;
	int vertL, vertR;
	int horizL, horizR;
	int heightL, heightR;

	virtual bool readData(std::ifstream& in);
};

struct BorderStorage : Storage
{
	Border* borders;

private:
	virtual bool allocateMemory();
	virtual bool readData(std::ifstream& in, int i);
};

struct BorderEmptyC : Border 
{
	double sat1;
	double sat2;
	bool readData(std::ifstream& in) override;
};
struct Border1C : Border{};
struct Border2C : Border{};
struct Border3C : Border
{
	double betta;
	bool readData(std::ifstream& in) override;
};

struct Border1Storage : Storage
{
	Border1C* borders;

private:
	virtual bool allocateMemory();
	virtual bool readData(std::ifstream& in, int i);
};

struct Border2Storage : Storage
{
	Border2C* borders;

private:
	virtual bool allocateMemory();
	virtual bool readData(std::ifstream& in, int i);
};

struct Border3Storage : Storage
{
	Border3C* borders;

private:
	virtual bool allocateMemory();
	virtual bool readData(std::ifstream& in, int i);
};

struct BorderEmptyStorage : Storage
{
	BorderEmptyC* borders;

private:
	virtual bool allocateMemory();
	virtual bool readData(std::ifstream& in, int i);
};

struct BorderStorageWell : Storage
{
	int* iWells;
	int nWells;
private:
	virtual bool allocateMemory();
	virtual bool readData(std::ifstream& in, int i);
};

struct ContainerBorders
{
	Border1Storage firstCondStor;
	BorderStorageWell firstCondStorWell;

	Border2Storage secondCondStor;
	BorderStorageWell secondCondStorWell;

	Border3Storage thirdCondStor;
	BorderStorageWell thirdCondStorWell;

	BorderEmptyStorage emptyCondStor;

	bool readBordersCond1(std::ifstream& in, std::string fileName);
	bool readBordersCond2(std::ifstream& in, std::string fileName);
	bool readBordersCond3(std::ifstream& in, std::string fileName);
	bool readBordersCondEmpty(std::ifstream& in, std::string fileName);

	bool readBorders(std::ifstream& in, std::string fileName);
	bool readBorders(std::string fileName);
};

struct Area
{
	double lambda;
	double gamma;
	int lXW, rXW; // Номера левой и правой границы области по X
	int lYW, rYW; // Номера левой и правой границы области по Y
	int lZW, rZW; // Номера левой и правой границы области по Z
	PhaseStorage phaseStorage;
};

struct AreaStorage: Storage
{
	Area* areas;
	int nPhases;
	virtual bool read(std::ifstream& in, std::string path, int count = -1) override;
	AreaStorage();
private:
	virtual bool allocateMemory();
	virtual bool readData(std::ifstream& in, int i);
};

struct AreaPhase
{
	int iArea;
	int iPhases;
	double PhaseSaturation;
};

struct AreaPhaseStorage: Storage
{
	AreaPhase* areaPhases;
private:
	virtual bool allocateMemory();
	virtual bool readData(std::ifstream& in, int i);
};

struct SepParam // параметры разбиения отрезка сетки
{
	int n; // количество отрезков разбиения
	double q; // коэффициент разбиения
};

struct SepParamStorage : Storage
{
	SepParam* sepParams;
	void multiplySeParams(int koeff);
private:
	virtual bool allocateMemory();
	virtual bool readData(std::ifstream& in, int i);
};

struct FinitElement
{
	double K;
	double FI;
	int ver[VER_NUM];
	int faces[FACES_NUM];
	signed char flowSign[FACES_NUM];
	double sqare;
	PhaseStorage phaseStorage;
	bool getFaceLocalNum(int iLocalFace, int faceVer[VER_NUM_FACE]);
	bool getOpposFaceLocalNum(int iLocalFace, int faceVer[VER_NUM_FACE]);
	bool getFaceGlobalNum(int iLocalFace, int faceVer[VER_NUM_FACE]);
};

struct FinitElementStore: protected BasicFunct3D, protected GausseIntegr3D
{
	FinitElement* finitElements;
	int nFinitElement;
	void clear();
	void calcSquares(Coord* coords);
	void write(std::ofstream &out);
protected:
	double integrFunct(Coord p) override;
};

struct CoordStorage : public Storage
{
	Coord* coords;
	bool copyElemsByInd(int* masInd, Coord* coordRes, int nResCoord);
	void write(std::ofstream &out);
protected:
	virtual bool allocateMemory() override;
	virtual bool readData(std::ifstream& in, int i) override;
};

struct CoordStorageSpecif : public CoordStorage
{
	virtual int getKnotIndex(int ix, int jy, int kz) = 0;
	virtual bool write(std::string fileName) = 0;
	virtual bool fillFinitElems(FinitElementStore& finitElementStore) = 0;
	bool read(std::string path, int count = -1) override;
	virtual bool read(std::ifstream& in, std::string path, int count = -1) = 0;
	virtual int finitElemNum() = 0;
	//bool read(std::ifstream& in, std::string path, int count = -1) override;
protected:
	virtual bool allocateMemory() override;
	virtual bool readData(std::ifstream& in, int i) override;
};

struct CoordStorageXYZ : public CoordStorageSpecif
{
	int nX, nY, nZ;

	int knotsInLine;
	int knotsInPlane;
	int nKnots;

	int nKnotsRes;

	bool init(int nX, int nY, int nZ);
	bool write(std::string fileName) override;
	bool fillFinitElems(FinitElementStore& finitElementStore) override;
	int getKnotIndex(int ix, int jy, int kz) override;
	bool read(std::ifstream& in, std::string path, int count = -1) override;
	void getAreaFinElemsNums(Area area, int* &indFinElem, int& nIndFinElems);
	void getAreaKnotsNums(Area area, int*& indFinElem, int& nIndFinElems);
	virtual int finitElemNum() override;
protected:
	//virtual bool allocateMemory() override;
	//virtual bool readData(std::ifstream& in, int i) override;
};

struct ContainerSepParam
{
	SepParamStorage sepX;
	SepParamStorage sepY;
	SepParamStorage sepZ;
	SepParamStorage sepW;
	int nSepX, nSepY, nSepZ, nSepW;

	ContainerSepParam(int nSepX, int nSepY, int nSepZ, int nSepW, std::string fileName);
	ContainerSepParam();

	bool readSepParams(int nSepX, int nSepY, int nSepZ, int nSepW, std::string fileName);
};

struct IndCoordStore
{
	int* coordInds; // massive index of main cooordinates 
	int indCount;

	bool calcCoordInd(SepParamStorage sepParamStorage);
	bool write(char seperator = ' ');
};

struct ContainerIXYZW
{
	IndCoordStore IXW, IYW, IZW;

	ContainerIXYZW(ContainerSepParam containerSepParam);
	ContainerIXYZW();

	void create(ContainerSepParam containerSepParam);
};

struct CrushedMeshCoordStorage : CoordStorageXYZ
{
	ContainerIXYZW IXYZW;

	CrushedMeshCoordStorage();
	CrushedMeshCoordStorage(CoordStorageXYZ XYZW, ContainerSepParam containerSepParam);

	void initDividesStore(CoordStorageXYZ XYZW, ContainerSepParam containerSepParam);
	bool fillFinitElems(FinitElementStore& finitElementStore, AreaStorage& AreaStorage);

private:
	void copyBasePoints(CoordStorageXYZ& XYZW);
	void copyLineOfBasePoints(SepParam sepParam, CoordStorageSpecif& XYZW, int pNext, int pCurr, int pStep = 1);
};

struct Face
{
	int knots[VER_NUM_FACE];
};

struct BorderFacesStore
{
	int* iFaces;
	int nFaces;
	PhaseStorage phaseStorage;
	void copyBorderFacesStore(BorderFacesStore borderFacesStore, int nFaces);
};

struct FaceStore
{
	Face* faces;
	int count;
	bool init(FinitElementStore& finitElementStore, int nCoord);
	bool createFaceStore(FinitElementStore& finitElementStore);
	bool createFaceStore(FinitElementStore& finitElementStore, BorderFacesStore& borderFacesStore);
	int findNeighboringFinitElem(int iFinElem, FinitElement& finitElement, int iLocalFace);
	void copyStore(FaceStore& facesStore, int nFaces);
private:

	int* ig;
	int* jg;
	int n;
};

struct FlowStore
{
	double* flows;
	int count;

	void init(FaceStore faceStore);
};

struct PhaseVolStore
{
	double* PhaseVol;
	int count;

	void init(FaceStore faceStore);
};

struct ContainerPhaseVol
{
	PhaseVolStore* phaseVolStore;
	int count;

	PhaseVolStore phaseVolStoreMix;

	void init(FaceStore faceStore, int nPhases);
};

struct FacePhaseStore
{
	int* iFace;
	Phase* phases;

	void init();
};

struct FinElementFaces
{
	int iFinElem;
	int iFaces[FACES_NUM];
};

struct FinElementFacesStore
{
	FinElementFaces* finElementFaces;
	int count;
};

struct WellArea
{
	int iCoordBeg;
	int iCoordEnd;

	int kSurfBeg;
	int kSurfEnd;

	double K;
	double phi;
	PhaseStorage phaseStorage;

};

struct WellAreaStore
{
	WellArea* wellAreas;
	int nWellAreas;
	void init(int cordEmptyInd[4], int kz, AreaStorage& AreaStorage, ContainerIXYZW IXYZW);
};

struct Well : CoordStorageSpecif
{
	Coord center;
	double r; // радиус скважины
	double R; // радиус области около скважины; 
	//Coord* coords;
	int* globI;

	
	int nCircleCoords; // Количество точек апроксимирующих окружность
	int nSirfCoords; // Количество точек в плоскости сечения скважины
	int nCoords; // Количество точек скважины

	int nCircle;
	int nSirfs;

	int buildWellGrid(CoordStorageXYZ XYZ, SepParam sepParam, AreaStorage& areaStorage, ContainerIXYZW IXYZW, int coordResStorage);
	bool readData(std::ifstream& in);
	bool write(std::string fileName) override;
	bool fillFinitElems(FinitElementStore& finitElementStore) override;
	virtual int finitElemNum() override;
	int getKnotIndex(int iCircleCoord, int jCircle, int kSurf) override;
	bool read(std::ifstream& in, std::string path, int count = -1) override;
	Area getXYZMeshInd();
	int getNumWellFaces();
	int getNumWellKnots();
private:
	int cordEmptyInd[4];
	WellAreaStore wellAreaStore;

	bool findBegApprCenterPlace(CoordStorageXYZ XYZ, int cordEmptyInd[4]);
	bool findCenterPlace(CoordStorageXYZ XYZ, int cordEmptyInd[4], int kz);

	bool is_x_less(double x, CoordStorageXYZ XYZ, int ix, int jy1, int jy2, int kz);
	bool is_y_less(double x, CoordStorageXYZ XYZ, int jy, int ix1, int ix2, int kz);
	bool is_x_greater(double x, CoordStorageXYZ XYZ, int ix, int jy1, int jy2, int kz);
	bool is_y_greater(double x, CoordStorageXYZ XYZ, int jy, int ix1, int ix2, int kz);

	void findWellInds(CoordStorageXYZ XYZ, int cordEmptyInd[4]);
	Coord findCrossPcircle(Coord p);
	int fillCoordLine(Coord coordOuter, Coord coordInner, SepParam sepParam, int iCoord, int iSurf, int iGlobAdd = 0);
	int fillLineFromPointToRadius(CoordStorageXYZ XYZ, int i, int k, int iGlob, int jGlobe, double z, int iGlobAdd, SepParam sepParamm);
};

struct WellStorage : Storage
{
	Well* wells;
	virtual bool readData(std::ifstream& in, int i) override;
	virtual bool allocateMemory() override;
	int iCoordMax;
	bool insertWells(CoordStorageXYZ XYZ, SepParam* sepParams, AreaStorage& areaStorage, ContainerIXYZW IXYZW);
};
/*
struct WellCordStore
{
	Well well;
	Coord* coords;

	bool buildWellGrid(CoordStorageXYZ XYZ);
};
*/
struct FaceHelper
{
	FaceStore& facesStore;
	
};

struct KnotsWithFirstConditionStorage
{
	int* IKnots;
	int kt1;
};

struct CalculationArea // пространственная сетка
{
	//StoreMeshKnots XYZW;
	//CrushedMeshCoordStorage XYZ;
	//AreaStorage areas;
	//AreaStorage areasOut;
	//ContainerBorders borders;
	//CoordStorageXYZ XYZW;
	//KnotsWithFirstConditionStorage knotsEmptyCond;
	//WellStorage wellStorage;
	CoordStorage coordsStore; // Узлы
	FinitElementStore finitElementStore; // Конечные элемента
	KnotsWithFirstConditionStorage knots1Cond;	// Узлы с первыми краевыми
	FlowStore flowStore; // Потоки через границы
	ContainerPhaseVol containerPhaseVol; // Объёмы перетекающих фаз
	FaceStore faceStore; // Грани
	FaceStore Faces2CondStore; // Грани со вторыми краевыми
	BorderFacesStore borderFacesStore; // Грани на границе расчётной области
	double endT; // Время до которого производиться расчёт
	int nPhases;

	CalculationArea(std::string fileNameCord, std::string fileNameSep, std::string fileNameBorders, bool writeInCordFile = false);
	CalculationArea();
private:
	bool fill1CondKnots(CrushedMeshCoordStorage XYZ, ContainerBorders borders, WellStorage wellStorage);
	bool fillCondFaces(CrushedMeshCoordStorage XYZ, ContainerBorders borders, WellStorage wellStorage);
	bool fillGeneralFinElems(CrushedMeshCoordStorage XYZ, WellStorage wellStorage, AreaStorage areas);
	bool fillGeneralCoords(CoordStorage XYZ, WellStorage wellStorage);
};
