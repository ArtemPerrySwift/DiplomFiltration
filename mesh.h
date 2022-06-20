#pragma once
#include <string>
#include <fstream>
#include "coord.h"
//#include "splain.h"
#define MESH_IS_INCLUDE

using namespace std;
namespace meshspace
{
	const int FACES_NUM = 6;
	const int VER_NUM = 8;
	struct MeshStorage
	{	
		virtual bool read(string path) = 0;
	protected:
		virtual bool readData(ifstream& in, int i) = 0;
		template <class DataType>
		bool readAny(DataType*& data, int& count, string path, string paramNameInFile, string paramName, string fileCutStorageNamee);
	};
	
	struct Knot
	{
		double x;
		double y;
	};

	struct KnotsStorage : MeshStorage
	{
		Knot* knots;
		Coord* coords;
		int kuzlov;
		virtual bool read(string path) override;
	private:
		bool readData(ifstream& in, int i) override;
	};

	struct FinalElementMaterialStorage : MeshStorage
	{
		int* finalElemMaterials;
		int ktr;
		bool read(string path) override;

	private:
		bool readData(ifstream& in, int i) override;
	};

	struct FinalElement
	{
		int ver1, ver2, ver3, ver4, ver5, ver6, ver7, ver8;
		int ver[VER_NUM];
		int faces[FACES_NUM];
		int ign0, ign1;
	};

	struct FinalElementStorage : MeshStorage
	{
		FinalElement* endElements;
		int ktr;
		bool read(string path) override;
	private:
		bool readData(ifstream& in, int i) override;
	};

	struct KnotsWithFirstConditionStorage : MeshStorage
	{
		int* IKnots;
		int kt1;

		bool read(string path) override;
	private:
		bool readData(ifstream& in, int i) override;
	};

	struct Material
	{
		int num;
		double mu;
		double toku;
		bool isMuConst;
		//splain muSplain;
	};

	struct MaterialStorage : MeshStorage
	{
		Material* materials;
		int nMat;
		Material findMaterialByNum(int num);
		bool read(string path) override;
	private:
		bool readData(ifstream& in, int i) override;
	};

	class FilePreparation
	{
	public:
		static bool prepareMaterialFile(string path);
	private:
		static void writeMaterialErrOfReading(string fileName, string errorObj, int iErrorString);
		static void closeStreams(ifstream& in1, ifstream& in2, ofstream& out);
	};

	struct Mesh
	{
		FinalElementMaterialStorage finalElementMaterialStorage;
		FinalElementStorage finalElementStorage;
		KnotsWithFirstConditionStorage knotsWithFirstConditionStorage;
		MaterialStorage materialStorage;
		KnotsStorage knotsStorage;

		Mesh();
		Mesh(string path);

	};
}

