#include <iostream>
#include <sstream>
#include "mesh.h"
#include "programlog.h"
//#include "fileparser.h"

using namespace programlog;
//using namespace fileparser;

namespace meshspace
{
	template <class DataType>
	bool MeshStorage::readAny(DataType*& data, int& count, string path, string paramNameInFile, string ofParamName, string fileCutStorageName)
	{
		string fileParamName = path + "inf2tr.dat";
		if (!fileParser::getParam(paramNameInFile, fileParamName, count))
		{
			writeErr(PARAM_READ_ERR + " - количество " + ofParamName + " из файла " + fileParamName);
			return false;
		}

		data = new DataType[count];
		string fileStorageName = path + fileCutStorageName;
		ifstream binfile_in(fileStorageName, ios::binary);

		if (!binfile_in.is_open())
		{
			writeErr(READ_ERR + " - " + ofParamName + " - " + FILEPARSER_OPEN_FILE_ERR + " - " + fileStorageName);
			return false;
		}

		for (int i = 0; i < count; i++)
			if (!readData(binfile_in, i))
			{
				ostringstream stream_num;
				stream_num << i;
				writeErr(READ_ERR + " " + ofParamName + " - данные в строке " + stream_num.str() + " либо повреждены, либо неправильно заданы в файле " + fileStorageName);
				binfile_in.close();
				return false;
			}
		
		binfile_in.close();
		return true;
	}

	bool FinalElementMaterialStorage::readData(ifstream& binfile_in, int i)
	{
		binfile_in.read((char*)&(finalElemMaterials[i]), sizeof(int));
		if (binfile_in.gcount() != sizeof(int)) return false;
		finalElemMaterials[i]--;
		return true;
	}
	bool KnotsStorage::readData(ifstream& binfile_in, int i)
	{
		binfile_in.read((char*)&(knots[i].x), sizeof(double));
		if (binfile_in.gcount() != sizeof(double)) return false;

		binfile_in.read((char*)&(knots[i].y), sizeof(double));
		if (binfile_in.gcount() != sizeof(double)) return false;

		return true;
	}

	bool FinalElementStorage::readData(ifstream& binfile_in, int i)
	{
		FinalElement endElementBuf;
		binfile_in.read((char*)&(endElementBuf.ver1), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;
		binfile_in.read((char*)&(endElementBuf.ver2), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;
		binfile_in.read((char*)&(endElementBuf.ver3), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;
		binfile_in.read((char*)&(endElementBuf.ver4), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;
		binfile_in.read((char*)&(endElementBuf.ign0), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;
		binfile_in.read((char*)&(endElementBuf.ign1), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;

		endElementBuf.ver1--;
		endElementBuf.ver2--;
		endElementBuf.ver3--;
		endElementBuf.ver4--;

		endElements[i] = endElementBuf;

		return true;
	}

	bool KnotsWithFirstConditionStorage::readData(ifstream& binfile_in, int i)
	{
		binfile_in.read((char*)&(IKnots[i]), sizeof(int));
		if (binfile_in.gcount() != sizeof(int)) false;
		IKnots[i]--;
		return true;
	}

	bool MaterialStorage::readData(ifstream& in, int i)
	{
		Material materialBuf;
		in >> materialBuf.num;  if (in.fail()) return false;
		in >> materialBuf.mu;   if (in.fail()) return false;
		in >> materialBuf.toku; if (in.fail()) return false;

		materialBuf.num--;
		materials[i] = materialBuf;

		return true;
	}

	bool FinalElementMaterialStorage::read(string path)
	{
		readAny(finalElemMaterials, ktr, path, "ktr", "номеров материалов", "nvkat2d.dat");
		return true;
	}

	bool KnotsStorage::read(string path)
	{
		readAny(knots, kuzlov, path, "kuzlov", "узлов", "rz.dat");
		return true;
	}

	bool FinalElementStorage::read(string path)
	{
		readAny(endElements, ktr, path, "ktr", "конечных элементов", "nvtr.dat");
		return true;
	}

	bool KnotsWithFirstConditionStorage::read(string path)
	{
		readAny(IKnots, kt1, path, "kt1", "узлов с первыми краевыми условиями", "l1.dat");
		return true;
	}

	bool MaterialStorage::read(string path)
	{
		readAny(materials, nMat, path, "nMat", "материалов", "mat.dat");
		ifstream muFile;
		for (int i = 0; i < nMat; i++)
		{
			ostringstream stream_num;
			stream_num << materials[i].num;
			muFile.open(path + "mu.00" + stream_num.str());
			if (!muFile.is_open())
				materials[i].isMuConst = true;
			else
			{
				muFile.close();
				materials[i].isMuConst = false;
				materials[i].muSplain.init(path + "mu.00" + stream_num.str());
			}
		}
		
		return true;
	}

	Material MaterialStorage::findMaterialByNum(int num)
	{
		for (int i = 0; i < nMat; i++)
			if (materials[i].num == num)
				return materials[i];

		Material materialBuf;
		materialBuf.num = -1;
		return materialBuf;
	}

	bool FilePreparation::prepareMaterialFile(string path)
	{
		string numMatFile = path + "mu";
		int nMat = fileParser::getNumLinesInFile(numMatFile);
		if (nMat < 0)
			writeErr("Ошибка при подготовке файлов - не удалось подсчитать количество материалов из-за проблем с файлом - " + numMatFile);

		fstream paramFileIn;
		paramFileIn.open(path + "inf2tr.dat");
		
		string filetext;
		if (paramFileIn.is_open())
		{
			string line;
			int indBeg;
			while (getline(paramFileIn, line))
			{
				filetext += line + "\n";
			}
			indBeg = filetext.find("nMat=");
			if (indBeg >= 0)
			{
				filetext.erase(indBeg, filetext.size() - indBeg);
			}
			//paramFileIn << filetext;
		}
		paramFileIn.close();

		ofstream paramFileOut;
		paramFileOut.open(path + "inf2tr.dat", ios::trunc);
		paramFileOut << filetext;
		paramFileOut.close();

		ofstream paramFile;
		

		paramFile.open(path + "inf2tr.dat", ios::app);

		paramFile << " nMat= " << nMat;
		paramFile.close();

		string readMatErr = READ_ERR + " материалов - не удалось прочитать";
		Material materialBuf;
		ifstream inMu, inToku;

		string muFileName = path + "mu";
		string tokuFileName = path + "toku";
		string materialFileName = path + "mat.dat";

		inMu.open(muFileName);
		if (!inMu.is_open())
		{
			writeErr(READ_ERR + " - mu - " + FILEPARSER_OPEN_FILE_ERR + " - " + muFileName);
			return false;
		}

		inToku.open(tokuFileName);
		if (!inMu.is_open())
		{
			writeErr(WRITE_ERR + " - " + FILEPARSER_OPEN_FILE_ERR + " - " + muFileName);
			return false;
		}

		ofstream matFile;
		matFile.open(materialFileName);


		for (int i = 0; i < nMat; i++)
		{
			inMu >> materialBuf.num;
			if (inMu.fail())
			{
				writeMaterialErrOfReading(muFileName, "номер", i);
				closeStreams(inMu, inToku, matFile);
				return false;
			}

			inToku >> materialBuf.num;
			if (inToku.fail())
			{
				writeMaterialErrOfReading(tokuFileName, "номер", i);
				closeStreams(inMu, inToku, matFile);
				return false;
			}

			inMu >> materialBuf.mu;
			if (inMu.fail())
			{
				writeMaterialErrOfReading(muFileName, "mu", i);
				closeStreams(inMu, inToku, matFile);
				return false;
			}

			inToku >> materialBuf.toku;
			if (inToku.fail())
			{
				writeMaterialErrOfReading(tokuFileName, "toku", i);
				closeStreams(inMu, inToku, matFile);
				return false;
			}

			matFile << materialBuf.num << " " << materialBuf.mu << " " << materialBuf.toku << endl;
			if (matFile.fail())
			{
				writeErr(WRITE_ERR + " при подготовке файлов - " + FILEPARSER_OPEN_FILE_ERR + " - " + muFileName);
				closeStreams(inMu, inToku, matFile);
				return false;
			}
		}

	}

	void FilePreparation::writeMaterialErrOfReading(string fileName, string errorObj, int iErrorString)
	{
		ostringstream stream_num;
		stream_num << iErrorString;
		writeErr(READ_ERR + " материалов при подготовке файлов - не удалось прочитать " + errorObj + " материала в файле " + fileName + " в строке " + stream_num.str());
	}

	void FilePreparation::closeStreams(ifstream& in1, ifstream& in2, ofstream& out)
	{
		in1.close();
		in2.close();
		out.close();
	}

	//Добавить необходимые проверки!
	Mesh::Mesh(string path)
	{
		FilePreparation::prepareMaterialFile(path);
		finalElementMaterialStorage.read(path);
		finalElementStorage.read(path);
		knotsWithFirstConditionStorage.read(path);
		materialStorage.read(path);
		knotsStorage.read(path);
	}

	Mesh::Mesh() {};
}

