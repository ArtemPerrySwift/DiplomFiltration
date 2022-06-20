#include "storage.h"
#include <sstream>
#include "programlog.h"

using namespace programlog;

const std::string READ_ERR = "Error of reading file";

bool Storage::read(std::string path, int count)
{
	std::string read_file_err = "Error of reading file " + path;
	std::ifstream binFile;
	binFile.open(path);

	if (!binFile.is_open())
	{
		writeErr(read_file_err + " - Cannot open the file");
		return false;
	}

	bool fl = read(binFile, path, count);
	binFile.close();
	return fl;
}

bool Storage::read(std::ifstream &in, std::string path, int count)
{
	std::string read_file_err = "Error of reading file " + path;

	if (count == -1)
	{
		in >> count;
		if (in.fail())
		{
			writeErr(read_file_err + " - Error of reading elements count");
			return false;
		}
	}

	if (count < 0)
	{
		writeErr(read_file_err + " - elements count cannot be < 0");
		return false;
	}
	this->count = count;
	allocateMemory();

	for (int i = 0; i < count; i++)
		if (!readData(in, i))
		{
			std::ostringstream stream_num;
			stream_num << i;
			writeErr(read_file_err + " - данные в строке " + stream_num.str() + " либо повреждены, либо неправильно заданы");
			in.close();
			return false;
		}

	return true;
}

bool Storage::readBin(std::string path, int count)
{
	std::string read_file_err = "Error of reading file " + path;
	std::ifstream binFile;
	binFile.open(path, std::ios::binary);

	if (!binFile.is_open())
	{
		writeErr(read_file_err + " - Cannot open the file");
		return false;
	}

	if (count == -1)
	{
		binFile.read((char*)&(count), sizeof(int));
		if (binFile.gcount() != sizeof(int))
		{
			writeErr(read_file_err + " - Error of reading elements count");
			return false;
		}
	}

	if (count < 0)
	{
		writeErr(read_file_err + " - elements count cannot be < 0");
		return false;
	}

	allocateMemory();

	for (int i = 0; i < count; i++)
		if (!readData(binFile, i))
		{
			std::ostringstream stream_num;
			stream_num << i;
			writeErr(read_file_err +  " - данные в строке " + stream_num.str() + " либо повреждены, либо неправильно заданы");
			binFile.close();
			return false;
		}

	binFile.close();
	return true;
}