#include <fstream>
#include <iostream>
#include <iomanip>
#include "programlog.h"
#include "timeinfo.h"

const std::string LogFileName = "ProgramLog.txt";

namespace programlog
{
	const std::string ProgressBar::EMPTY = " ";
	const char ProgressBar::x = -37;

	void ProgressBar::begin()
	{
		oldLocale = setlocale(LC_ALL, NULL);
		setlocale(LC_ALL, "C");
		//std::cout << std::cout.fixed;
		//std::cout.setf(std::cout.fixed);
		//std::cout.precision(5);

		//std::cout << "ggdgsdgggfgfhfh\b\b\b";
		//std::cout << std::endl;
		//getchar();
		int i;
		std::cout << "[";
		for (i = 0; i < nAdd; i++)
			std::cout << " ";

		std::cout << "] ";
		std::cout.setf(std::cout.fixed);
		std::cout.precision(2);
		std::cout.width(7);
		std::cout << 0.0;
		std::cout << " %";
		//getchar();
	}
	void ProgressBar::returnInBeg()
	{
		for (int i = 0; i < nAdd + 11; i++)
			std::cout << "\b";
		//getchar();
	}
	void ProgressBar::showProgress(double perc)
	{
		returnInBeg();
		//std::cout << "\r[";
		int iCur = perc >= 100.0 ? nAdd : perc / hAdd;
		if (perc > 999.99) perc = 999.99;
		if (perc < -999.99) perc = -999.99;
		if (iCur < 0) iCur = 0;
		int i;		
		//std::cout << "\b";
		for (i = 0; i < iCur; i++)
			std::cout << x;

		for (; i < nAdd; i++)
			std::cout << " ";

		std::cout << "] ";
		//std::cout.setf(std::cout.fixed);
		//std::cout.precision(2);
		//std::cout.width(6);
		std::cout.width(7);
		std::cout << perc;
		std::cout << " %";
	}

	void operator << (ProgressBar progressBar, const double perc)
	{
		progressBar.showProgress(perc);
	}

	void ProgressBar::end()
	{
		setlocale(LC_ALL, oldLocale.c_str());
	}

	ProgressBar::ProgressBar()
	{
		oldLocale = EMPTY;
	}

	ProgressBar::~ProgressBar()
	{
		if(oldLocale != EMPTY)
			setlocale(LC_ALL, oldLocale.c_str());
	}
	void writeErrInLogFile(std::string error)
	{
		std::ofstream errFile;
		errFile.open(LogFileName, std::ios::app);

		errFile << timeinfo::getCurrentTime() << " [Error] - " << error << std::endl;
	}

	void writeErrInConcole(std::string error)
	{
		std::cout << timeinfo::getCurrentTime() << " [Error] - " << error << std::endl;
	}

	void writeErr(std::string error)
	{
		writeErrInLogFile(error);
		writeErrInConcole(error);
	}
}
