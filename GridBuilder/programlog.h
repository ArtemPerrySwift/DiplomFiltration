#pragma once
#include <string>

namespace programlog
{
	class ProgressBar
	{
	private:
		static const std::string EMPTY;
		static const char x;
		const int nAdd = 40;
		const double hAdd = 100.0 / nAdd;		
		int iCurrent;
		std::string oldLocale; // Предыдущий языковой формат
		void returnInBeg();
	public:
		//void operator << (const double perc);
		ProgressBar();
		~ProgressBar();
		void begin();
		void showProgress(double perc);
		void end();
	};
	void writeErrInLogFile(std::string error);
	void writeErrInConcole(std::string error);
	void writeErr(std::string error);
}
