#define _CRT_SECURE_NO_WARNINGS
#include <chrono>
#include <ctime>
#include "timeinfo.h"

std::string timeinfo::getCurrentTime()
{
	auto chronoCurrentTime = std::chrono::system_clock::now();
	time_t CurrentTime = std::chrono::system_clock::to_time_t(chronoCurrentTime);

	return ctime(&CurrentTime);
}