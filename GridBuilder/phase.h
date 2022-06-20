#pragma once
#include "storage.h"
#include <sstream>

struct Component
{
	double molMass;
};

struct ComponentStorage : Storage
{
	Component* components;

private:
	virtual bool allocateMemory() override;
	virtual bool readData(std::ifstream& in, int i) override;
};

struct PhaseComponent
{
	int iComponent; // номер компонента
	double iComponProp; // доля компонента в фазе
};

struct PhaseComponentStorage : Storage
{
	PhaseComponent* phaseComponents;

private:
	virtual bool allocateMemory() override;
	virtual bool readData(std::ifstream& in, int i) override;
};

struct Phase
{
	double saturation;
	double dynamicVisc; // динамическая вязкость
	//PhaseComponentStorage phaseComponentStorage;
	//PhaseComponent* phaseComponents;
	//int PhaseCompCount;
};

struct PhaseStorage : Storage
{
	Phase* phases;
	bool initStorage(Phase* phases, int nPhases)
	{
		count = nPhases;
		this->phases = new Phase[nPhases];
		for (int i = 0; i < nPhases; i++)
			this->phases[i] = phases[i];
		return true;
	}

private:
	virtual bool allocateMemory() override;
	virtual bool readData(std::ifstream& in, int i) override;
};

