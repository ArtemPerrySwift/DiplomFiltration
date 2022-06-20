#include "phase.h"

bool ComponentStorage::allocateMemory()
{
	components = new Component[count];
	return true;
}

bool ComponentStorage::readData(std::ifstream& in, int i)
{
	Component componentBuf;
	in >> componentBuf.molMass;
	components[i] = componentBuf;
	return true;
}

bool PhaseComponentStorage::allocateMemory()
{
	phaseComponents = new PhaseComponent[count];
	return true;
}

bool PhaseComponentStorage::readData(std::ifstream& in, int i)
{
	PhaseComponent phaseComponentBuf;
	in >> phaseComponentBuf.iComponent >> phaseComponentBuf.iComponProp;
	phaseComponents[i] = phaseComponentBuf;
	return true;
}


bool PhaseStorage::allocateMemory()
{
	phases = new Phase[count];
	return true;
}

bool PhaseStorage::readData(std::ifstream& in, int i)
{
	std::ostringstream stream_num;
	stream_num << i;
	Phase phaseBuf;
	in >> phaseBuf.dynamicVisc >> phaseBuf.saturation;
	//phaseBuf.phaseComponentStorage.read("phase" + stream_num.str() + ".txt");
	phases[i] = phaseBuf;
	return true;
}

