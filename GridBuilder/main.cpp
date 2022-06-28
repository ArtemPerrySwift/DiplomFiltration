#include "gridbuilder.h"
#include "fem.h"
#include "flows.h"
#include "satculcer.h"
#include "array.h"
#include "balance.h"
#include <iostream>
#include <string>
#include <direct.h>
#include <iomanip>
#include "DifferentEquParams.h"

const int N_TIMES = 100;
void Output2D(int it, double* q, CalculationArea calculationArea, double t);

int main()
{
	CalculationArea calculationArea("CoordAndAreas.txt", "SepParams.txt", "Borders.txt");
	//std::ofstream out;

	FEM fem;
	
	//fem.init(calculationArea);
	double* q = new double[calculationArea.coordsStore.count];
	arrayspace::fill_vec(q, calculationArea.coordsStore.count, 0);
	//fem.getSolutWeights(q);
	/*
	int qN = calculationArea.XYZ.count;
	/*
	int qN = calculationArea.coordsStore.count;
	
	std::std::cout << "Solution" << std::std::endl << std::std::endl;
	for (int i = 0; i < qN; i++)
	{
		std::std::cout << q[i] << std::std::endl;
	}
	*/
	
	
	
	FlowCulcer flowCulcer;
	SatCulcer satCulcer;

	std::ofstream outPress, outFlow, outSat;
	fem.init(calculationArea);
	outPress.open("Pressure.txt");
	outPress << "Time " << 0.0 << std::endl;
	fem.getSolutWeights(q);
	fem.printSolution(outPress);
	outPress.close();
	outPress << std::endl;
	
	Output2D(0, q, calculationArea, 0);
	//return 0;
	std::cout << setprecision(15) << endl;
	Coord p(2, 0.0, 0.5);
	std::cout << "r               Solution " << std::endl;
	for (double x = 2.0; x < 499.0; x += 10)
	{
		p.x = x;
		std::cout << p.x << " " << p.y << " " << p.z << " " << fem.countSolution(p) /* << " " << DifferentEquParams::u1(p.x, p.y, p.z) << " " << fem.countSolution(p) - DifferentEquParams::u1(p.x, p.y, p.z)*/ << std::endl;
	}
	std::cout << endl;
	
	flowCulcer.init(calculationArea, q);
	flowCulcer.calcFlows();
	outFlow.open("Flows.txt");
	outFlow << "Time " << 0.0 << std::endl;
	flowCulcer.printFlows(outFlow);
	outFlow << std::endl;

	Balance balance;
	balance.init(calculationArea);
	balance.balanceFlows();
	
	outFlow << "Time " << 0.0 << std::endl;
	flowCulcer.printFlows(outFlow);
	outFlow << std::endl;
	outFlow.close();
	//flowCulcer.init(calculationArea, q);
	//double* flows = flowCulcer.getflows();
	//flowCulcer.calcFlows();
	
	double t = 0, endT = calculationArea.endT;
	/*Расчёт насыщенностей*/
	satCulcer.init(calculationArea, endT - t);
	t += satCulcer.reculcSat();
	outSat.open("Saturations.txt");
	outSat << "Time " << t << std::endl;
	satCulcer.printSat(outSat);
	outSat << std::endl;
	outSat.close();
	//return 0;
	FinitElemCulcer finitElemCulcer;
	int nElems = calculationArea.finitElementStore.nFinitElement;
	FinitElement* finitElements = calculationArea.finitElementStore.finitElements;
	bool isElemFound = false;
	double ans;
	//std::cout.scientific;
	//endT = 1.7e-4;
	for(int i = 1; t < endT && i < N_TIMES; i++)
	{
		//std::cout << std::scientific << "Time " << t << " begin to count" << std::endl;
		/*Расчёт давления*/
		//fem.init(calculationArea);
		outPress.open("Pressure.txt", std::ios_base::app);
		outPress << "Time " << t << std::endl;
		fem.printSolution(outPress);
		outPress.close();
		outPress << std::endl;
		fem.getSolutWeights(q);
		
		/*Расчёт потоков через границу*/
		flowCulcer.init(calculationArea, q);
		flowCulcer.calcFlows();
		outFlow.open("Flows.txt", std::ios_base::app);
		outFlow << "Time " << t << std::endl;
		balance.balanceFlows();
		flowCulcer.printFlows(outFlow);
		outFlow << std::endl;
		outFlow.close();

		/*Расчёт насыщенностей*/
		//satCulcer.init(calculationArea, endT - t);
		t += satCulcer.reculcSat();
		outSat.open("Saturations.txt", std::ios_base::app);
		outSat << "Time " << t << std::endl;
		satCulcer.printSat(outSat);
		outSat << std::endl;
		outSat.close();
		//std::cout << std::scientific << "Time " << t << " end to count" << std::endl;
		//if (t > 0.5) return 0;
		for (double x = 0.1; x < 499.0; x += 10)
		{
			p.y = x;
			isElemFound = false;
			for (int i = 0; i < nElems && !isElemFound; i++)
			{
				finitElemCulcer.init(finitElements[i], calculationArea.coordsStore, q);
				ans = finitElemCulcer.countFunct(p, isElemFound);
				if (isElemFound)
				{
					std::cout << x << " " << finitElements[i].phaseStorage.phases[0].saturation << " " << finitElements[i].phaseStorage.phases[1].saturation << std::endl;//return ans;
				}
			}

		}
		Output2D(i, q, calculationArea, t);
	}
	
	for (double x = 0.1; x < 499.0; x += 10)
	{
		p.x = x;
		isElemFound = false;
		for (int i = 0; i < nElems && !isElemFound; i++)
		{
			finitElemCulcer.init(finitElements[i], calculationArea.coordsStore, q);
			ans = finitElemCulcer.countFunct(p, isElemFound);
			if (isElemFound)
			{
				std::cout << p.x << " " << finitElements[i].phaseStorage.phases[0].saturation << " " << finitElements[i].phaseStorage.phases[1].saturation << std::endl;//return ans;
			}
		}
		
	}
	


	//qN = calculationArea.faceStore.count;
	/*
	int nX = calculationArea.XYZ.nX;
	int midX = calculationArea.XYZ.nX / 2;
	int midY = calculationArea.XYZ.nY / 2;
	int midZ = calculationArea.XYZ.nZ / 2;
	Coord* coords = calculationArea.XYZ.coords;
	for (int i = midX; i < nX; i++)
	{
		std::cout << coords[calculationArea.XYZ.getKnotIndex(i, midY, midZ)].x << " " << q[calculationArea.XYZ.getKnotIndex(i, midY, midZ)] << std::endl;
	}
	std::cout << std::endl;
	*/
	/*
	std::cout << "FLOWS" << std::endl << std::endl;
	for (int i = 0; i < qN; i++)
	{
		std::cout << flows[i] << " ";
	}
	*/
	return 0;
}

/*Переопределение операторов для вывода и ввода в бинарный файл*/
__forceinline std::ostream& operator < (std::ostream& file, const double& data)
{
	file.write((char*)&data, sizeof(data));
	return file;
}
__forceinline std::ostream& operator < (std::ostream& file, const int& data)
{
	file.write((char*)&data, sizeof(data));
	return file;
}

__forceinline std::istream& operator > (std::istream& file, double& data)
{
	file.read((char*)&data, sizeof(data));
	return file;
}

__forceinline std::istream& operator > (std::istream& file, int& data)
{
	file.read((char*)&data, sizeof(data));
	return file;
}

void Output2D(int it, double* q, CalculationArea calculationArea, double t)
{
	const int SCALE_PRINT = 1; // Не знаю что это такое
	std::cout << "!!!!!!!!!!!!!Scale = " << SCALE_PRINT << std::endl;
	int kuslov = calculationArea.coordsStore.count; //  Количество узлов сетки
	int kolel = calculationArea.finitElementStore.nFinitElement; // Количество конечных элементов
	int kolVerInel = VER_NUM; // Количество вершин в конечном элементе

	FinitElement* finitElems = calculationArea.finitElementStore.finitElements;
	Coord* coords = calculationArea.coordsStore.coords;
	int i, j;
	const int* nKnotElems = calculationArea.faceStore.getIg();/*calculationArea.faceStore.ig*/;
	double* qPhase1 = new double[kuslov];
	double* qPhase2 = new double[kuslov];
	std::string pathInput = "output2D_temperature";
	if (it == 0)
	{
		/*Создаём директорию для выходных файлов, которые будут использоваться построителем сетки 
		Необходима библиотека #include <direct.h>*/
		system(std::string("rmdir /s /q " + pathInput).c_str());
		_mkdir(pathInput.c_str());


		/*Эту часть можно оставить без изменений*/
		std::string path = pathInput + "/inftry.dat";
		std::ofstream ofp;
		ofp.open(path, std::ios::binary);
		ofp << "\tISLAU=	0 INDKU1=	 0 INDFPO=	1" << std::endl;
		ofp << "KUZLOV= " << kuslov << "  KPAR= " << kolel << "    KT1= 0   KTR2= 0   KTR3= 0" << std::endl;
		ofp << "KISRS1= 0 KISRS2= 0 KISRS3= 0   KBRS= 0" << std::endl;
		ofp << "\tKT7= 0   KT10= 0   KTR4= 0  KTSIM= 0" << std::endl;
		ofp << "\tKT6= 0" << std::endl;
		ofp.close();
		ofp.clear();

		path = pathInput + "/nver.dat";
		ofp.open(path, std::ios::binary);
		/*Заполняем информацию о конечных элементах*/
		for (int i = 0; i < kolel; i++)
		{
			/*Заплняем информацию об i-ом конечном элементе*/
			for (int j = 0; j < kolVerInel; j++)
				ofp < (finitElems[i].ver[j] + 1);
			

			/*Непонятно зачем но без этого вроде как работать не будет, так что это необходиомо оставить*/
			for (int j = 0; j < 6; j++)
				ofp < 2;
			
		}
		ofp.close();
		ofp.clear();

		path = pathInput + "/xyz.dat";
		ofp.open(path, std::ios::binary);
		/*Запрлняем информацию о координатах сетки*/
		for (int i = 0; i < kuslov; i++)
		{
			ofp < coords[i].x;
			ofp < coords[i].y;
			ofp < coords[i].z;
		}

		ofp.close();
		ofp.clear();

		path = pathInput + "/nvkat.dat";
		ofp.open(path, std::ios::binary);
		/*Ставиться номер материала на конечном элементе*/
		for (int i = 0; i < kolel; i++)
		{
			ofp < 1;
		}
		ofp.close();
		ofp.clear();

		path = pathInput + "\\smtr";
		ofp.open(path, std::ios::binary);
		/*Непонятно зачем но без этого вроде как работать не будет, так что это необходиомо оставить*/
		for (int uz = 0; uz < kuslov; uz++)
		{
			ofp < 1;
		}
		ofp.close();
		ofp.clear();


		path = pathInput + "\\fields.cnf";
		ofp.open(path, std::ios::binary);
		/*Непонятно зачем но без этого вроде как работать не будет, так что это необходиомо оставить*/
		ofp << "4" << std::endl;
		ofp << std::endl;
		ofp << "Pressure" << std::endl;
		ofp << "Displacement - X" << std::endl;
		ofp << "Displacement - Y" << std::endl;
		ofp << "Displacement - Z" << std::endl;
		ofp << std::endl;
		ofp << std::endl;

		std::ofstream ouf;
		path = pathInput + "\\times_main";
		/*Заполняем временные слои вроде как (пока что временной слой один)*/
		ouf.open(path);
		ouf << N_TIMES << std::scientific << std::setprecision(15) << std::endl;
		//ouf << 0.0 << std::endl;
		//ouf << 0.1 << std::endl;
		/*
		* На случай нескольких временных слоёв
		ouf << TimeGrid.size() << scientific << setprecision(15) << std::endl;
		for (int i = 0; i < TimeGrid.size(); i++)
			ouf << (double)TimeGrid[i] / (3600.0 * 24.0) << std::endl;
		*/
		ouf.close();
		ouf.clear();
		
	}
	

	arrayspace::fill_vec(qPhase1, kuslov, 0);
	arrayspace::fill_vec(qPhase2, kuslov, 0);

	for (int i = 0; i < kolel; i++)
	{
		/*Заплняем информацию об i-ом конечном элементе*/
		for (int j = 0; j < kolVerInel; j++)
		{
			qPhase1[finitElems[i].ver[j]] += finitElems[i].phaseStorage.phases[0].saturation;
			qPhase2[finitElems[i].ver[j]] += finitElems[i].phaseStorage.phases[1].saturation;
		}

	}

	for (int i = 0; i < kuslov; i++)
	{
		qPhase1[i] /= nKnotElems[i + 1] - nKnotElems[i];
		qPhase2[i] /= nKnotElems[i + 1] - nKnotElems[i];
	}

	std::string path = pathInput + "\\sx." + std::to_string(it);
	std::ofstream ofp3(path, std::ios::binary);
	/*Заполняем значения давления на it временном слое*/
	for (int i = 0; i < kuslov; i++) {
		ofp3 < q[i];
	}
	ofp3.close();
	ofp3.clear();

	path = pathInput + "\\sy." + std::to_string(it);
	ofp3.open(path, std::ios::binary);
	for (int i = 0; i < kuslov; i++) {
		ofp3 < qPhase1[i];
	}
	ofp3.close();
	ofp3.clear();

	path = pathInput + "\\sz." + std::to_string(it);
	ofp3.open(path, std::ios::binary);
	for (int i = 0; i < kuslov; i++) {
		ofp3 < qPhase2[i];
	}
	ofp3.close();
	ofp3.clear();

	std::ofstream ouf;
	path = pathInput + "\\times_main";
	/*Заполняем временные слои вроде как (пока что временной слой один)*/
	ouf.open(path, std::ios::app);
	ouf << t << std::endl;

	/*
	for (int i = 0; i < mesh.koluz; i++) {
		ofp3 < T[i];
	}
	*/
	path = pathInput + "/materials";
	std::ofstream ofp;
	ofp.open(path, std::ios::binary);
	ofp << 1;
	ofp << 1;
	ofp << 1;

	ofp3.close();
	ofp3.clear();
	/*
	path = pathInput + "\\sx." + std::to_string(it + 1);
	ofp3.open(path, std::ios::binary);
	for (int i = 0; i < kuslov; i++) {
		ofp3 < q[i];
	}
	ofp3.close();
	ofp3.clear();
	*/
	delete[] qPhase1;
	delete[] qPhase2;
}