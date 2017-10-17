#include <iostream>
#include <string>
#include <Eigen>
#include "1DPeriodicCrystal.h"

using namespace Eigen;
using namespace std;

int main(int argc, char* argv[])
{
//	cin.get();

//	if(argc != 4) {
//		cout << "Not enough arguments" << endl;
//		return -1;
//	}
//

//	SpinWaveProblem1D A = SpinWaveProblem1D(std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]));
	SpinWaveProblem1D A = SpinWaveProblem1D(1, 400, 500, 6, 100, 2);
//	SpinWaveProblem1D(int N, int kStpes, int omegaSteps, int precision, int P, int searchIterations)
}

