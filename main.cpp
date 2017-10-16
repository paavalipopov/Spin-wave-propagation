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
	SpinWaveProblem1D A = SpinWaveProblem1D(3, 100, 4000, 1e-2, 500);
}

