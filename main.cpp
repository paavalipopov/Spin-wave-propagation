#include <iostream>
#include <string>
#include <Eigen>
#include "1DPeriodicCrystal.h"

using namespace Eigen;
using namespace std;

int main(int argc, char* argv[])
{
	SpinWaveProblem1D A = SpinWaveProblem1D(1, 200, 2000, 4, 100, 8);
//	(int N, int kStpes, int omegaSteps, int precision, int P, int searchIterations)

}

