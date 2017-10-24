#ifndef PERIODICCRYSTAL1D
#define PERIODICCRYSTAL1D

#include <vector>
#include <Eigen>

class SpinWaveProblem1D {
public:
	SpinWaveProblem1D(int N, int kStpes, int omegaSteps,
					int precision, int P, int searchIterations);

private:
	int N;
	int kSteps, omegaSteps;
	double omegaDelta;
	double a, d;
	double gamma;
	double H, M1x4pi, M2x4pi;
	double omegaH, omegaM1, omegaM2;
	double mu1(double omega);
	double mu2(double omega);
	std::complex<double> M(int n, double omega);
	double b(int m);

	int Debug = 1;
	int P, searchIterations, precision;



	void fillMatrix1(double k, double omega);
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
	void fixEigens();
	void refillMatrix1();
	void goThroughGrid();
	std::complex<double> findDeterminant(double k, double omega);
	double checkNull(double& k, double omega1, double omega2,
					double& startingDet, std::string mode, double& omegaCheck);

	bool isMinimum(std::vector<std::complex<double> >& determinants, int i, int depth);
	bool isMinimum(std::vector<double>& determinants, int i);

	Eigen::MatrixXcd Matrix1;
	Eigen::MatrixXcd eigenVectors;
	Eigen::VectorXcd eigenValues;
};

#endif // PERIODICCRYSTAL1D
