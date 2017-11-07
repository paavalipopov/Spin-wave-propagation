#include "1DPeriodicCrystal.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace Eigen;

SpinWaveProblem1D::SpinWaveProblem1D(int N, int kSteps, int omegaSteps,
									int precision, int P, int searchIterations){
	this->N = N;
	this->P = P;
	this->searchIterations = searchIterations;
	this->precision = precision;
	this->kSteps = kSteps;
	this->omegaSteps = omegaSteps;

	a = 0.1;
	d = 0.001;
	gamma = 8791208;

	H = 500;
	M1x4pi = 1750;
	M2x4pi = 1350;

	omegaH = gamma * H;
	omegaM1 = gamma * M1x4pi;
	omegaM2 = gamma * M2x4pi;
	omegaDelta = (sqrt(omegaH*(omegaH + omegaM1)) - omegaH);

	goThroughGrid();

}

void SpinWaveProblem1D::goThroughGrid() {
	double kMax = M_PI / a * 500;
	double omegaCheck;
	vector<complex<double> > determinants;
	vector<double> suspiciousOmega;


	ofstream fout1;
	ofstream fout2;
	ofstream fout3;

	fout1.open("all results");
	fout1 << "k*a/2pi omega |detF| ln(|detF|) These are results for grid: k: " << kSteps << " points, omega starting: " << omegaSteps << " points, precision = " << precision << endl;

	fout2.open("results");
	fout2 << "k*a/2pi omega |detF| ln(|detF|) These are results for grid: k: " << kSteps << " points, omega starting: " << omegaSteps << " points, precision = " << precision << endl;

	fout3.open("raw mins");
	fout3 << "k*a/2pi omega |detF| ln(|detF|) These are results for grid: k: " << kSteps << " points, omega starting: " << omegaSteps << " points, precision = " << precision << endl;


	for(double k = -kMax; k < kMax; k += kMax/kSteps) {
		cout << "Calculating layer k = " << k << endl;

		for(double omega = omegaH + 1;
		omega < omegaH + omegaDelta*1.07; omega += omegaDelta/omegaSteps) {

			determinants.push_back(findDeterminant(k, omega));

			fout1 << k/b(1) << " " << omega << " " << abs(determinants.back()) << " " << log(abs(determinants.back())) << endl;

		}

		cout << "Looking for mins" << endl;
		for(int i = 1; i < determinants.size()-1; i++) {
			if(isMinimum(determinants, i, 1)) {
				suspiciousOmega.push_back(checkNull(k, omegaH + omegaDelta/omegaSteps*i, omegaH + omegaDelta/omegaSteps*(i+2), kMax, "init", omegaCheck));

				fout3 << k/b(1) << " " << omegaH + omegaDelta/omegaSteps * (i+1) << " " << determinants[i].real() << " " << determinants[i].imag() << " " << abs(determinants[i]) << endl;
			}
		}



		for(int i = 0; i < suspiciousOmega.size(); i++) {
				if(suspiciousOmega[i] != 0) {
					fout2 << k/b(1) << " " << suspiciousOmega[i] << " " << abs(findDeterminant(k, suspiciousOmega[i])) << "\n";

					fout1 << k/b(1) << " " << suspiciousOmega[i] << " " << abs(findDeterminant(k, suspiciousOmega[i])) << " " << log(abs(findDeterminant(k, suspiciousOmega[i]))) << endl;
				}
		}

		cout << "k = " << k << " done" << endl;
		determinants.clear();
		suspiciousOmega.clear();
	}

	fout1.close();
	fout2.close();
	fout3.close();
}

std::complex<double> SpinWaveProblem1D::findDeterminant(double k, double omega) {
	fillMatrix1(k, omega);
	ces.compute(Matrix1);
	fixEigens();
	refillMatrix1("sin");
	complex<double> det1 = Matrix1.determinant();
	refillMatrix1("cos");
	complex<double> det2 = Matrix1.determinant();
	if(abs(det1) > abs(det2))
		return det2;
	else
		return det1;
}

double SpinWaveProblem1D::checkNull(double& k, double omega1, double omega2, double& startingDet, string mode, double& omegaCheck) {
	if(mode == "abs") {
		if((omega2 - omega1) <= omegaDelta/omegaSteps / (pow(P, searchIterations) + P)) {
			return 0;
		}

		vector<double> determinants;
		for(double omega = omega1; omega < omega2 + 0.5 * ((omega2-omega1) / P);
		omega += (omega2-omega1) / P) {

			determinants.push_back(abs(findDeterminant(k, omega)));
		}


		for(int i = 0;  i < P+1; i++) {
			if((determinants[i]/startingDet) < pow(10, -1 * precision)) {
				cout << "Found omega = " << omega1 + (omega2-omega1) / P * (1+i)
				<< " with determinant " << determinants[i] << endl;

				return omega1 + (omega2-omega1) / P * (1+i);
			}
		}

		for(int i = 1;  i < P+1; i++) {
			if(isMinimum(determinants, i)) {
				if(determinants[i-1] < determinants[i+1] &&
					determinants[i-1] <= determinants[0] &&
					determinants[i-1] <= determinants[P]) {

					omegaCheck = checkNull(k, omega1 + (omega2 - omega1) / P * (i-1),
											omega1 + (omega2 - omega1) / P * i,
												startingDet, "abs", omegaCheck);

					if(omegaCheck != 0)
						return omegaCheck;
				}


				else if(determinants[i-1] > determinants[i+1] &&
							determinants[i+1] <= determinants[0] &&
							determinants[i+1] <= determinants[P]) {

					omegaCheck = checkNull(k, omega1 + (omega2 - omega1) / P * i,
												omega1 + (omega2 - omega1) / P * (i+1),
												startingDet, "abs", omegaCheck);

					if(omegaCheck != 0)
						return omegaCheck;
				}
			}
		}
		return 0;
	}

	else if(mode == "init") {
		complex<double> det1 = findDeterminant(k, omega1);
		complex<double> det2 = findDeterminant(k, omega2);
		double det = min(abs(det1), abs(det2));
		cout << "Checking (" << omega1 << ", " << omega2 <<
		"), starting determinant = " << det << endl;

		return checkNull(k, omega1, omega2, det, "abs", omegaCheck);
	}

	else
		return 0;
}

bool SpinWaveProblem1D::isMinimum(vector<complex<double> >& determinants, int i, int depth) {
	if(depth < 1){
		return true;
	}

	if( abs(determinants[i]) < abs(determinants[i+depth]) &&
		abs(determinants[i]) < abs(determinants[i-depth]) ) {

		return isMinimum(determinants, i, depth-1);
	}

	else
		return false;
}

bool SpinWaveProblem1D::isMinimum(vector<double>& determinants, int i) {
	for(int j = 0; j < determinants.size(); j++) {
		if(j != i)
			if(determinants[j] < determinants[i])
				return false;
	}

	return true;
}

double SpinWaveProblem1D::mu1(double omega) {
	return (omegaH * (omegaH + omegaM1) - omega*omega) / (omegaH*omegaH - omega*omega);
}

double SpinWaveProblem1D::mu2(double omega) {
	return (omegaH * (omegaH + omegaM2) - omega*omega) / (omegaH*omegaH - omega*omega);
}

std::complex<double> SpinWaveProblem1D::M(int n, double omega) {
	if (n == 0)
		return (mu1(omega) + mu2(omega)) * 0.5;
	else if(n % 2 == 0)
		return 0;
	else
		return 1i * (mu2(omega) - mu1(omega)) / (M_PI * n);
}

double SpinWaveProblem1D::b(int m) {
	return m * 2.0 * M_PI / a;
}

void SpinWaveProblem1D::fillMatrix1(double k, double omega) {
	Matrix1 = MatrixXcd(2*N + 1, 2*N + 1);

	for(int i = -N; i <= N; i++) {
		for(int j = -N; j <= N; j++) {
			Matrix1(i + N, j + N) =  -1.0 * M(i - j, omega) * (k + b(i)) * (k + b(j));

			if(Debug == 1) {
				if(abs(Matrix1(i + N, j + N).real()) < pow(10, -1*precision))
					Matrix1(i + N, j + N).real() = 0;

				if(abs(Matrix1(i + N, j + N).imag()) < pow(10, -1*precision))
					Matrix1(i + N, j + N).imag() = 0;
			}
		}
	}

//	cout << "Matrix for eigensolve: " << endl << Matrix1 << endl;
//	cin.get();

	return;
}

void SpinWaveProblem1D::fixEigens() {
	if(Debug == 1) {
		eigenValues = ces.eigenvalues();
		for(int i = 0; i < 2*N + 1; i++) {
			eigenValues(i).imag() = 0;
		}

		eigenVectors = ces.eigenvectors();
		for(int i = 0; i < 2*N + 1; i++) {
			for(int j = 0; j < 2*N + 1; j++) {
				if(abs(eigenVectors(i, j).real()) < pow(10, -1*precision))
					eigenVectors(i, j).real() = 0;

				if(abs(eigenVectors(i, j).imag()) < pow(10, -1*precision))
					eigenVectors(i, j).imag() = 0;
			}
		}
	}
	else {
		eigenValues = ces.eigenvalues();
		eigenVectors = ces.eigenvectors();
	}

//	cout << "Its eigenvectors: " << endl << eigenVectors << endl;
//	cout << "Its eigenvalues: " << endl << eigenValues << endl;
//	cin.get();


	return;
}

void SpinWaveProblem1D::refillMatrix1(string mode) {
	Matrix1 = MatrixXcd(2*N + 1, 2*N + 1);

	for(int i = 0; i < 2*N + 1; i++) {
		for(int j = 0; j < 2*N + 1; j++) {
			if(mode == "sin") {
				Matrix1(i, j) = eigenValues(j) * eigenVectors(i, j) * sin( sqrt(eigenValues(j)) * d);
			}
			if(mode == "cos") {
				Matrix1(i, j) = eigenValues(j) * eigenVectors(i, j) * cos( sqrt(eigenValues(j)) * d);
			}


			if(Debug == 1) {
				if(abs(Matrix1(i, j).real()) < pow(10, -1*precision))
					Matrix1(i, j).real() = 0;

				if(abs(Matrix1(i, j).imag()) < pow(10, -1*precision))
					Matrix1(i, j).imag() = 0;
			}
		}
	}


//	cout << "Matrix for determinant: " << endl << Matrix1 << endl;
//	cin.get();

	return;
}


