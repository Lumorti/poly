#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <iomanip>
#include <random>
#include <fstream>
#include <chrono>
#include <math.h>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>

// Standard cpp entry point
int main(int argc, char ** argv) {

	// Get the problem from the args
	int d = 2;
	int n = 2;
	int k = 0;
	int iters = 1000;
	//if (argc > 1) {
		//d = std::stoi(argv[1]);
	//}
	//if (argc > 2) {
		//n = std::stoi(argv[2]);
	//}
	//if (argc > 3) {
		//k = std::stoi(argv[3]);
	//}
	//if (argc > 4) {
		//iters = std::stoi(argv[4]);
	//}

	// Load the data from file
	std::cout << "Loading file..." << std::endl;
	//std::string fileName = "matrices/d" + std::to_string(d) + "n" + std::to_string(n) + "k" + std::to_string(k) + ".csv"; 
	std::string fileName = "cert.csv";
	std::cout << fileName << std::endl;
	std::ifstream infile(fileName);
	std::vector<Eigen::Triplet<float>> coefficients;
	float val;
	int x, y;
	int arrayWidth = 0;
	int arrayHeight = 0;
	int numOrig = 0;
	int numVars = 0;
	//infile >> arrayHeight >> arrayWidth >> numOrig >> numVars;
	while (infile >> x >> y >> val) {
		coefficients.push_back(Eigen::Triplet<float>(x, y, val));
		if (x > arrayHeight-1) {
			arrayHeight = x+1;
		}
		if (y > arrayWidth-1) {
			arrayWidth = y+1;
		}
	}
	std::cout << "Matrix is " << arrayHeight << " by " << arrayWidth << " (" << coefficients.size() << " non-zero elements)" << std::endl;

	// Create a sparse Eigen array
	Eigen::SparseMatrix<float,Eigen::ColMajor> ACols(arrayHeight, arrayWidth);
	ACols.setFromTriplets(coefficients.begin(), coefficients.end());
	ACols = ACols.transpose();

	// Create the b vector
	Eigen::SparseVector<float> bVec(ACols.rows());
	bVec.coeffRef(0) = 1;

	// Set up the solver
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float>> solver;
	solver.setMaxIterations(iters);
	solver.setTolerance(1e-15);

	// Solve the system
	solver.compute(ACols);
	Eigen::VectorXf sol = solver.solve(bVec);

	// Output the results
	std::cout << "num iterations = " << solver.iterations() << std::endl;
	std::cout << "solver error = " << solver.error() << std::endl;
	std::cout << "real error = " << (ACols*sol - bVec).norm() << std::endl;
	Eigen::VectorXf del = ACols*sol - bVec;
	std::cout << "difference in first term = " << std::abs(del[0]) << std::endl;
	//std::cout << std::endl;
	//std::cout << sol << std::endl;
	//std::cout << std::endl;
	//std::cout << ACols*sol << std::endl;
	
	return 0;

}
