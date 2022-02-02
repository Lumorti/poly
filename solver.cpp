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
#include <Eigen/../unsupported/Eigen/KroneckerProduct>
#include <Eigen/../unsupported/Eigen/MatrixFunctions>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>

// Pretty print a generic 1D vector
template <typename type> void prettyPrint(std::string pre, std::vector<type> arr) {

	// Used fixed precision
	int precision = 2;
	std::cout << std::fixed << std::showpos << std::setprecision(precision);

	// For the first line, add the pre text
	std::cout << pre << " { ";

	// For the x values, combine them all on one line
	for (int x=0; x<arr.size(); x++) {
		std::cout << arr[x];
		if (x < arr.size()-1) {
			std::cout << ", ";
		}
	}

	// Output the row
	std::cout << "}" << std::endl;

	// Reset things for normal output
	std::cout << std::defaultfloat << std::noshowpos;

}

// Pretty print a general 2D dense Eigen array
template <typename type>
void prettyPrint(std::string pre, Eigen::Matrix<type, -1, -1> arr) {

	// Used fixed precision
	int precision = 2;
	std::cout << std::fixed << std::showpos << std::setprecision(precision);

	// Loop over the array
	std::string rowText;
	for (int y=0; y<arr.rows(); y++) {

		// For the first line, add the pre text
		if (y == 0) {
			rowText = pre + "{";

		// Otherwise pad accordingly
		} else {
			rowText = "";
			while (rowText.length() < pre.length()+1) {
				rowText += " ";
			}
		}

		// Spacing
		std::cout << rowText << " { ";

		// For the x values, combine them all on one line
		for (int x=0; x<arr.cols(); x++) {
			std::cout << std::setw(precision+2) << arr(y,x);
			if (x < arr.cols()-1) {
				std::cout << ", ";
			}
		}

		// Output the row
		std::cout << "}";
		if (y < arr.rows()-1) {
			std::cout << ",";
		} else {
			std::cout << " } ";
		}
		std::cout << std::endl;

	}

	// Reset things for normal output
	std::cout << std::defaultfloat << std::noshowpos;

}

// Pretty print a general 1D Eigen vector
void prettyPrint(std::string pre, Eigen::VectorXf arr) {

	// Extract the dense array and then call the routine as normal
	prettyPrint(pre, Eigen::MatrixXf(arr));

}

// Pretty print a general 2D sparse Eigen array
template <typename type>
void prettyPrint(std::string pre, Eigen::SparseMatrix<type> arr) {

	// Extract the dense array and then call the routine as normal
	prettyPrint(pre, Eigen::Matrix<type,-1,-1>(arr));

}

// Pretty print a general 2D sparse Eigen array
template <typename type>
void prettyPrint(std::string pre, Eigen::SparseMatrix<type, 1> arr) {

	// Extract the dense array and then call the routine as normal
	prettyPrint(pre, Eigen::Matrix<type,-1,-1>(arr));

}

// Standard cpp entry point
int main(int argc, char ** argv) {

	// Get the problem from the args
	int d = 2;
	int n = 2;
	int k = 0;
	int iters = 100;
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}
	if (argc > 3) {
		k = std::stoi(argv[3]);
	}
	if (argc > 4) {
		iters = std::stoi(argv[4]);
	}

	// Load the data from file
	std::cout << "Loading file..." << std::endl;
	std::string fileName = "matrices/d" + std::to_string(d) + "n" + std::to_string(n) + "k" + std::to_string(k) + ".csv"; 
	std::cout << fileName << std::endl;
	std::ifstream infile(fileName);
	std::vector<Eigen::Triplet<float>> coefficients;
	//std::vector<Eigen::Triplet<float>> coefficientsMinusOrig;
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
	//Eigen::SparseMatrix<float,Eigen::RowMajor> ARows(arrayHeight, arrayWidth);
	//ARows.setFromTriplets(coefficients.begin(), coefficients.end());



	//int canRemove = 0;
	//Eigen::SparseMatrix<float,Eigen::ColMajor> blankVec(ACols.rows(), 1);
	//for (int i=0; i<ACols.cols(); i++) {
		//if (ACols.col(i).nonZeros() <= 1) {
			//canRemove++;
			//ACols.col(i) = blankVec.col(0);
		//}
	//}
	//std::cout << "removed " << canRemove << std::endl;
	//std::cout << "Matrix is now " << ACols.rows() << " by " << ACols.cols() << " (" << ACols.nonZeros() << " non-zero elements)" << std::endl;

	// Create augmented matrix and compare ranks
	//Eigen::SparseMatrix<float,Eigen::ColMajor> aug = ACols;
	//aug.conservativeResize(ACols.rows(), ACols.cols()+1);
	//aug.insert(0, ACols.cols()) = 1;
	//aug.makeCompressed();
	//prettyPrint("col = ", ACols);
	//prettyPrint("aug = ", aug);
	//int rank1 = Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int>>(ACols).rank();
	//std::cout << "rank of coefficient matrix = " << rank1 << std::endl;
	//int rank2 = Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int>>(aug).rank();
	//std::cout << "rank of augmented matrix = " << rank2 << std::endl;
	//std::cout << "system is consistent? " << (rank1 < rank2) << std::endl;
	//std::cout << "(if yes, polynomial system not possible)" << std::endl;
	
	// testing https://arxiv.org/pdf/0801.3788.pdf TODO
	// 2 2 1 is apparently incon
	Eigen::SparseVector<float> bVec(ACols.rows());
	bVec.coeffRef(0) = 1;
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float>> solver;
	solver.setMaxIterations(iters);
	solver.setTolerance(1e-15);
	solver.compute(ACols);
	Eigen::VectorXf sol = solver.solve(bVec);
	std::cout << "num iterations = " << solver.iterations() << std::endl;
	std::cout << "solver error = " << solver.error() << std::endl;
	std::cout << "real error = " << (ACols*sol - bVec).norm() << std::endl;
	Eigen::VectorXf del = ACols*sol - bVec;
	std::cout << "difference in first term = " << std::abs(del[0]) << std::endl;
	//std::cout << std::endl;
	//std::cout << sol << std::endl;
	//std::cout << std::endl;
	//std::cout << ACols*sol << std::endl;
	//prettyPrint("", ACols*sol);
	//prettyPrint("", bVec);
	
	//Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int>> test(ACols);
	//Eigen::VectorXf testSol = test.solve(bVec);
	//std::cout << (ACols*testSol - bVec).norm() << std::endl;

	//int squareSize = std::max(arrayHeight, arrayWidth);
	//ACols.conservativeResize(squareSize, squareSize);
	//bVec.conservativeResize(squareSize);
	//Eigen::VectorXf opt = Eigen::VectorXf::Zero(squareSize);
	//float gam = 0.1 / ACols.squaredNorm();
	//for (int i=0; i<1000; i++) {
		//opt += gam * (ACols*opt - bVec);
	//}
	//std::cout << (ACols*opt - bVec).squaredNorm() << std::endl;
	//std::cout << std::endl;
	//std::cout << opt << std::endl;
	//std::cout << ACols*opt << std::endl;

	return 0;

}
