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
	std::cout << std::noshowpos;

}

// Pretty print a general 2D sparse Eigen array
template <typename type>
void prettyPrint(std::string pre, Eigen::SparseMatrix<type> arr) {

	// Extract the dense array and then call the routine as normal
	prettyPrint(pre, Eigen::Matrix<type,-1,-1>(arr));

}

// Standard cpp entry point
int main (int argc, char ** argv) {

	// Load the data from file
	std::cout << "Loading file..." << std::endl;
	std::string fileName = "matrices/d" + std::string(argv[1]) + "n" + std::string(argv[2]) + "k" + std::string(argv[3]) + ".csv"; 
	std::cout << fileName << std::endl;
	std::ifstream infile(fileName);
	std::vector<Eigen::Triplet<float>> coefficients;
	float val;
	int x, y;
	int arrayWidth = 0;
	int arrayHeight = 0;
	int numOrig = 0;
	int numVars = 0;
	infile >> arrayHeight >> arrayWidth >> numOrig >> numVars;
	while (infile >> x >> y >> val) {
		coefficients.push_back(Eigen::Triplet<float>(x, y, val));
	}
	std::cout << "Matrix is " << arrayHeight << " by " << arrayWidth << " (" << coefficients.size() << " non-zero elements)" << std::endl;

	// Create a sparse Eigen array
	Eigen::SparseMatrix<float> A(arrayHeight, arrayWidth);
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	// Do it without sparsity for now
	Eigen::MatrixXf B = Eigen::MatrixXf(A);

	//prettyPrint("before = ", B);

	// Perform row reduction without swaps TODO
	std::cout << "Solving..." << std::endl;
	std::vector<bool> rowAllowed(arrayHeight, true);
	float factor = 1;
	for (int i=0; i<arrayWidth; i++) {

		// Pick the last allowed non-zero value
		int pivot = -1;
		for (int k=arrayHeight-1; k>0; k--) {
			if (B(k,i) != 0 && rowAllowed[k]) {
				pivot = k;
				rowAllowed[k] = false;
				break;
			}
		}

		// If there is a possible pivot
		if (pivot >= 0) {
			for (int k=arrayHeight-1; k>0; k--) {
				if (k != pivot && B(k,i) != 0) {
					factor = -B(k,i) / B(pivot,i);
					B.row(k) += factor*B.row(pivot);
				}
			}
		}

	}

	// Count the rows with all zero
	int numZero = 0;
	for (int k=0; k<numOrig; k++) {
		if (B.row(k).isZero(0)) {
			numZero += 1;
		}
	}

	//prettyPrint("after = ", B);
	std::cout << "num zero = " << numZero << std::endl;
	std::cout << "num indep eqns = " << numOrig-numZero << std::endl;
	std::cout << "num vars = " << numVars << std::endl;

}
