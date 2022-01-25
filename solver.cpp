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
	std::cout << std::noshowpos;

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
	std::cout << std::noshowpos;

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
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}
	if (argc > 3) {
		k = std::stoi(argv[3]);
	}

	// Load the data from file
	std::cout << "Loading file..." << std::endl;
	std::string fileName = "matrices/d" + std::to_string(d) + "n" + std::to_string(n) + "k" + std::to_string(k) + ".csv"; 
	std::cout << fileName << std::endl;
	std::ifstream infile(fileName);
	std::vector<Eigen::Triplet<float>> coefficients;
	std::vector<Eigen::Triplet<float>> coefficientsMinusOrig;
	float val;
	int x, y;
	int arrayWidth = 0;
	int arrayHeight = 0;
	int numOrig = 0;
	int numVars = 0;
	infile >> arrayHeight >> arrayWidth >> numOrig >> numVars;
	while (infile >> x >> y >> val) {
		coefficients.push_back(Eigen::Triplet<float>(x, y, val));
		if (x > numOrig) {
			coefficientsMinusOrig.push_back(Eigen::Triplet<float>(x, y, val));
		}
	}
	std::cout << "Matrix is " << arrayHeight << " by " << arrayWidth << " (" << coefficients.size() << " non-zero elements)" << std::endl;

	// Create a sparse Eigen array
	Eigen::SparseMatrix<float,Eigen::ColMajor> ACols(arrayHeight, arrayWidth);
	ACols.setFromTriplets(coefficients.begin(), coefficients.end());
	//Eigen::SparseMatrix<float,Eigen::RowMajor> ARows(arrayHeight, arrayWidth);
	//Eigen::SparseMatrix<float,Eigen::ColMajor> AColsMinusOrig(arrayHeight, arrayWidth);
	//ARows.setFromTriplets(coefficients.begin(), coefficients.end());
	//AColsMinusOrig.setFromTriplets(coefficientsMinusOrig.begin(), coefficientsMinusOrig.end());

	// Create augmented matrix and compare ranks TODO
	//Eigen::SparseMatrix<float,Eigen::ColMajor> aug = ACols;
	//aug.conservativeResize(arrayHeight, arrayWidth+1);
	//aug.insert(0, arrayWidth) = 1;
	//aug.makeCompressed();

	//prettyPrint("", ACols);
	//prettyPrint("", aug);

	// testing https://arxiv.org/pdf/0801.3788.pdf
	// 2 2 1 is apparently incon
	ACols = ACols.transpose();
	Eigen::VectorXf bVec = Eigen::VectorXf::Zero(ACols.rows());
	bVec(0) = 1;
	//int rank1 = Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int>>(ACols).rank();
	//std::cout << "rank of coefficient matrix = " << rank1 << std::endl;
	//int rank2 = Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int>>(aug).rank();
	//std::cout << "rank of augmented matrix = " << rank2 << std::endl;
	//std::cout << "system is consistent? " << (rank1 < rank2) << std::endl;
	//std::cout << "(if yes, polynomial system not possible)" << std::endl;
	//Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int>> solver;;
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float>> solver;
	solver.compute(ACols);
	Eigen::VectorXf sol = solver.solve(bVec);
	std::cout << "num iterations = " << solver.iterations() << std::endl;
	std::cout << "|Ax-b|^2 error = " << solver.error() << std::endl;
	std::cout << (ACols*sol - bVec).squaredNorm() << std::endl;
	//prettyPrint("", ACols*sol);
	//prettyPrint("", bVec);

	return 0;

	// Perform row reduction without swaps TODO
	//std::cout << "Solving..." << std::endl;
	//std::vector<bool> rowAllowed(arrayHeight, true);
	//float factor = 1;
	//for (long int j=0; j<numOrig; j++) {
		//std::cout << j << " / " << numOrig << std::endl;

		//// Get the non-zero elements of this row
		//std::vector<bool> colAllowed(arrayWidth, true);
		//std::vector<int> colsToProcess;
		//for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator itRow(ARows, j); itRow; ++itRow) {
			//colsToProcess.push_back(itRow.col());
		//}

		//// For each non-zero element of this row
		//while (colsToProcess.size() > 0 && ARows.row(j).nonZeros() < 100) {
			//int currentCol = colsToProcess.front();
			//colsToProcess.erase(colsToProcess.begin());
			//if (!colAllowed[currentCol]) {
				//continue;
			//}
			//colAllowed[currentCol] = false;

			//// Try to first pick a pivot from outside the orig vectors
			//int pivotRow = -1;
			//int pivotVal = -1;
			//for (Eigen::SparseMatrix<float>::InnerIterator itCol(ACols, currentCol); itCol; ++itCol) {
				//if (itCol.row() >= numOrig && rowAllowed[itCol.row()]) {
					//pivotRow = itCol.row();
					//pivotVal = itCol.value();
					//break;
				//}
			//}

			//// If can't find something, use one of the orig vectors
			//if (pivotRow == -1) {
				//for (Eigen::SparseMatrix<float>::InnerIterator itCol(ACols, currentCol); itCol; ++itCol) {
					//if (itCol.row() != j) {
						//pivotRow = itCol.row();
						//pivotVal = itCol.value();
						//break;
					//}
				//}
			//}

			//// If there is another possible pivot
			//if (pivotRow >= 0) {

				//// Stop this row from being used as the pivot for a different column
				//rowAllowed[pivotRow] = false;

				//// For all rows where this elemenet in non-zero
				//for (Eigen::SparseMatrix<float>::InnerIterator itCol(ACols, currentCol); itCol; ++itCol) {

					//// Modify all the other rows
					//if (itCol.row() != pivotRow && itCol.row() != j) {
						//ARows.row(itCol.row()) = (ARows.row(itCol.row()) - (itCol.value() / pivotVal) * ARows.row(pivotRow)).pruned(1e-5);

					//// Modify the top row, updating the list of things to do TODO
					//} else if (itCol.row() != pivotRow) {
						//ARows.row(itCol.row()) = (ARows.row(itCol.row()) - (itCol.value() / pivotVal) * ARows.row(pivotRow)).pruned(1e-5);
						//for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator itRow(ARows, pivotRow); itRow; ++itRow) {
							//colsToProcess.push_back(itRow.col());
						//}
					//}

				//}

				//// Recalculate the col array TODO do smaller updates
				//ACols = Eigen::SparseMatrix<float,Eigen::ColMajor>(ARows);

			//}

		//}
		
	//}

	//// Count the rows with all zero
	//int numZero = 0;
	//for (int l=0; l<numOrig; l++) {
		//if (ARows.row(l).squaredNorm() <= 1e-3) {
			//numZero += 1;
		//}
	//}

	////prettyPrint("after = ", B);
	//std::cout << "num zero = " << numZero << std::endl;
	//std::cout << "num indep eqns = " << numOrig-numZero << std::endl;
	//std::cout << "num vars = " << numVars << std::endl;

}
