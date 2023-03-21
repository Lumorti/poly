#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Get the problem from the args
	std::srand(time(0));
	int d = 2;
	int n = 4;
	std::string task = "infeasible";
	bool useFull = false;
	int maxIters = -1;
	int verbosity = 1;
	int numToSplit = 0;
	int cores = 4;
	double testParam = 43;
	double stabilityTerm = 1e-13;
	float alpha = 0.9;
	bool firstIsComputational = true;
	bool secondIsUniform = true;
	bool firstElementIsOne = true;
	std::string solver = "mosek";
	std::string level = "1f";
	std::string fileName = "";
	for (int i=0; i<argc; i++) {
		std::string arg = argv[i];
		if (arg == "-h") {
			std::cout << " -d [int]    set the dimension" << std::endl;
			std::cout << " -n [int]    set the number of bases" << std::endl;
			std::cout << " -l [str]    set the level for the relaxation e.g. 1+2f,3p" << std::endl;
			std::cout << " -i [int]    set max iterations (-1 for no limit)" << std::endl;
			std::cout << " -t [dbl]    set the test parameter" << std::endl;
			std::cout << " -a [dbl]    set the scaling for the feasible check" << std::endl;
			std::cout << " -b [dbl]    set the stability term to be added to the Hessian" << std::endl;
			std::cout << " -v [int]    set the verbosity level (0,1,2)" << std::endl;
			std::cout << " -p [int]    set the number of vars to split initially" << std::endl;
			std::cout << " -c [int]    set the number of cores to use" << std::endl;
			std::cout << " -o [str]    log points to a csv file" << std::endl;
			std::cout << " -1          don't assume the first basis is the computational" << std::endl;
			std::cout << " -2          don't assume the first vector of the second basis is uniform" << std::endl;
			std::cout << " -3          don't assume the first element of each is one" << std::endl;
			std::cout << " -s          use scs as the SDP solver insead of mosek" << std::endl;
			std::cout << " -w          use whole bases, not partial" << std::endl;
			std::cout << " -f          try to find a feasible point instead of proving infeasiblity" << std::endl;
			return 0;
		} else if (arg == "-d" && i+1 < argc) {
			d = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-n" && i+1 < argc) {
			n = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-i" && i+1 < argc) {
			maxIters = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-p" && i+1 < argc) {
			numToSplit = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-t" && i+1 < argc) {
			testParam = std::stod(argv[i+1]);
			i++;
		} else if (arg == "-b" && i+1 < argc) {
			stabilityTerm = std::stod(argv[i+1]);
			i++;
		} else if (arg == "-a" && i+1 < argc) {
			alpha = std::stod(argv[i+1]);
			i++;
		} else if (arg == "-l" && i+1 < argc) {
			level = argv[i+1];
			i++;
		} else if (arg == "-v" && i+1 < argc) {
			verbosity = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-c" && i+1 < argc) {
			cores = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-o" && i+1 < argc) {
			fileName = argv[i+1];
			i++;
		} else if (arg == "-s") {
			solver = "scs";
		} else if (arg == "-w") {
			useFull = true;
		} else if (arg == "-1") {
			firstIsComputational = false;
		} else if (arg == "-2") {
			secondIsUniform = false;
		} else if (arg == "-3") {
			firstElementIsOne = false;
		} else if (arg == "-f") {
			task = "feasible";
		} else if (arg == "-i") {
			task = "infeasible";
		}

	}

	// Useful quantities
	int numVarsNonConj = n*d*d;
	int numVars = 2*numVarsNonConj+1000;
	int conjDelta = numVarsNonConj;

	// Different "bases"
	std::vector<std::vector<int>> dLimits;
	if (useFull) {
		dLimits.push_back(std::vector<int>(n, d));
	} else {
		if (d == 2) {
			dLimits.push_back({2, 1, 1, 1});
		} else if (d == 3) {
			dLimits.push_back({3, 1, 1, 1, 1});
		} else if (d == 4) {
			//dLimits.push_back({4, 1, 1, 1, 1, 1});
			dLimits.push_back({4, 2, 1, 1, 1, 1});
		} else if (d == 5) {
			//dLimits.push_back({5, 2, 1, 1, 1, 1, 1});
			dLimits.push_back({5, 3, 1, 1, 1, 1, 1});
		} else if (d == 6) {
			//dLimits.push_back({4, 4, 3, 3});
			//dLimits.push_back({4, 4, 4, 3});
			//dLimits.push_back({4, 4, 4, 4});
			dLimits.push_back({5, 5, 5, 5});
			//dLimits.push_back({6, 3, 3, 3});
			//dLimits.push_back({6, 6, 6, 6});
			//dLimits.push_back({4, 4, 3, 3}); // 14
			//dLimits.push_back({4, 4, 4, 3}); // 15
		} else if (d == 7) {

			// 18
			//dLimits.push_back({5, 4, 4, 1, 1, 1, 1, 1, 1});
			//dLimits.push_back({5, 5, 3, 1, 1, 1, 1, 1, 1});
			//dLimits.push_back({6, 5, 2, 1, 1, 1, 1, 1, 1});
			//dLimits.push_back({7, 5, 1, 1, 1, 1, 1, 1, 1});

			// 23
			dLimits.push_back({3, 3, 3, 3, 3, 2, 2, 2, 2});
			dLimits.push_back({4, 4, 4, 4, 3, 1, 1, 1, 1});
			dLimits.push_back({5, 4, 4, 3, 3, 1, 1, 1, 1});
			dLimits.push_back({6, 4, 3, 3, 3, 1, 1, 1, 1});
			dLimits.push_back({6, 4, 4, 3, 2, 1, 1, 1, 1});

			// 24
			dLimits.push_back({4, 4, 4, 4, 4, 1, 1, 1, 1});
			dLimits.push_back({5, 4, 4, 4, 3, 1, 1, 1, 1});
			dLimits.push_back({6, 4, 4, 3, 3, 1, 1, 1, 1});
			dLimits.push_back({7, 4, 3, 3, 3, 1, 1, 1, 1});
			dLimits.push_back({7, 4, 4, 3, 2, 1, 1, 1, 1});

		} else if (d == 8) {
			//dLimits.push_back({8, 5, 1, 1, 1, 1, 1, 1, 1, 1});
			dLimits.push_back({8, 6, 1, 1, 1, 1, 1, 1, 1, 1});
		} else if (d == 9) {
			dLimits.push_back({3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2});
		} else if (d == 10) {
			dLimits.push_back({10, 5, 5, 4}); // 24
			dLimits.push_back({6, 6, 6, 6}); // 24
			dLimits.push_back({7, 6, 6, 6}); // 25
			dLimits.push_back({10, 5, 5, 5}); // 25
		} else if (d == 11) {
			//dLimits.push_back({11, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
			dLimits.push_back({11, 9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
		} else if (d == 12) {
			//dLimits.push_back({12, 6, 6, 5}); // 29
			dLimits.push_back({8, 7, 7, 7}); // 29
			dLimits.push_back({8, 8, 7, 7}); // 30
			//dLimits.push_back({12, 6, 6, 6}); // 30
		} else if (d == 13) {
			//dLimits.push_back({13, 10, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
			dLimits.push_back({13, 11, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
		} else if (d == 14) {
			//dLimits.push_back({14, 7, 7, 6});
			dLimits.push_back({14, 7, 7, 7});
		} else if (d == 15) {
			//dLimits.push_back({15, 7, 7, 7});
			dLimits.push_back({10, 10, 10, 7});
		} else if (d == 16) {
			dLimits.push_back({16, 14, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
		} else if (d == 17) {
			dLimits.push_back({17, 15, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
		} else if (d == 18) {
			dLimits.push_back({18, 9, 9, 9});
		} else if (d == 19) {
			dLimits.push_back({19, 17, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
		} else if (d == 20) {
			dLimits.push_back({20, 10, 10, 10});
		} else {
			dLimits.push_back(std::vector<int>(n, d));
		}

		// If going past infeasibility, pad with zeros
		while (dLimits[0].size() < n) {
			dLimits[0].push_back(0);
		}

	}

	for (int i2=0; i2<dLimits.size(); i2++) {

		// TODO seesaw between PSD matrix and mag cons
		std::vector<int> setSizes = dLimits[i2];
		int numVectorsTotal = 0;
		for (int i=0; i<n; i++) {
			numVectorsTotal += setSizes[i];
		}
		std::cout << setSizes << " " << numVectorsTotal << std::endl;

		Eigen::MatrixXcd X = Eigen::MatrixXcd::Random(numVectorsTotal, numVectorsTotal);
		Eigen::MatrixXd idealMags = Eigen::MatrixXd::Constant(numVectorsTotal, numVectorsTotal, 1.0/std::sqrt(d));
		int deltaInd = 0;
		for (int i=0; i<n; i++) {
			for (int j=0; j<setSizes[i]; j++) {
				for (int k=0; k<setSizes[i]; k++) {
					if (j == k) {
						idealMags(j+deltaInd, k+deltaInd) = 1;
					} else {
						idealMags(j+deltaInd, k+deltaInd) = 0;
					}
				}
			}
			deltaInd += setSizes[i];
		}

		double prevDelta = -1000;
		double lowestError = 10000;
		for (int i=0; i<1000000; i++) {

			// Perform eigenvalue decomposition
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(X);
			Eigen::MatrixXcd eigenvectors = es.eigenvectors();
			Eigen::VectorXcd eigenvalues = es.eigenvalues();
			double minEigen = eigenvalues[0].real();

			// Set negative eigenvalues to zero
			int rank = 0;
			for (int j=0; j<eigenvalues.size(); j++) {
				eigenvalues(j) = std::max(eigenvalues(j).real(), 0.0);
				if (rank >= d) {
					eigenvalues(j) = 0.0;
				} else if (eigenvalues(j).real() > 1e-10) {
					rank++;
				}

			}

			// Reconstruct the semidefinite matrix
			X = eigenvectors * eigenvalues.asDiagonal() * eigenvectors.adjoint();

			// Check for convergence 
			double error = (X.cwiseAbs() - idealMags).norm();
			lowestError = std::min(error, lowestError);
			if (verbosity >= 2) {
				std::cout << i << " " << error << " " << rank << "         \n" << std::flush;
			} else {
				std::cout << i << " " << error << " " << rank << "         \r" << std::flush;
			}
			if (error < 1e-13) {
				break;
			}

			// Correct all the magnitudes
			Eigen::MatrixXcd XCopy = X;
			for (int j=0; j<X.rows(); j++) {
				for (int k=j; k<X.cols(); k++) {
					double angle = std::arg(X(j, k));
					X(j, k) = std::polar(idealMags(j, k), angle);
					X(k, j) = std::polar(idealMags(j, k), -angle);
				}
			}

			// Check for stalling 
			double delta = (X-XCopy).norm();
			if (std::abs((delta - prevDelta) / delta) < 1e-5) {
				X = Eigen::MatrixXcd::Random(numVectorsTotal, numVectorsTotal);
			}
			prevDelta = delta;

		}
		std::cout << std::endl;
		std::cout << std::endl;
		if (verbosity >= 2) {
			std::cout << X << std::endl;
			std::cout << std::endl;
			std::cout << X.cwiseAbs() << std::endl;
			std::cout << std::endl;
		}

		if (n == 4 && setSizes[3] == 5) {

			// Perform SVD on G.
			Eigen::JacobiSVD<Eigen::MatrixXcd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);

			Eigen::VectorXcd S = svd.singularValues();

			// Define the dimensionality of the reduced space.
			int reducedDimension = 6;
			Eigen::MatrixXcd U_reduced = svd.matrixU().leftCols(reducedDimension);
			Eigen::VectorXcd S_reduced = S.head(reducedDimension);
			Eigen::MatrixXcd lowerDimVectors = U_reduced * S_reduced.asDiagonal();

			std::cout << std::endl;
			std::cout << S << std::endl;
			std::cout << std::endl;

			std::cout << std::endl;
			std::cout << lowerDimVectors*lowerDimVectors.transpose() << std::endl;
			std::cout << std::endl;

			// Print the lower-dimensional vectors.
			for (int i=0; i<lowerDimVectors.rows(); i++) {
				std::cout << "Vector " << i+1 << ": " << lowerDimVectors.row(i) << std::endl;
			}


			//Eigen::LLT<Eigen::MatrixXcd> llt(X);

			//Eigen::MatrixXcd L = llt.matrixL();

			//for (int k=0; k<L.cols(); k++) {
				//std::cout << L.col(k).transpose() << std::endl;
			//}
			//std::cout << L << std::endl;

			//Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(X);
			//Eigen::VectorXcd eigenvalues = es.eigenvalues();
			//double minEigen = eigenvalues[0].real();
			//std::cout << minEigen << std::endl;

			//std::cout << std::endl;
			//std::cout << std::endl;

			//Eigen::MatrixXcd mags = Eigen::MatrixXcd::Constant(numVectorsTotal, numVectorsTotal, 0.0);
			//for (int k=0; k<L.cols(); k++) {
				//for (int l=0; l<L.cols(); l++) {
					//mags(k,l) = L.col(k).conj().dot(L.col(l));
				//}
			//}
			//std::cout << mags << std::endl;

		}

	}

	return 0;

	// For each different restriction
	for (int i2=0; i2<dLimits.size(); i2++) {

		// List the bases and the variables indices
		if (verbosity >= 2) {
			std::cout << "---------------------" << std::endl;
			std::cout << "Basis Indices:" << std::endl;
			std::cout << "---------------------" << std::endl;
			int nextInd = 0;
			for (int i=0; i<n; i++) {
				std::cout << std::endl;
				for (int k=0; k<d; k++) {
					std::cout << "{";
					for (int m=0; m<d; m++) {
						if (k < dLimits[i2][i]) {
							if ((m == 0 && firstElementIsOne) || (i == 0 && firstIsComputational) || (secondIsUniform && i == 1 && k == 0)) {
								std::cout << "-" ;
							} else {
								std::cout << nextInd << "+i" << nextInd + conjDelta;
							}
							if (m < d-1) {
								std::cout << ", ";
							}
						}
						nextInd++;
					}
					std::cout << "}";
					if (k < d-1) {
						std::cout << ", ";
					}
				}
				std::cout << std::endl << std::endl;
			}
		}

		// The list of equations to fill
		std::vector<Polynomial<double>> eqns;

		// If the first basis is computational, don't do normal equations
		int startBasis = 0;
		if (firstIsComputational) {
			startBasis = 1;
		}

		// Generate equations, here iterating over all the vectors
		int newVarInd = 2*numVarsNonConj;
		for (int i=startBasis; i<n; i++) {
			for (int k=0; k<dLimits[i2][i]; k++) {

				// We assume the first vector of the second basis is uniform, so skip it
				if (secondIsUniform && i == 1 && k == 0) {
					continue;
				}

				// Second vector (not repeating the first)
				for (int j=i; j<n; j++) {
					for (int l=0; l<dLimits[i2][j]; l++) {

						// Prevent repeat equations
						if (i == j && l <= k) {
							continue;
						}

						// If we're assuming that the first element of each vector is 1
						Polynomial<std::complex<double>> eqn(numVars);
						int startInd = 0;
						if (firstElementIsOne) {
							startInd = 1;
							eqn.addTerm(1.0/d);
						}

						// (a+ib)*(c-id)
						for (int m=startInd; m<d; m++) {
							int var1 = i*d*d + k*d + m;
							int var2 = j*d*d + l*d + m;
							int var3 = var1 + conjDelta;
							int var4 = var2 + conjDelta;
							eqn.addTerm(1, {var1, var2});
							eqn.addTerm(1, {var3, var4});
							eqn.addTerm(1i, {var1, var4});
							eqn.addTerm(-1i, {var2, var3});
						}

						// For the MUB-ness equations
						if (i != j) {

							// Constrain that this should equal a new complex number
							eqn.addTerm(-1, {newVarInd});
							eqn.addTerm(-1i, {newVarInd+1});

							// Constrain that this new complex should have mag 1/sqrt(d)
							Polynomial<double> extraEqn(numVars);
							extraEqn.addTerm(1, {newVarInd,newVarInd});
							extraEqn.addTerm(1, {newVarInd+1,newVarInd+1});
							extraEqn.addTerm(-1.0/d, {});
							eqns.push_back(extraEqn);
							newVarInd += 2;

						}

						// Split into real and imag parts (both should be 0)
						Polynomial<double> eqnReal = std::real<double>(eqn);
						Polynomial<double> eqnImag = std::imag<double>(eqn);

						// Add these equations if they're not empty
						if (eqnReal.size() > 0) {
							eqns.push_back(eqnReal);
						}
						if (eqnImag.size() > 0) {
							eqns.push_back(eqnImag);
						}

					}
				}

				// All should have mag 1/sqrt(d) (since we're setting first basis to the comp)
				if (firstIsComputational) {
					int startInd = 0;
					if (firstElementIsOne) {
						startInd = 1;
					}
					for (int m=startInd; m<dLimits[i2][0]; m++) {
						int var1 = i*d*d + k*d + m;
						int var2 = var1 + conjDelta;
						Polynomial<double> extraEqn(numVars);
						extraEqn.addTerm(1, {var1,var1});
						extraEqn.addTerm(1, {var2,var2});
						extraEqn.addTerm(-1.0/d, {});
						eqns.push_back(extraEqn);
					}
				}

				// Vs the uniform vector (e.g. |1/sqrt(d) * sum of basis| = 1/sqrt(d))
				if (secondIsUniform) {
					Polynomial<std::complex<double>> extraEqn(numVars);
					int startInd = 0;
					if (firstElementIsOne) {
						extraEqn.addTerm(1.0/d);
						startInd = 1;
					}
					for (int m=startInd; m<d; m++) {
						int var1 = i*d*d + k*d + m;
						int var2 = var1 + conjDelta;
						extraEqn.addTerm(1/std::sqrt(d), {var1});
						extraEqn.addTerm(1i/std::sqrt(d), {var2});
					}
					extraEqn = std::conj(extraEqn)*extraEqn;
					eqns.push_back(std::real<double>(extraEqn));
				}

			}
		}
		int ogEqns = eqns.size();

		// Find the symmetries of the problem
		std::vector<std::unordered_map<int,int>> syms;

		// Ordering within each basis
		int base = 2;
		for (int i=1; i<n; i++) {
			for (int j=0; j<dLimits[i2][i]-1; j++) {
				if (i == 1 && j == 0) {
					continue;
				}
				std::unordered_map<int,int> newSym;
				for (int k=1; k<d; k++) {
					int indLeft = i*d*d+j*d+k; 
					int indRight = i*d*d+(j+1)*d+k; 
					newSym[indLeft] = indRight;
					newSym[indLeft+conjDelta] = indRight+conjDelta;
				}
				syms.push_back(newSym);
			}
		}

		// Ordering of the bases
		for (int i=1; i<n-1; i++) {
			int ind1 = 0;
			if (i == 1) {
				ind1 = 1;
			}
			if (ind1 > dLimits[i2][i]-1) {
				continue;
			}
			std::unordered_map<int,int> newSym;
			for (int k=1; k<d; k++) {
				int indLeft = i*d*d+ind1*d+k;
				int indRight = (i+1)*d*d+0*d+k;
				newSym[indLeft] = indRight;
				newSym[indLeft+conjDelta] = indRight+conjDelta;
			}
			syms.push_back(newSym);
		}

		// Ordering of the elements in the vectors
		for (int k=1; k<d-1; k++) {
			std::unordered_map<int,int> newSym;
			for (int i=1; i<n; i++) {
				for (int j=0; j<dLimits[i2][i]; j++) {
					if (i == 1 && j == 0) {
						continue;
					}
					int indLeft = i*d*d+j*d+k;
					int indRight = i*d*d+j*d+k+1;
					newSym[indLeft] = indRight;
					newSym[indLeft+conjDelta] = indRight+conjDelta;
				}
			}
			syms.push_back(newSym);
		}

		// Could also take the conjugate of everything
		std::unordered_map<int,int> newSym;
		for (int i=1; i<n; i++) {
			for (int j=0; j<dLimits[i2][i]; j++) {
				if (i == 1 && j == 0) {
					continue;
				}
				for (int k=1; k<d; k++) {
					int indLeft = i*d*d+j*d+k;
					newSym[indLeft+conjDelta] = indLeft+conjDelta;
				}
			}
		}
		syms.push_back(newSym);

		// Convert these symmetries into constraints with break them
		std::vector<Polynomial<double>> orderingCons;
		for (int i=0; i<syms.size(); i++) {

			// If the indices are the same, this means the sum should be positive
			if (syms[i].begin()->first == syms[i].begin()->second) {
				Polynomial<double> newCon(numVars);
				int pow = 0;
				for (auto const &pair: syms[i]) {
					newCon.addTerm(1, {pair.first});
					pow++;
					break;
				}
				orderingCons.push_back(newCon);

			// Otherwise we're comparing the two sums
			} else {
				Polynomial<double> newCon(numVars);
				int pow = 0;
				for (auto const &pair: syms[i]) {
					newCon.addTerm(std::pow(2, pow), {pair.first});
					newCon.addTerm(-std::pow(2, pow), {pair.second});
					pow++;
				}
				orderingCons.push_back(newCon);
			}

		}

		// Combine these equations into a single object
		PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, orderingCons);

		// Use as few indices as possible
		std::unordered_map<int,int> reducedMap = prob.getMinimalMap();
		prob = prob.replaceWithVariable(reducedMap);
		if (verbosity >= 2) {
			std::cout << "---------------------" << std::endl;
			std::cout << "Index Mapping: " << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << std::endl;
			std::cout << reducedMap << std::endl;
			std::cout << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << "Final Problem: " << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << prob << std::endl;
		}
		int numVectors = 0;
		for (int i=0; i<n; i++) {
			numVectors += dLimits[i2][i];
		}
		std::cout << dLimits[i2] << " vars: " << prob.maxVariables << ", vectors: " << numVectors << std::endl;

		// If told to find a feasible point
		if (task == "feasible") {
			std::cout << std::scientific;
			std::vector<double> x = prob.findFeasibleEqualityPoint(-1, alpha, 1e-12, maxIters, cores, verbosity, 1.0/std::sqrt(d), stabilityTerm);
			double maxVal = -1000;
			for (int i=0; i<prob.conZero.size(); i++) {
				maxVal = std::max(maxVal, std::abs(prob.conZero[i].eval(x)));
			}
			std::cout << "max viol = " << maxVal << std::endl;

		// If told to prove the search space is infeasible
		} else if (task == "infeasible") {

			// If we're using a higher level mat, add higher-order cons
			int ogCons = prob.conZero.size();
			if (level.find("2") != std::string::npos) {
				for (int i=0; i<ogCons; i++) {
					prob.conZero.push_back(prob.conZero[i]*prob.conZero[i]);
				}
			}

			// If we're using a higher level mat, add higher-order cons
			if (level.find("3") != std::string::npos) {
				for (int i=0; i<ogCons; i++) {
					prob.conZero.push_back(prob.conZero[i]*prob.conZero[i]*prob.conZero[i]);
				}
			}

			// If we're using a higher level mat, add higher-order cons
			if (level.find("4") != std::string::npos) {
				for (int i=0; i<ogCons; i++) {
					prob.conZero.push_back(prob.conZero[i]*prob.conZero[i]*prob.conZero[i]*prob.conZero[i]);
				}
			}

			// Try to prove infeasiblity
			if (solver == "mosek") {
				prob.proveInfeasible(maxIters, level, 1.0/std::sqrt(d), fileName, verbosity, numToSplit);
			} else if (solver == "scs") {
				prob.proveInfeasibleSCS(maxIters, level, 1.0/std::sqrt(d), fileName, verbosity, numToSplit);
			}
			
		}

		// Add a newline if we're trying multiple set sizes
		if (i2 < dLimits.size()-1) {
			std::cout << std::endl;
		}

	}

	return 0;

}
	
