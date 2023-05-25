#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Get the problem from the args
	int d = 2;
	int n = 4;
	int task = 0;
	int maxIters = -1;
	int verbosity = 1;
	int numToSplit = 0;
	int cores = 4;
	int numBits = 4;
	double stabilityTerm = 1e-20;
	float alpha = 0.99;
	double tolerance = 1e-12;
	bool debugFlag = false;
	bool firstIsComputational = true;
	bool secondIsUniform = true;
	bool firstElementIsOne = true;
	std::string givenSizes = "";
	std::string solver = "mosek";
	std::string level = "1f";
	std::string fileName = "";
	for (int i=0; i<argc; i++) {
		std::string arg = argv[i];
		if (arg == "-h") {
			std::cout << " -d [int]    set the dimension" << std::endl;
			std::cout << " -n [int]    set the number of bases" << std::endl;
			std::cout << " -m [int]    change the mode:" << std::endl;
			std::cout << "             0 = try to find a feasible point" << std::endl;
			std::cout << "             1 = try to prove infeasiblity" << std::endl;
			std::cout << "             2 = perform a redundancy analysis" << std::endl;
			std::cout << "             3 = testing binarization" << std::endl;
			std::cout << " -N [str]    set the basis sizes e.g. 2,1,1,1" << std::endl;
			std::cout << " -l [str]    set the level for the relaxation e.g. 1+2f,3p" << std::endl;
			std::cout << " -i [int]    set max iterations (-1 for no limit)" << std::endl;
			std::cout << " -t [dbl]    set the tolerance" << std::endl;
			std::cout << " -a [dbl]    set the scaling for the feasible check" << std::endl;
			std::cout << " -b [dbl]    set the stability term to be added to the Hessian" << std::endl;
			std::cout << " -v [int]    set the verbosity level (0,1,2)" << std::endl;
			std::cout << " -p [int]    set the number of vars to split initially" << std::endl;
			std::cout << " -c [int]    set the number of cores to use" << std::endl;
			std::cout << " -o [str]    log points to a csv file" << std::endl;
			std::cout << " -B [int]    set number of bits to use for the binarization" << std::endl;
			std::cout << " -r          use a random seed" << std::endl;
			std::cout << " -1          don't assume the first basis is the computational" << std::endl;
			std::cout << " -2          don't assume the first vector of the second basis is uniform" << std::endl;
			std::cout << " -3          don't assume the first element of each is one" << std::endl;
			std::cout << " -s          use SCS as the SDP solver insead of Mosek" << std::endl;
			std::cout << " -0          debug flag" << std::endl;
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
		} else if (arg == "-b" && i+1 < argc) {
			stabilityTerm = std::stod(argv[i+1]);
			i++;
		} else if (arg == "-t" && i+1 < argc) {
			tolerance = std::stod(argv[i+1]);
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
		} else if (arg == "-B" && i+1 < argc) {
			numBits = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-c" && i+1 < argc) {
			cores = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-o" && i+1 < argc) {
			fileName = argv[i+1];
			i++;
		} else if (arg == "-N" && i+1 < argc) {
			givenSizes = argv[i+1];
			i++;
		} else if (arg == "-m" && i+1 < argc) {
			task = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-s") {
			solver = "scs";
		} else if (arg == "-r") {
			std::srand(time(0));
		} else if (arg == "-1") {
			firstIsComputational = false;
		} else if (arg == "-2") {
			secondIsUniform = false;
		} else if (arg == "-3") {
			firstElementIsOne = false;
		} else if (arg == "-0") {
			debugFlag = true;
		}

	}

	// Different "bases"
	std::vector<int> basisSizes;
	if (givenSizes.size() > 0) {

		// Load the basis sizes
		std::string currentNum = "";
		for (int i=0; i<givenSizes.size(); i++) {
			if (givenSizes[i] == ',') {
				basisSizes.push_back(std::stoi(currentNum));
				currentNum = "";
			} else {
				currentNum += givenSizes[i];
			}
		}
		basisSizes.push_back(std::stoi(currentNum));
		n = basisSizes.size();

	// Otherwise use sizes based on dimension
	} else {
		if (d == 2) {
			basisSizes = {2, 1, 1, 1};
		} else if (d == 3) {
			basisSizes = {3, 1, 1, 1, 1};
		} else if (d == 4) {
			basisSizes = {4, 2, 1, 1, 1, 1};
		} else if (d == 5) {
			basisSizes = {5, 2, 2, 1, 1, 1, 1};
		} else if (d == 6) {
			basisSizes = {6, 6, 3, 1};
		} else if (d == 7) {
			basisSizes = {7, 2, 2, 2, 1, 1, 1, 1, 1};
		} else if (d == 8) {
			basisSizes = {8, 2, 2, 2, 2, 1, 1, 1, 1, 1};
		} else if (d == 9) {
			basisSizes = {3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2};
		} else if (d == 10) {
			basisSizes = {6, 6, 6, 6};
		} else if (d == 11) {
			basisSizes = {3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2};
		} else if (d == 12) {
			basisSizes = {12, 6, 6, 6};
		} else if (d == 13) {
			basisSizes = {3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
		} else if (d == 14) {
			basisSizes = {14, 7, 7, 7};
		} else if (d == 15) {
			basisSizes = {10, 10, 10, 7};
		} else if (d == 16) {
			basisSizes = {16, 14, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		} else if (d == 17) {
			basisSizes = {17, 15, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		} else if (d == 18) {
			basisSizes = {18, 9, 9, 9};
		} else if (d == 19) {
			basisSizes = {19, 17, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		} else if (d == 20) {
			basisSizes = {20, 10, 10, 10};
		} else {
			basisSizes = std::vector<int>(n, d);
		}

	}

	// Useful quantities
	int numVarsNonConj = n*d*d;
	int numVars = 2*numVarsNonConj+1000;
	int conjDelta = numVarsNonConj;

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
					if (k < basisSizes[i]) {
						if (i == 0 && firstIsComputational) {
							if (i == m) {
								std::cout << 1;
							} else {
								std::cout << 0;
							}
						} else if ((m == 0 && firstElementIsOne) || (secondIsUniform && i == 1 && k == 0)) {
							std::cout << std::setprecision(3) << 1.0/std::sqrt(d);
						} else {
							std::cout << "x_" << nextInd << "+i*x_" << nextInd + conjDelta;
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
		for (int k=0; k<basisSizes[i]; k++) {

			// We assume the first vector of the second basis is uniform, so skip it
			if (secondIsUniform && i == 1 && k == 0) {
				continue;
			}

			// Second vector (not repeating the first)
			for (int j=i; j<n; j++) {
				for (int l=0; l<basisSizes[j]; l++) {

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
				for (int m=startInd; m<basisSizes[0]; m++) {
					int var1 = i*d*d + k*d + m;
					int var2 = var1 + conjDelta;
					Polynomial<double> extraEqn(numVars);
					extraEqn.addTerm(1, {var1,var1});
					extraEqn.addTerm(1, {var2,var2});
					extraEqn.addTerm(-1.0/d, {});
					eqns.push_back(extraEqn);
				}
			}

			// If we don't have implict normalization, need to add
			if (!firstIsComputational || basisSizes[0] < d) {
				Polynomial<std::complex<double>> extraEqn(numVars);
				int startInd = 0;
				if (firstElementIsOne) {
					extraEqn.addTerm(1.0/d);
					startInd = 1;
				}
				for (int m=startInd; m<d; m++) {
					int var1 = i*d*d + k*d + m;
					int var2 = var1 + conjDelta;
					extraEqn.addTerm(1.0, {var1, var1});
					extraEqn.addTerm(1.0, {var2, var2});
				}
				extraEqn.addTerm(-1.0);
				eqns.push_back(std::real<double>(extraEqn));
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
				if (i != 1) {
					extraEqn.addTerm(-1.0/d);
				}
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
		for (int j=0; j<basisSizes[i]-1; j++) {
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
		if (ind1 > basisSizes[i]-1) {
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
			for (int j=0; j<basisSizes[i]; j++) {
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
		for (int j=0; j<basisSizes[i]; j++) {
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

	// Convert these symmetries into constraints which break them
	std::vector<Polynomial<double>> orderingCons;
	if (task == 1) {
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

	}

	// Combine these equations into a single object
	PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, orderingCons);

	// Use as few indices as possible
	std::unordered_map<int,int> reducedMap = prob.getMinimalMap();
	if (verbosity >= 2) {
		std::cout << "---------------------" << std::endl;
		std::cout << "Original Problem: " << std::endl;
		std::cout << "---------------------" << std::endl;
		std::cout << prob << std::endl;
	}
	prob = prob.replaceWithVariable(reducedMap);
	if (verbosity >= 2) {
		std::cout << "---------------------" << std::endl;
		std::cout << "Reduced Problem: " << std::endl;
		std::cout << "---------------------" << std::endl;
		std::cout << prob << std::endl;
		std::cout << "---------------------" << std::endl;
		std::cout << "Summary:" << std::endl;
		std::cout << "---------------------" << std::endl;
	}
	int numVectors = 0;
	for (int i=0; i<n; i++) {
		numVectors += basisSizes[i];
	}
	if (verbosity >= 1) {
		std::cout << "set sizes: " << basisSizes << ", vars: " << prob.maxVariables << ", vectors: " << numVectors << ", equations: " << eqns.size() << std::endl;
	}

	// If told to find a feasible point
	if (task == 0) {

		// Heading if verbose
		if (verbosity >= 2) {
			std::cout << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << "Optimization:" << std::endl;
			std::cout << "---------------------" << std::endl;
		}

		// Find a feasible point of the equality constraints
		std::cout << std::scientific;
		std::vector<double> x;
		x = prob.findFeasibleEqualityPoint(-1, alpha, tolerance, maxIters, cores, verbosity, 1.0/std::sqrt(d), stabilityTerm);
		double maxVal = -1000;
		for (int i=0; i<prob.conZero.size(); i++) {
			maxVal = std::max(maxVal, std::abs(prob.conZero[i].eval(x)));
		}
		std::cout << "max viol = " << maxVal << std::endl;

		// If verbose, also display the overlaps
		if (verbosity >= 2) {

			std::vector<std::vector<std::vector<std::complex<double>>>> bases(n, std::vector<std::vector<std::complex<double>>>(d, std::vector<std::complex<double>>(d, 0.0)));
			std::cout << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << "Result:" << std::endl;
			std::cout << "---------------------" << std::endl;
			int nextInd = 0;
			for (int i=0; i<n; i++) {
				std::cout << std::endl;
				for (int k=0; k<d; k++) {
					std::cout << "{";
					for (int m=0; m<d; m++) {
						if (k < basisSizes[i]) {
							if (i == 0 && firstIsComputational) {
								bases[i][k][m] = int(m == k);
							} else if ((m == 0 && firstElementIsOne) || (secondIsUniform && i == 1 && k == 0)) {
								bases[i][k][m] = 1.0/std::sqrt(d);
							} else {
								bases[i][k][m] = x[reducedMap[nextInd]] + 1i*x[reducedMap[nextInd+conjDelta]];
							}
							std::cout << bases[i][k][m];
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

			std::cout << "---------------------" << std::endl;
			std::cout << "Abs of Gram matrix:" << std::endl;
			std::cout << "---------------------" << std::endl;
			std::vector<std::vector<double>> gramMat(n*d, std::vector<double>(n*d));
			for (int i1=0; i1<n; i1++) {
				for (int k1=0; k1<d; k1++) {
					for (int i2=0; i2<n; i2++) {
						for (int k2=0; k2<d; k2++) {
							std::complex<double> sum = 0;
							for (int m=0; m<d; m++) {
								sum += std::conj(bases[i1][k1][m])*bases[i2][k2][m];
							}
							gramMat[i1*d+k1][i2*d+k2] = std::abs(sum);
						}
					}
				}
			}
			std::cout << std::setprecision(4) << std::endl;
			std::cout << gramMat << std::endl;

		}

	// If told to identify redundant constraints
	} else if (task == 2) {

		// See if there are any constraints that can be removed without LoG
		std::vector<int> redundantInds;
		for (int i=0; i<prob.conZero.size(); i++) {

			// Get the con to make
			Polynomial<double> toMatch = prob.conZero[i];

			// Erase various constraints
			PolynomialProblem<double> probCopy = prob;
			probCopy.conZero.erase(probCopy.conZero.begin()+i);

			// Repeat a bunch to make sure it's not a fluke
			int zeroCount = 0;
			for (int j=0; j<50; j++) {
				auto res = probCopy.findFeasibleEqualityPoint(-1, alpha, tolerance, maxIters, cores, 0, 1.0/std::sqrt(d), stabilityTerm);
				auto error = std::abs(toMatch.eval(res));
				if (error < 1e-4) {
					zeroCount += 2;
				}
			}

			// If it's highly redundant, add it to the list
			std::cout << "con " << i << " is " << zeroCount << "% redundant" << std::endl;
			if (zeroCount > 90) {
				redundantInds.push_back(i);
			}

		}

		std::cout << "redundant list = " << redundantInds << std::endl;

		// Starting with small sets, growing to larger
		for (int num=2; num<redundantInds.size(); num++) {

			// Try a bunch of random selections of the redundants
			for (int l=0; l<20; l++) {

				// Remove this many of the redundants
				PolynomialProblem<double> probCopy = prob;
				std::vector<int> redundantCopy = redundantInds;
				std::vector<int> removedInds;
				for (int k=0; k<num; k++) {
					int r = int(redundantCopy.size()*(double(rand())/(RAND_MAX)));
					removedInds.push_back(redundantCopy[r]);
					redundantCopy.erase(redundantCopy.begin()+r);
				}

				// Only copy the ones we want to keep
				probCopy.conZero = {};
				std::vector<Polynomial<double>> removedCons;
				for (int k=0; k<prob.conZero.size(); k++) {
					if (std::find(removedInds.begin(), removedInds.end(), k) == removedInds.end()) {
						probCopy.conZero.push_back(prob.conZero[k]);
					} else {
						removedCons.push_back(prob.conZero[k]);
					}
				}

				// Repeat a bunch to make sure it's not a fluke
				int zeroCount = 0;
				for (int j=0; j<50; j++) {
					auto res = probCopy.findFeasibleEqualityPoint(-1, alpha, tolerance, maxIters, cores, 0, 1.0/std::sqrt(d), stabilityTerm);
					auto error = 0;
					for (int k=0; k<removedCons.size(); k++) {
						error += std::abs(removedCons[k].eval(res));
					}
					if (error < 1e-6) {
						zeroCount += 2;
					}
				}

				std::cout << num << " " << l << " " << zeroCount << " " << removedInds << std::endl;

			}

		}

	// If told to turn it into a binary problem TODO
	} else if (task == 3 || task == 4) {

		// Binarize
		prob.toBinaryProblem(numBits, 1.0/std::sqrt(d));

		if (verbosity >= 2) {
			std::cout << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << "Binarized Problem: " << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << prob << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << "Summary:" << std::endl;
			std::cout << "---------------------" << std::endl;
		}
		if (verbosity >= 1) {
			std::cout << "binary vars: " << prob.maxVariables << ", equations: " << eqns.size() << std::endl;
		}
		if (verbosity >= 2) {
			std::cout << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << "Optimisation:" << std::endl;
			std::cout << "---------------------" << std::endl;
		}

		// Brute force if it's small enough
		if (task == 4) {
			auto res = prob.bruteForce();
			std::cout << res.first << std::endl;
			std::cout << res.second << std::endl;
		}

	// If told to prove the search space is infeasible
	} else if (task == 1) {

		// Heading if verbose
		if (verbosity >= 2) {
			std::cout << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << "Optimization:" << std::endl;
			std::cout << "---------------------" << std::endl;
		}

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

	return 0;

}
	
