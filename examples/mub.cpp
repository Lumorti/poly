#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Get the problem from the args
	std::srand(time(0));
	int d = 2;
	int n = 4;
	std::string task = "infeasible";
	bool useFull = false;
	bool useQuadratic = true;
	bool removeLinear = true;
	int maxIters = -1;
	int verbosity = 1;
	double testParam = 43;
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
			std::cout << " -v [int]    set the verbosity level (0,1,2)" << std::endl;
			std::cout << " -2          use quadratic equations" << std::endl;
			std::cout << " -4          use quartic equations" << std::endl;
			std::cout << " -w          use whole bases, not partial" << std::endl;
			std::cout << " -f          try to find a feasible point" << std::endl;
			std::cout << " -r          don't use linear reductions" << std::endl;
			std::cout << " -o [str]    log points to a csv file" << std::endl;
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
		} else if (arg == "-t" && i+1 < argc) {
			testParam = std::stod(argv[i+1]);
			i++;
		} else if (arg == "-l" && i+1 < argc) {
			level = argv[i+1];
			i++;
		} else if (arg == "-v" && i+1 < argc) {
			verbosity = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-o" && i+1 < argc) {
			fileName = argv[i+1];
			i++;
		} else if (arg == "-4") {
			useQuadratic = false;
		} else if (arg == "-2") {
			useQuadratic = true;
		} else if (arg == "-w") {
			useFull = true;
		} else if (arg == "-f") {
			task = "feasible";
		} else if (arg == "-i") {
			task = "infeasible";
		} else if (arg == "-r") {
			removeLinear = false;
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
			dLimits.push_back({2, 1, 1, 1, 1, 1, 1, 1});
		} else if (d == 3) {
			dLimits.push_back({3, 1, 1, 1, 1, 1, 1, 1});
		} else if (d == 4) {
			//dLimits.push_back({4, 4, 3, 3, 3, 3});
			//dLimits.push_back({4, 4, 4, 3, 3, 3});
			//dLimits.push_back({4, 4, 4, 4, 3, 3});
			//dLimits.push_back({4, 4, 4, 4, 4, 3});
			dLimits.push_back({4, 2, 1, 1, 1, 1});
		} else if (d == 5) {
			dLimits.push_back({5, 2, 2, 1, 1, 1, 1});
		} else if (d == 6) {
			//dLimits.push_back({6, 2, 2, 2, 2, 1, 1, 1});
			dLimits.push_back({6, 5, 3, 1}); 
			//dLimits.push_back({6, 3, 3, 3}); 
		} else if (d == 7) {
			dLimits.push_back({7, 2, 2, 2, 2, 2, 1, 1, 1});
		} else if (d == 8) {
			dLimits.push_back({8, 2, 2, 2, 2, 2, 2, 1, 1});
		} else if (d == 9) {
			dLimits.push_back({9, 2, 2, 2, 2, 2, 2, 2, 1});
		} else if (d == 10) {
			dLimits.push_back({10, 7, 3, 1});
			dLimits.push_back({10, 8, 3, 1});
			dLimits.push_back({10, 9, 3, 1});
		} else {
			dLimits.push_back(std::vector<int>(n, d));
		}
	}

	// For each different restriction
	for (int i2=0; i2<dLimits.size(); i2++) {

		// List the bases and the variables indices
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
						if (m == 0 || i == 0 || (i == 1 && k == 0)) {
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

		// The list of equations to fill
		std::vector<Polynomial<double>> eqns;

		// Generate equations
		int newVarInd = 2*numVarsNonConj;
		for (int i=1; i<n; i++) {
			for (int j=i; j<n; j++) {
				for (int k=0; k<dLimits[i2][i]; k++) {
					for (int l=0; l<dLimits[i2][j]; l++) {

						// Prevent repeat equations
						if (i == j && l <= k) {
							continue;
						}

						// (a+ib)*(c-id)
						Polynomial<std::complex<double>> eqn(numVars);
						for (int m=0; m<d; m++) {
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

							// If we want the final equations to be quadratic 
							if (useQuadratic) {

								// Constrain that this should equal a new complex number
								eqn.addTerm(-1, {newVarInd});
								eqn.addTerm(-1i, {newVarInd+1});

								// Constrain that this new complex should have mag 1/d
								Polynomial<double> extraEqn(numVars);
								extraEqn.addTerm(1, {newVarInd,newVarInd});
								extraEqn.addTerm(1, {newVarInd+1,newVarInd+1});
								extraEqn.addTerm(-1.0/d, {});
								eqns.push_back(extraEqn);
								newVarInd += 2;

							} else {

								// Get the magnitude of this
								eqn = std::conj(eqn)*eqn;

								// Which should be equal to 1/d
								eqn.addTerm(-1.0/d, {});

							}

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
			}
		}

		// All should have mag 1/sqrt(d) (since we're setting first basis to the comp)
		int ogEqns = eqns.size();
		for (int i=d*d; i<numVarsNonConj; i++) {
			for (int j=0; j<ogEqns; j++) {
				if (eqns[j].contains(i)) {
					Polynomial<double> extraEqn(numVars);
					extraEqn.addTerm(1, {i,i});
					extraEqn.addTerm(1, {i+conjDelta,i+conjDelta});
					extraEqn.addTerm(-1.0/d, {});
					eqns.push_back(extraEqn);
					break;
				}
			}
		}

		// Can assume the first basis is the computational
		std::vector<int> indsToReplace;
		std::vector<double> valsToReplace;
		for (int i=0; i<d; i++) {
			for (int j=0; j<d; j++) {
				indsToReplace.push_back(i*d+j);
				indsToReplace.push_back(i*d+j+conjDelta);
				if (i == j) {
					valsToReplace.push_back(1);
				} else {
					valsToReplace.push_back(0);
				}
				valsToReplace.push_back(0);
			}
		}

		// Can assume first vector of second basis is uniform
		for (int j=0; j<d; j++) {
			indsToReplace.push_back(d*d+j);
			indsToReplace.push_back(d*d+j+conjDelta);
			valsToReplace.push_back(1/std::sqrt(d));
			valsToReplace.push_back(0);
		}

		// Can assume first value of each vector is 1 / sqrt(d)
		for (int i=0; i<n; i++) {
			for (int j=0; j<d; j++) {
				indsToReplace.push_back(i*d*d+j*d);
				indsToReplace.push_back(i*d*d+j*d+conjDelta);
				valsToReplace.push_back(1/std::sqrt(d));
				valsToReplace.push_back(0);
			}
		}

		// Perform the replacement
		for (int i=0; i<eqns.size(); i++) {
			eqns[i] = eqns[i].replaceWithValue(indsToReplace, valsToReplace);
		}

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

		// Convert these symmetries into constraints with break them TODO
		std::vector<Polynomial<double>> orderingCons;
		for (int i=0; i<syms.size(); i++) {
			Polynomial<double> newCon(numVars);
			int pow = 0;
			for (auto const &pair: syms[i]) {
				newCon.addTerm(std::pow(2, pow), {pair.first});
				newCon.addTerm(-std::pow(2, pow), {pair.second});
				pow++;
			}
			orderingCons.push_back(newCon);
		}

		// Combine these equations into a single object
		PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, orderingCons);

		// Remove variables using linear equalities if possible
		if (removeLinear) {
			prob = prob.removeLinear();
		}
		
		// Try to simplify the equations a bit
		for (int i=0; i<prob.conZero.size(); i++) {

			// Only consider non-normalization equations
			if (prob.conZero[i].size() > 3) {

				// Loop over the monoms of this
				std::vector<std::string> mons = prob.conZero[i].getMonomials();
				int digitsPerInd = prob.conZero[i].digitsPerInd;
				for (int m=0; m<mons.size(); m++) {

					// Find a squared term
					if (mons[m].size() == 2*digitsPerInd && mons[m].substr(0,digitsPerInd) == mons[m].substr(digitsPerInd,digitsPerInd)) {
						int realInd = std::stoi(mons[m].substr(0,digitsPerInd));
						double realCoeff = prob.conZero[i][mons[m]];

						// Now find the imag of this 
						for (int m2=0; m2<mons.size(); m2++) {
							if (mons[m2].size() == 2*digitsPerInd && mons[m2].substr(0,digitsPerInd) == mons[m2].substr(digitsPerInd,digitsPerInd) && std::stoi(mons[m2].substr(0,digitsPerInd)) == realInd+conjDelta) {
								int imagInd = std::stoi(mons[m2].substr(0,digitsPerInd));
								double imagCoeff = prob.conZero[i][mons[m2]];

								// Add a scaled version of the corresponding normal equation to simplify
								if (std::abs(realCoeff-imagCoeff) < 1e-6) {
									Polynomial<double> adjustment(numVars);
									adjustment.addTerm(realCoeff/d, {});
									adjustment.addTerm(-realCoeff, {realInd, realInd});
									adjustment.addTerm(-realCoeff, {imagInd, imagInd});
									prob.conZero[i] += adjustment;
									break;
								}

							}
						}

					}

				}

			}

			// Get rid of any almost zeros
			prob.conZero[i] = prob.conZero[i].prune();

		}

		// Use as few indices as possible
		std::unordered_map<int,int> reducedMap = prob.getMinimalMap();
		prob = prob.replaceWithVariable(reducedMap);
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
		std::cout << dLimits[i2] << " " << prob.maxVariables << std::endl;

		// If told to find a feasible point
		if (task == "feasible") {

			// Find a upper bound
			std::cout << std::scientific;
			std::vector<double> x = prob.findFeasiblePoint(-1, 0.1, 1e-12, 1000000, 4, false, 1.0/std::sqrt(d));
			double maxVal = -1000;
			for (int i=0; i<prob.conZero.size(); i++) {
				maxVal = std::max(maxVal, std::abs(prob.conZero[i].eval(x)));
			}
			std::cout << "max viol = " << maxVal << std::endl;

			// Reverse the map
			std::vector<double> origX(numVars, 0);
			for (auto const &pair: reducedMap) {
				origX[pair.first] = x[pair.second];
			}
			double maxVal2 = -1000;
			for (int i=0; i<eqns.size(); i++) {
				maxVal2 = std::max(maxVal2, std::abs(eqns[i].eval(origX)));
			}
			std::cout << "max viol of orig = " << maxVal2 << std::endl;

			// Put the bases in a nicer form
			std::vector<std::vector<std::vector<std::complex<double>>>> bases(n, std::vector<std::vector<std::complex<double>>>(d, std::vector<std::complex<double>>(d, 0)));
			nextInd = 0;
			for (int i=0; i<n; i++) {
				for (int k=0; k<d; k++) {
					for (int m=0; m<d; m++) {
						if (i > 0 && m == 0) {
							bases[i][k][m] = 1/std::sqrt(d);
						} else if (i == 1 && k == 0) {
							bases[i][k][m] = 1/std::sqrt(d);
						} else if (i == 0 && m == k) {
							bases[i][k][m] = 1;
						} else {
							bases[i][k][m] = origX[nextInd] + 1i*origX[nextInd+conjDelta];
						}
						nextInd++;
					}
				}
			}

			// List the bases 
			std::cout << "In a nicer form:" << std::endl;
			for (int i=0; i<n; i++) {
				std::cout << bases[i] << std::endl;
			}

			// Calculate overlaps
			if (!removeLinear) {
				for (int i=0; i<n; i++) {
					for (int j=0; j<n; j++) {
						std::vector<std::vector<double>> overlap(d, std::vector<double>(d, -1));
						for (int k=0; k<dLimits[i2][i]; k++) {
							for (int l=0; l<dLimits[i2][j]; l++) {
								std::complex<double> beforeSquare = 0;
								for (int m=0; m<d; m++) {
									beforeSquare += std::conj(bases[i][k][m]) * bases[j][l][m];
								}
								overlap[k][l] = std::real(std::conj(beforeSquare) * beforeSquare);
							}
						}
						std::cout << "overlap between " << i << " and " << j << ":" << std::endl;
						std::cout << overlap << std::endl;
					}
				}
			}

		// If told to prove the search space is infeasible
		} else if (task == "infeasible") {

			// If we're using a higher level mat, add higher-order cons
			if (level.find("2") != std::string::npos) {
				int ogCons = prob.conZero.size();
				for (int i=0; i<ogCons; i++) {
					prob.conZero.push_back(prob.conZero[i]*prob.conZero[i]);
				}
			}

			// If we're using a higher level mat, add higher-order cons
			if (level.find("3") != std::string::npos) {
				int ogCons = prob.conZero.size();
				for (int i=0; i<ogCons; i++) {
					prob.conZero.push_back(prob.conZero[i]*prob.conZero[i]*prob.conZero[i]);
				}
			}

			// If we're using a higher level mat, add higher-order cons
			if (level.find("4") != std::string::npos) {
				int ogCons = prob.conZero.size();
				for (int i=0; i<ogCons; i++) {
					prob.conZero.push_back(prob.conZero[i]*prob.conZero[i]*prob.conZero[i]*prob.conZero[i]);
				}
			}

			// Try to prove infeasiblity
			std::vector<double> params = {1};
			if (d == 3) {
				params[0] = 43;
			}
			prob.proveInfeasible(maxIters, level, 1.0/std::sqrt(d), fileName, verbosity);
			//prob.proveInfeasibleRadial(maxIters, level, 1.0/std::sqrt(d), fileName, params, verbosity);
			
		}

	}

	return 0;

}
	
