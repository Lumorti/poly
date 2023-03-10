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
	float alpha = 0.9;
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
			std::cout << " -v [int]    set the verbosity level (0,1,2)" << std::endl;
			std::cout << " -p [int]    set the number of vars to split initially" << std::endl;
			std::cout << " -c [int]    set the number of cores to use" << std::endl;
			std::cout << " -o [str]    log points to a csv file" << std::endl;
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
			//dLimits.push_back({3, 2, 1, 1, 1});
		} else if (d == 4) {
			dLimits.push_back({4, 2, 1, 1, 1, 1});
		} else if (d == 5) {
			dLimits.push_back({5, 2, 2, 1, 1, 1, 1});
		} else if (d == 6) {
			dLimits.push_back({6, 5, 3, 1}); 
		} else if (d == 7) {
			dLimits.push_back({7, 2, 2, 2, 2, 1, 1, 1, 1});
		} else if (d == 8) {
			dLimits.push_back({8, 2, 2, 2, 2, 2, 1, 1, 1, 1});
		} else if (d == 9) {
			dLimits.push_back({9, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1});
		} else if (d == 10) {
			dLimits.push_back({10, 9, 3, 1});
		} else if (d == 11) {
			dLimits.push_back({11, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1});
		} else if (d == 12) {
			dLimits.push_back({12, 11, 3, 1});
		} else if (d == 13) {
			dLimits.push_back({13, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1});
		} else if (d == 14) {
			dLimits.push_back({14, 13, 3, 1});
		} else if (d == 15) {
			dLimits.push_back({15, 14, 3, 1});
		} else {
			dLimits.push_back(std::vector<int>(n, d));
		}

		// If going past infeasibility, pad with zeros
		while (dLimits[0].size() < n) {
			dLimits[0].push_back(0);
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

		// Generate equations, here iterating over all the vectors
		int newVarInd = 2*numVarsNonConj;
		for (int i=1; i<n; i++) {
			for (int k=0; k<dLimits[i2][i]; k++) {

				// We assume the first vector of the second basis is uniform, so skip it
				if (i == 1 && k == 0) {
					continue;
				}

				// Second vector (not repeating the first)
				for (int j=i; j<n; j++) {
					for (int l=0; l<dLimits[i2][j]; l++) {

						// Prevent repeat equations
						if (i == j && l <= k) {
							continue;
						}

						// (a+ib)*(c-id)
						Polynomial<std::complex<double>> eqn(numVars);
						for (int m=1; m<d; m++) {
							int var1 = i*d*d + k*d + m;
							int var2 = j*d*d + l*d + m;
							int var3 = var1 + conjDelta;
							int var4 = var2 + conjDelta;
							eqn.addTerm(1, {var1, var2});
							eqn.addTerm(1, {var3, var4});
							eqn.addTerm(1i, {var1, var4});
							eqn.addTerm(-1i, {var2, var3});
						}

						// Assuming the first element of both is 1/sqrt(d)
						eqn.addTerm(1.0/d);

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
				for (int m=1; m<d; m++) {
					int var1 = i*d*d + k*d + m;
					int var2 = var1 + conjDelta;
					Polynomial<double> extraEqn(numVars);
					extraEqn.addTerm(1, {var1,var1});
					extraEqn.addTerm(1, {var2,var2});
					extraEqn.addTerm(-1.0/d, {});
					eqns.push_back(extraEqn);
				}

				// Vs the uniform vector (e.g. |1/sqrt(d) * sum of basis| = 1/sqrt(d))
				Polynomial<std::complex<double>> extraEqn(numVars);
				for (int m=1; m<d; m++) {
					int var1 = i*d*d + k*d + m;
					int var2 = var1 + conjDelta;
					extraEqn.addTerm(1/std::sqrt(d), {var1});
					extraEqn.addTerm(1i/std::sqrt(d), {var2});
				}
				extraEqn.addTerm(1.0/d);
				extraEqn = std::conj(extraEqn)*extraEqn;
				for (auto const &pair: extraEqn.coeffs) {
					if (pair.first.size() == 2*extraEqn.digitsPerInd && pair.first.substr(0, extraEqn.digitsPerInd) == pair.first.substr(extraEqn.digitsPerInd, extraEqn.digitsPerInd)) {
						extraEqn.coeffs[pair.first] = 0.0;
						extraEqn.addTerm(0.5/(d*d));
					}
				}
				extraEqn.addTerm(-1.0/d);
				eqns.push_back(std::real<double>(extraEqn));

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
		std::cout << "---------------------" << std::endl;
		std::cout << "Index Mapping: " << std::endl;
		std::cout << "---------------------" << std::endl;
		std::cout << std::endl;
		std::cout << reducedMap << std::endl;
		std::cout << std::endl;
		prob = prob.replaceWithVariable(reducedMap);
		std::cout << "---------------------" << std::endl;
		std::cout << "Final Problem: " << std::endl;
		std::cout << "---------------------" << std::endl;
		std::cout << prob << std::endl;
		std::cout << dLimits[i2] << " " << prob.maxVariables << std::endl;

		// If told to find a feasible point
		if (task == "feasible") {
			std::cout << std::scientific;
			std::vector<double> x = prob.findFeasibleEqualityPoint(-1, alpha, 1e-12, maxIters, cores, verbosity, 1.0/std::sqrt(d));
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

	}

	return 0;

}
	
