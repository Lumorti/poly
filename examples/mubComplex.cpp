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
	std::string level = "1f";
	std::string fileName = "";
	for (int i=0; i<argc; i++) {
		std::string arg = argv[i];
		if (arg == "-h") {
			std::cout << " -d [int]    set the dimension" << std::endl;
			std::cout << " -n [int]    set the number of bases" << std::endl;
			std::cout << " -l [str]    set the level for the relaxation e.g. 1+2f,3p" << std::endl;
			std::cout << " -i [int]    set max iterations (-1 for no limit)" << std::endl;
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
		} else if (arg == "-l" && i+1 < argc) {
			level = argv[i+1];
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
							std::cout << nextInd;
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
		std::cout << "conjugate delta = " << conjDelta << std::endl;
		std::cout << std::endl;

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

						// z_i* z_j
						Polynomial<double> eqn(numVars);
						Polynomial<double> eqnConj(numVars);
						for (int m=0; m<d; m++) {
							int var1 = i*d*d + k*d + m;
							int var2 = j*d*d + l*d + m;
							eqn.addTerm(1, {var2+conjDelta, var1});
							eqnConj.addTerm(1, {var2, var1+conjDelta});
						}

						// For the MUB-ness equations
						if (i != j) {

							// If we want the final equations to be quadratic 
							if (useQuadratic) {

								// Constrain that this should equal a new complex number
								eqn.addTerm(-1, {newVarInd});
								eqnConj.addTerm(-1, {newVarInd+conjDelta});

								// Constrain that this new complex should have mag 1/d
								Polynomial<double> extraEqn(numVars);
								extraEqn.addTerm(1, {newVarInd,newVarInd+conjDelta});
								extraEqn.addTerm(-1.0/d, {});
								eqns.push_back(extraEqn);
								newVarInd += 2;

							} else {

								// Get the magnitude of this
								eqn = eqnConj*eqn;

								// Which should be equal to 1/d
								eqn.addTerm(-1.0/d, {});
								eqnConj = eqn;

							}

						}

						// Add this equation if it's not empty
						if (eqn.size() > 0) {
							eqns.push_back(eqn);
							eqns.push_back(eqnConj);
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
					extraEqn.addTerm(1, {i+conjDelta,i});
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
					valsToReplace.push_back(1);
				} else {
					valsToReplace.push_back(0);
					valsToReplace.push_back(0);
				}
			}
		}

		// Can assume first vector of second basis is uniform
		for (int j=0; j<d; j++) {
			indsToReplace.push_back(d*d+j);
			indsToReplace.push_back(d*d+j+conjDelta);
			valsToReplace.push_back(1/std::sqrt(d));
			valsToReplace.push_back(1/std::sqrt(d));
		}

		// Can assume first value of each vector is 1 / sqrt(d)
		for (int i=0; i<n; i++) {
			for (int j=0; j<d; j++) {
				indsToReplace.push_back(i*d*d+j*d);
				indsToReplace.push_back(i*d*d+j*d+conjDelta);
				valsToReplace.push_back(1/std::sqrt(d));
				valsToReplace.push_back(1/std::sqrt(d));
			}
		}

		// Perform the replacement
		for (int i=0; i<eqns.size(); i++) {
			eqns[i] = eqns[i].replaceWithValue(indsToReplace, valsToReplace);
		}

		// Find the symmetries of the problem
		std::vector<std::unordered_map<int,int>> syms;

		// Ordering within each basis
		//int base = 2;
		//for (int i=1; i<n; i++) {
			//for (int j=0; j<dLimits[i2][i]-1; j++) {
				//if (i == 1 && j == 0) {
					//continue;
				//}
				//std::unordered_map<int,int> newSym;
				//for (int k=1; k<d; k++) {
					//int indLeft = i*d*d+j*d+k; 
					//int indRight = i*d*d+(j+1)*d+k; 
					//newSym[indLeft] = indRight;
					//newSym[indLeft+conjDelta] = indRight+conjDelta;
				//}
				//syms.push_back(newSym);
			//}
		//}

		// Ordering of the bases
		//for (int i=1; i<n-1; i++) {
			//int ind1 = 0;
			//if (i == 1) {
				//ind1 = 1;
			//}
			//if (ind1 > dLimits[i2][i]-1) {
				//continue;
			//}
			//std::unordered_map<int,int> newSym;
			//for (int k=1; k<d; k++) {
				//int indLeft = i*d*d+ind1*d+k;
				//int indRight = (i+1)*d*d+0*d+k;
				//newSym[indLeft] = indRight;
				//newSym[indLeft+conjDelta] = indRight+conjDelta;
			//}
			//syms.push_back(newSym);
		//}

		// Ordering of the elements in the vectors
		//for (int k=1; k<d-1; k++) {
			//std::unordered_map<int,int> newSym;
			//for (int i=1; i<n; i++) {
				//for (int j=0; j<dLimits[i2][i]; j++) {
					//if (i == 1 && j == 0) {
						//continue;
					//}
					//int indLeft = i*d*d+j*d+k;
					//int indRight = i*d*d+j*d+k+1;
					//newSym[indLeft] = indRight;
					//newSym[indLeft+conjDelta] = indRight+conjDelta;
				//}
			//}
			//syms.push_back(newSym);
		//}

		// Convert these symmetries into constraints with break them
		std::vector<Polynomial<double>> orderingCons;
		//for (int i=0; i<syms.size(); i++) {
			//Polynomial<double> newCon(numVars);
			//int pow = 0;
			//for (auto const &pair: syms[i]) {
				//newCon.addTerm(std::pow(2, pow), {pair.first});
				//newCon.addTerm(-std::pow(2, pow), {pair.second});
				//pow++;
			//}
			//orderingCons.push_back(newCon);
		//}

		// Combine these equations into a single object
		PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, orderingCons);

		// Remove variables using linear equalities if possible
		if (removeLinear) {
			prob = prob.removeLinear();
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

		// Nullstellensatz TODO
		prob.useNull(1);

	}

	return 0;

}
	
