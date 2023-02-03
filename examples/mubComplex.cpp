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
	int level = 1;
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
			level = std::stoi(argv[i+1]);
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

		// Add positivity cons for conj * orig
		std::vector<Polynomial<double>> consPositive;
		//consPositive.push_back(Polynomial<double>(numVars, 1, {9, 25}));
		//consPositive.push_back(Polynomial<double>(numVars, 1, {13, 29}));
		//consPositive.push_back(Polynomial<double>(numVars, 1, {36, 52}));

		//eqns = {};
		//eqns.push_back(Polynomial<double>(2, "1*{00}-1*{}"));
		//eqns.push_back(Polynomial<double>(2, "1*{11}+1*{}"));
		//eqns.push_back(Polynomial<double>(2, "1*{01}+1*{}"));
		//PolynomialProblem<double> prob(Polynomial<double>(2), eqns, {});

		//std::cout << eqns << std::endl;

		// Combine these equations into a single object
		PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, consPositive);

		// Remove variables using linear equalities if possible
		if (removeLinear) {
			prob = prob.removeLinear();
		}
		
		// Use as few indices as possible
		std::unordered_map<int,int> reducedMap = prob.getMinimalMap();
		prob = prob.replaceWithVariable(reducedMap);
		int newNumVars = prob.obj.maxVariables;
		std::cout << "---------------------" << std::endl;
		std::cout << "Index Mapping: " << std::endl;
		std::cout << "---------------------" << std::endl;
		std::cout << std::endl;
		std::cout << reducedMap << std::endl;
		std::cout << std::endl;
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{0}+1*{1}+1*{4}"));
		//prob.conZero.push_back(Polynomial<double>(numVars, "1*{0}+1*{1}+1*{2}+1*{3}-1*{4}+1*{5}"));
		//prob.conZero.push_back(Polynomial<double>(numVars, "1*{0}+1*{1}+1*{2}+1*{3}+1*{4}-1*{5}"));
		//prob.conZero.push_back(Polynomial<double>(numVars, "-1*{0}-1*{1}-1*{2}-1*{3}+1*{4}-1*{5}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{0}+1*{1}+1*{2}+1*{3}-1*{4}+1*{5}"));
		for (int i=0; i<prob.obj.maxVariables; i++) {
			prob.conPositive.push_back(Polynomial<double>(newNumVars, 1, {i}));
		}
		//auto mons = prob.getMonomials();
		//for (int i=0; i<mons.size(); i++) {
			//if (mons[i].size() > 0) {
				//if (double(rand())/RAND_MAX > 0.5) {
					//if (double(rand())/RAND_MAX > 0.5) {
						//prob.conPositive.push_back(Polynomial<double>(newNumVars, 1, mons[i]));
					//} else {
						//prob.conPositive.push_back(Polynomial<double>(newNumVars, -1, mons[i]));
					//}
				//}
			//}
		//}
		//for (int i=0; i<prob.obj.maxVariables; i++) {
			//for (int j=i; j<prob.obj.maxVariables; j++) {
				//prob.conPositive.push_back(Polynomial<double>(newNumVars, 1, {i,j}));
			//}
		//}
		//for (int i=0; i<prob.obj.maxVariables; i++) {
			//for (int j=i; j<prob.obj.maxVariables; j++) {
				//for (int k=j; k<prob.obj.maxVariables; k++) {
					//prob.conPositive.push_back(Polynomial<double>(numVars, 1, {i,j,k}));
				//}
			//}
		//}
		prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 0 8}"));
		prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 1 9}"));
		prob.conPositive.push_back(Polynomial<double>(numVars, "1*{12}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 0 6}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 1 7}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 2 8}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 3 9}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 410}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 511}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1215}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1316}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1417}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 1}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 2}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 3}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 4}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 5}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 6}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 7}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 8}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 9}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{10}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{11}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{12}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{13}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{14}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{15}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{16}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{17}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 0 0}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 1 1}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 2 2}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 3 3}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 4 4}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 5 5}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 6 6}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 7 7}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 8 8}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 9 9}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1010}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1111}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1212}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1313}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1414}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1515}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1616}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1717}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1215}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1316}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{1417}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 0 6}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 1 7}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 2 8}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 3 9}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 410}"));
		//prob.conPositive.push_back(Polynomial<double>(numVars, "1*{ 511}"));
		std::cout << "---------------------" << std::endl;
		std::cout << "Final Problem: " << std::endl;
		std::cout << "---------------------" << std::endl;
		std::cout << prob << std::endl;
		std::cout << dLimits[i2] << " " << prob.maxVariables << std::endl;

		// Nullstellensatz TODO
		prob.useNull(level, 1.0/std::sqrt(d), maxIters);

	}

	return 0;

}
	
