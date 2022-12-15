#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Get the problem from the args
	int d = 2;
	int n = 4;
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}

	// Whether to finish with a quadratic or quartic set of equations
	bool useQuadratic = true;

	// Useful quantities
	int numVarsNonConj = n*d*d;
	int numVars = 2*numVarsNonConj+100000;
	int conjDelta = numVarsNonConj;
	double rt2 = 1.0/std::sqrt(2.0);

	// Different "bases"
	std::vector<std::vector<int>> dLimits;
	if (d == 2) {
		dLimits.push_back({2, 1, 1, 1});
	} else if (d == 3) {
		dLimits.push_back({3, 1, 1, 1, 1});
	} else if (d == 4) {
		dLimits.push_back({4, 2, 1, 1, 1, 1});
	} else if (d == 5) {
		dLimits.push_back({5, 2, 2, 1, 1, 1, 1});
	} else if (d == 6) {
		//dLimits.push_back({6, 2, 2, 2, 1, 1, 1, 1});
		dLimits.push_back({6, 5, 3, 1});
		//dLimits.push_back({6, 3, 3, 3});
	} else if (d == 7) {
		dLimits.push_back({7, 2, 2, 2, 2, 1, 1, 1, 1});
	} else {
		dLimits.push_back(std::vector<int>(n, d));
	}

	// For each different restriction
	for (int i2=0; i2<dLimits.size(); i2++) {

		// The list of equations to fill
		std::vector<Polynomial<double>> eqns;

		// Generate equations
		int newVarInd = 2*numVarsNonConj;
		for (int i=1; i<n; i++) {
			for (int j=i; j<n; j++) {
				for (int k=0; k<dLimits[i2][i]; k++) {
					for (int l=0; l<dLimits[i2][j]; l++) {

						// Prevent repeat equations
						if (i == j && l < k) {
							continue;
						}

						// (a+ib)*(c+id)
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

						// For the normalisation equations
						if (i == j && k == l) {
							eqn.addTerm(-1, {});
							continue;
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

		// All should have mag 1 (since we're setting first basis to the comp)
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
			std::cout << eqns[i] << std::endl;
		}

		// The symmetries of the problem TODO
		std::vector<std::unordered_map<int,int>> syms;

		// For any order of bases
		std::vector<std::vector<int>> basisOrders;
		std::vector<int> v1(n);
		for (int i=0; i<v1.size(); i++) {
			v1[i] = i;
		}
		do {
			basisOrders.push_back(v1);
			std::cout << "possible basis order: " << v1 << std::endl;
		} while (std::next_permutation(v1.begin()+1, v1.end()));

		// For any order of elements within a vector
		std::vector<std::vector<int>> intravectorOrder;
		std::vector<int> v3(d);
		for (int i=0; i<v3.size(); i++) {
			v3[i] = i;
		}
		do {
			intravectorOrder.push_back(v3);
			std::cout << "possible intravector order: " << v3 << std::endl;
		} while (std::next_permutation(v3.begin(), v3.end()));

		// For any order of vectors within each basis
		std::vector<std::vector<int>> vectorOrder;
		std::vector<int> v2(d);
		for (int i=0; i<v2.size(); i++) {
			v2[i] = i;
		}
		do {
			vectorOrder.push_back(v2);
			std::cout << "possible vector order: " << v2 << std::endl;
		} while (std::next_permutation(v2.begin(), v2.end()));

		// Loop over the basis orders
		for (int i=0; i<basisOrders.size(); i++) {
			std::cout << "basis order: " << i << std::endl;

			// Loop over the orders within all vectors
			for (int j=0; j<intravectorOrder.size(); j++) {
				std::cout << "   intravector order: " << j << std::endl;

				// The different orders for each vector set
				for (int k=0; k<std::pow(vectorOrder.size(), n-1); k++) {
					std::vector<int> decomp(n-1);
					int toConvert = k;
					int base = vectorOrder.size();
					int ind = 0;
					while (toConvert) {
						decomp[ind] = toConvert % base;
						toConvert /= base;
						ind++;
					}

					std::cout << "       vector order: " << decomp << std::endl;

				}

			}

		}

		std::cout << syms << std::endl;
		return 0;

		// Combine these equations into a single object
		PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, {}, syms);
		prob = prob.removeLinear();
		prob = prob.collapse();
		std::cout << prob << std::endl;
		std::cout << dLimits[i2] << " " << prob.maxVariables << std::endl;

		return 0;

		// Find a lower bound
		prob.proveInfeasible(1000);

		// Find a upper bound
		//prob.findFeasiblePoint(1, 0.5, 1e-13, 100000, 4, false);

	}

	return 0;

}
	
