#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Get the problem from the args
	int d = 2;
	int n = 2;
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}

	// Useful quantities
	int numVarsNonConj = n*d*d;
	int numVars = 2*numVarsNonConj;
	int conjDelta = numVarsNonConj;
	double rt2 = 1.0/std::sqrt(2.0);

	// Different "bases"
	std::vector<std::vector<int>> dLimits;
	if (d == 6 && n == 4) {
		dLimits.push_back({6, 3, 3, 3});
	} else if (d == 5 && n == 4) {
		dLimits.push_back({5, 3, 3, 3});
	} else if (d == 7 && n == 4) {
		dLimits.push_back({7, 3, 3, 3});
	} else if (d == 10 && n == 4) {
		dLimits.push_back({10, 3, 3, 3});
	} else if (d == 2 && n == 4) {
		dLimits.push_back({1, 1, 1, 1});
	} else {
		dLimits.push_back(std::vector<int>(n, d));
	}

	// For each different restriction
	for (int i2=0; i2<dLimits.size(); i2++) {

		// The list of equations to fill
		std::vector<Polynomial<double>> eqns;

		// Generate equations
		for (int i=0; i<n; i++) {
			for (int j=i; j<n; j++) {
				for (int k=0; k<dLimits[i2][i]; k++) {
					for (int l=0; l<dLimits[i2][j]; l++) {

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
							eqn.addTerm(-1.0, {});
						}

						// For the MUB-ness equations
						if (i != j) {
							eqn = eqn.conjugate() * eqn;
							eqn.addTerm(-1.0/d, {});
						}

						// Split into real and imag parts (both should be 0)
						Polynomial<double> eqnReal = eqn.real();
						Polynomial<double> eqnImag = eqn.imag();

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

		// Can assume the first basis is the computational TODO
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
		for (int i=0; i<eqns.size(); i++) {
			eqns[i] = eqns[i].substitute(indsToReplace, valsToReplace);
		}

		// Collapse to the minimum number of vars
		std::vector<int> originalInds;
		for (int i=0; i<numVars; i++) {
			originalInds.push_back(i);
		}
		std::vector<int> newInds(numVars, -1);
		int nextInd = 0;
		for (int i=0; i<numVars; i++) {
			for (int j=0; j<eqns.size(); j++) {
				if (eqns[j].contains(i)) {
					newInds[i] = nextInd;
					nextInd++;
					break;
				}
			}
		}
		numVars = nextInd+1;
		for (int i=0; i<eqns.size(); i++) {
			eqns[i] = eqns[i].changeVariables(originalInds, newInds);
			eqns[i].maxVariables = numVars;
		}

		// Combine these to create a single polynomial
		Polynomial<double> poly(numVars);
		for (int i=0; i<eqns.size(); i++) {
			poly += eqns[i]*eqns[i];
		}

		// Minimize number of equations needed TODO
		std::cout << dLimits[i2] << "    " << eqns.size() << "   " << numVars << std::endl;
		
		// Find a lower bound TODO
		PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, {});
		prob.proveInfeasible();

		// Find a root of this using an auxilary variable
		//std::vector<double> x = poly.findRoot(0, 0.5, 1e-10, 10000, 16, false);
		//std::cout << "   Testing this x = " << poly.eval(x) << std::endl;

	}

	return 0;

}
	
