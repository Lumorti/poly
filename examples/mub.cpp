#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Get the problem from the args
	int d = 2;
	int n = 4;
	std::string param = "";
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}
	if (argc > 3) {
		param = argv[3];
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
	//dLimits.push_back(std::vector<int>(n, d));
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
		//dLimits.push_back({6, 6, 6, 6}); 
	} else if (d == 7) {
		dLimits.push_back({7, 2, 2, 2, 2, 2, 1, 1, 1});
	} else {
		dLimits.push_back(std::vector<int>(n, d));
	}

	// For each different restriction
	for (int i2=0; i2<dLimits.size(); i2++) {

		// List the bases and the variables indices
		std::cout << "Basis indices:" << std::endl;
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

		// Add an ordering to all bases
		std::vector<Polynomial<double>> orderingCons;

		// Ordering within each basis
		for (int i=1; i<n; i++) {
			for (int j=0; j<dLimits[i2][i]-1; j++) {
				if (i == 1 && j == 0) {
					continue;
				}
				Polynomial<double> leftSide(numVars);
				Polynomial<double> rightSide(numVars);
				int leftPow = 0;
				int rightPow = 0;
				for (int k=1; k<d; k++) {
					leftSide.addTerm(std::pow(2,leftPow), {i*d*d+j*d+k});
					leftSide.addTerm(std::pow(2,leftPow+1), {i*d*d+j*d+k+conjDelta});
					leftPow += 2;
					rightSide.addTerm(std::pow(2,rightPow), {i*d*d+(j+1)*d+k});
					rightSide.addTerm(std::pow(2,rightPow+1), {i*d*d+(j+1)*d+k+conjDelta});
					rightPow += 2;
				}
				orderingCons.push_back(leftSide - rightSide);
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
			Polynomial<double> leftSide(numVars);
			Polynomial<double> rightSide(numVars);
			int leftPow = 0;
			int rightPow = 0;
			for (int k=1; k<d; k++) {
				leftSide.addTerm(std::pow(2,leftPow), {i*d*d+ind1*d+k});
				leftSide.addTerm(std::pow(2,leftPow+1), {i*d*d+ind1*d+k+conjDelta});
				leftPow += 2;
				rightSide.addTerm(std::pow(2,rightPow), {(i+1)*d*d+0*d+k});
				rightSide.addTerm(std::pow(2,rightPow+1), {(i+1)*d*d+0*d+k+conjDelta});
				rightPow += 2;
			}
			orderingCons.push_back(leftSide - rightSide);
		}

		// Ordering of the elements in the vectors
		for (int k=1; k<d-1; k++) {
			Polynomial<double> leftSide(numVars);
			Polynomial<double> rightSide(numVars);
			int leftPow = 0;
			int rightPow = 0;
			for (int i=1; i<n; i++) {
				for (int j=0; j<dLimits[i2][i]; j++) {
					if (i == 1 && j == 0) {
						continue;
					}
					leftSide.addTerm(std::pow(2,leftPow), {i*d*d+j*d+k});
					leftSide.addTerm(std::pow(2,leftPow+1), {i*d*d+j*d+k+conjDelta});
					leftPow += 2;
					rightSide.addTerm(std::pow(2,rightPow), {i*d*d+j*d+k+1});
					rightSide.addTerm(std::pow(2,rightPow+1), {i*d*d+j*d+k+1+conjDelta});
					rightPow += 2;
				}
			}
			orderingCons.push_back(leftSide - rightSide);
		}

		// Combine these equations into a single object
		PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, orderingCons);
		prob = prob.removeLinear();
		std::unordered_map<int,int> reducedMap = prob.getMinimalMap();
		prob = prob.replaceWithVariable(reducedMap);
		std::cout << std::endl;
		std::cout << prob << std::endl;
		std::cout << dLimits[i2] << " " << prob.maxVariables << std::endl;

		// If told to find a feasible point
		if (param == "f") {

			// Find a upper bound
			std::vector<double> x = prob.findFeasiblePoint(-1, 0.1, 1e-15, 10000, 4, false, 1.0/std::sqrt(d));
			double maxVal = -1000;
			for (int i=0; i<prob.conZero.size(); i++) {
				maxVal = std::max(maxVal, std::abs(prob.conZero[i].eval(x)));
			}
			std::cout << "max viol = " << maxVal << std::endl;
			//std::cout << "resulting x = " << x << std::endl;

			//// Reverse the map
			std::vector<double> origX(numVars, 0);
			for (auto const &pair: reducedMap) {
				origX[pair.first] = x[pair.second];
			}
			double maxVal2 = -1000;
			for (int i=0; i<eqns.size(); i++) {
				maxVal2 = std::max(maxVal2, std::abs(eqns[i].eval(origX)));
			}
			std::cout << "max viol of orig = " << maxVal2 << std::endl;
			//std::cout << "orig x = " << origX << std::endl;

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
			std::cout << std::fixed << std::setprecision(15);
			for (int i=0; i<n; i++) {
				std::cout << bases[i] << std::endl;
			}

			// Calculate overlaps
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

		} else {

			// Find a lower bound
			std::cout << std::scientific;
			prob.proveInfeasible(100000000);

		}

	}

	return 0;

}
	
