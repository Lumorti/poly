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
	dLimits.push_back(std::vector<int>(n, d));
	//if (d == 2) {
		//dLimits.push_back({2, 1, 1, 1, 1, 1, 1, 1});
	//} else if (d == 3) {
		//dLimits.push_back({3, 2, 1, 1, 1, 1, 1, 1});
	//} else if (d == 4) {
		////dLimits.push_back({4, 4, 3, 3, 3, 3});
		////dLimits.push_back({4, 4, 4, 3, 3, 3});
		////dLimits.push_back({4, 4, 4, 4, 3, 3});
		////dLimits.push_back({4, 4, 4, 4, 4, 3});
		//dLimits.push_back({4, 4, 4, 4, 4, 4});
	//} else if (d == 5) {
		//dLimits.push_back({5, 2, 2, 2, 1, 1, 1});
	//} else if (d == 6) {
		////dLimits.push_back({6, 2, 2, 2, 2, 1, 1, 1});
		////dLimits.push_back({6, 4, 4, 4}); 
		//dLimits.push_back({6, 6, 6, 6}); 
	//} else if (d == 7) {
		//dLimits.push_back({7, 2, 2, 2, 2, 2, 1, 1, 1});
	//} else {
		//dLimits.push_back(std::vector<int>(n, d));
	//}

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
						//if (i == j && (l < k || k == 1)) {
						if (i == j) {
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
		//std::cout << std::endl << "Equations before simplification: " << std::endl;
		for (int i=0; i<eqns.size(); i++) {
			eqns[i] = eqns[i].replaceWithValue(indsToReplace, valsToReplace);
			//std::cout << eqns[i] << std::endl;
		}

		// The symmetries of the problem
		//std::vector<std::unordered_map<int,int>> syms;

		//// For any order of bases
		//// (apart from the first, which we set to the computational)
		//std::vector<std::vector<int>> basisOrders;
		//std::vector<int> v1(n);
		//for (int i=0; i<v1.size(); i++) {
			//v1[i] = i;
		//}
		//do {
			//basisOrders.push_back(v1);
			//std::cout << "possible basis order: " << v1 << std::endl;
		//} while (std::next_permutation(v1.begin()+1, v1.end()));

		//// For any order of elements within a vector
		//// (apart from the first, which we set to be one anyway)
		//std::vector<std::vector<int>> intravectorOrder;
		//std::vector<int> v3(d);
		//for (int i=0; i<v3.size(); i++) {
			//v3[i] = i;
		//}
		//do {
			//intravectorOrder.push_back(v3);
			//std::cout << "possible intravector order: " << v3 << std::endl;
		//} while (std::next_permutation(v3.begin()+1, v3.end()));

		//// For any order of vectors within each basis
		//std::vector<std::vector<int>> vectorOrder;
		//std::vector<int> v2(d);
		//for (int i=0; i<v2.size(); i++) {
			//v2[i] = i;
		//}
		//do {
			//vectorOrder.push_back(v2);
			//std::cout << "possible vector order: " << v2 << std::endl;
		//} while (std::next_permutation(v2.begin(), v2.end()));

		//// Loop over the basis orders
		//for (int i=0; i<basisOrders.size(); i++) {

			//// Loop over the orders within all vectors
			//for (int j=0; j<intravectorOrder.size(); j++) {

				//// The different orders for each vector set
				//for (int k=0; k<std::pow(vectorOrder.size(), n-1); k++) {

					//// Count in a different base to loop through
					//std::vector<int> decomp(n-1);
					//int toConvert = k;
					//int base = vectorOrder.size();
					//int ind = 0;
					//while (toConvert) {
						//decomp[ind] = toConvert % base;
						//toConvert /= base;
						//ind++;
					//}
					//decomp.insert(decomp.begin(), 0);

					//// Construct the index ordering (symmetry map)
					//std::unordered_map<int,int> sym;

					//// For each basis, vector, then within each vector
					//int nextInd = 0;
					//for (int l1=0; l1<basisOrders[i].size(); l1++) {
						//for (int l2=0; l2<vectorOrder[decomp[l1]].size(); l2++) {
							//for (int l3=0; l3<intravectorOrder[j].size(); l3++) {

								//// Map the original ind to the new ind
								//sym[nextInd] = d*d*basisOrders[i][l1] + d*vectorOrder[decomp[l1]][l2] + intravectorOrder[j][l3]; 
								//sym[nextInd+conjDelta] = sym[nextInd]+conjDelta; 
								//nextInd++;

							//}
						//}
					//}

					//// TODO
					////std::cout << i << "  " << j << "  " << decomp << std::endl;
					////std::cout << sym << std::endl;
					////std::cout << std::endl;

					//// Add this symmetry
					//syms.push_back(sym);

				//}

			//}

		//}

		//std::cout << syms << std::endl;
		//std::cout << syms.size() << std::endl;
		//return 0;
		
		// Add an ordering to all bases TODO
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
				//std::cout << i << "-" << j << ">" << i << "-" << j+1 << std::endl;
				for (int k=1; k<d; k++) {
					leftSide.addTerm(std::pow(2,leftPow), {i*d*d+j*d+k});
					leftSide.addTerm(std::pow(2,leftPow+1), {i*d*d+j*d+k+conjDelta});
					leftPow += 2;
					rightSide.addTerm(std::pow(2,rightPow), {i*d*d+(j+1)*d+k});
					rightSide.addTerm(std::pow(2,rightPow+1), {i*d*d+(j+1)*d+k+conjDelta});
					rightPow += 2;
				}
				//std::cout << leftSide << " > " << rightSide << std::endl;
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
			//std::cout << i << "-" << ind1 << " > " << i+1 << "-" << 0 << std::endl;
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
			//std::cout << leftSide << " > " << rightSide << std::endl;
			orderingCons.push_back(leftSide - rightSide);
		}

		// Ordering of the elements in the vectors
		for (int k=1; k<d-1; k++) {
			//std::cout << "vec el" << k << " > vec el" << k+1 << std::endl;
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
			//std::cout << leftSide << " > " << rightSide << std::endl;
			orderingCons.push_back(leftSide - rightSide);
		}

		//std::cout << "With sym cons: " << std::endl << std::endl;
		//for (auto con : orderingCons) {
			//std::cout << con << std::endl << std::endl;
		//}

		//return 0;

		// Combine these equations into a single object
		PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, orderingCons);
		//prob = prob.removeLinear();
		std::unordered_map<int,int> reducedMap = prob.getMinimalMap();
		prob = prob.replaceWithVariable(reducedMap);
		std::cout << std::endl;
		std::cout << prob << std::endl;
		std::cout << dLimits[i2] << " " << prob.maxVariables << std::endl;

		// Find a upper bound
		std::vector<double> x = prob.findFeasiblePoint(-1, 0.1, 1e-15, 100000, 4, false, 1.0/std::sqrt(d));
		double maxVal = -1000;
		for (int i=0; i<prob.conZero.size(); i++) {
			maxVal = std::max(maxVal, std::abs(prob.conZero[i].eval(x)));
		}
		std::cout << "max viol = " << maxVal << std::endl;
		std::cout << "resulting x = " << x << std::endl;

		// Reverse the map
		std::vector<double> origX(numVars, 0);
		for (auto const &pair: reducedMap) {
			origX[pair.first] = x[pair.second];
		}
		//std::cout << "orig x = " << origX << std::endl;
		double maxVal2 = -1000;
		for (int i=0; i<eqns.size(); i++) {
			maxVal2 = std::max(maxVal2, std::abs(eqns[i].eval(origX)));
		}
		std::cout << "max viol of orig = " << maxVal2 << std::endl;

		// List the bases and the variables indices
		std::cout << "In a nicer form:" << std::endl;
		std::cout << std::fixed << std::setprecision(15);
		nextInd = 0;
		for (int i=0; i<n; i++) {
			std::cout << std::endl;
			for (int k=0; k<d; k++) {
				std::cout << "{";
				for (int m=0; m<d; m++) {
					if (k < dLimits[i2][i]) {
						if (m == 0 || i == 0 || (i == 1 && k == 0)) {
							std::cout << "-" ;
						} else {
							std::cout << origX[nextInd] << "+i(" << origX[nextInd+conjDelta] << ")";
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

		// Find a lower bound
		//prob.proveInfeasible(100000000);

	}

	return 0;

}
	
