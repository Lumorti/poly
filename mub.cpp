#include "poly.h"
#include <time.h>

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Sometimes we want this, sometimes we don't
	std::srand(std::time(NULL));

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
	int numVars = 2*n*d*d;
	int conjDelta = numVarsNonConj;
	double rt2 = 1.0/std::sqrt(2.0);

	// Known ideal
	std::vector<double> idealX(numVars);
	if (d == 2 && n == 2) {
		idealX = {1, 0, 0, 1, rt2, rt2, rt2, -rt2, 0, 0, 0, 0, 0, 0, 0, 0};
	} else if (d == 2 && n == 3) {
		idealX = {1, 0, 0, 1, rt2, rt2, rt2, -rt2, rt2, 0, rt2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rt2, 0, -rt2};
	}

	// The list of equations to fill
	std::vector<Polynomial<double>> eqns;

	// Generate equations
	std::cout << "Generating equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			for (int k=0; k<d; k++) {
				for (int l=0; l<d; l++) {

					// Get the equation before squaring
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

					// Both the real and imag parts should be 0
					eqns.push_back(Polynomial<double>(eqn.real()));
					eqns.push_back(Polynomial<double>(eqn.imag()));

				}
			}
		}
	}

	// TODO
	if (d == 2) {

		std::vector<int> indsToReplace = {};
		std::vector<double> valsToReplace = {};

		// 1 for 2
		if (d == 2 && n == 2) {
			indsToReplace = {0,1,2,3,  8,9,10,11};
			valsToReplace = {1, 0, 0, 1,    0, 0, 0, 0};
			//indsToReplace = {0,1,2,3,4,5,6,7,  8,9,10,11,12,13,14,15};
			//valsToReplace = {1, 0, 0, 1, rt2, rt2, rt2, -rt2, 0, 0, 0, 0, 0, 0, 0, 0};

		// 2 for 3
		} else if (d == 2 && n == 3) {
			indsToReplace = {0,1,2,3,  4,5,6,7,     12,13,14,15,16,17,18,19};
			valsToReplace = {1,0,0,1,  rt2,rt2,rt2,-rt2,      0,0,0,0,    0,0,0,0};
			//indsToReplace = {0,1,2,3,  4,5,6,7, 8,9,10,11,    12,13,14,15,16,17,18,19,20,21,22,23};
			//valsToReplace = {1,0,0,1,  rt2,rt2,rt2,-rt2,  rt2,0,rt2,0,    0,0,0,0,    0,0,0,0,    0,rt2,0,-rt2};

		// 3 for 4
		} else if (d == 2 && n == 4) {
			indsToReplace = {0,1,2,3,  4,5,6,7,   8,9,10,11,        16,17,18,19,  20,21,22,23,  24,25,26,27};
			valsToReplace = {1,0,0,1,  rt2,rt2,rt2,-rt2,  rt2,0,rt2,0,    0,0,0,0,    0,0,0,0,    0,rt2,0,-rt2};

		}

		for (int i=0; i<eqns.size(); i++) {
			eqns[i] = eqns[i].substitute(indsToReplace, valsToReplace);
			std::cout << "0 = " << eqns[i] << std::endl;
		}

	}

	// Combine these to create a single polynomial
	std::cout << "Creating single polynomial..." << std::endl;
	Polynomial<double> poly(numVars+1);
	for (int i=0; i<eqns.size(); i++) {
		poly += eqns[i]*eqns[i];
	}

	std::vector<int> varList = poly.getVars();
	std::vector<int> mapTo = {};
	for (int i=0; i<varList.size(); i++) {
		mapTo.push_back(i);
	}
	poly = poly.changeVars(varList, mapTo); 

	std::cout << poly << std::endl;

	// Find a root of this
	std::cout << "Attempting to find a root..." << std::endl;
	std::vector<double> x = poly.findRoot(0, 0.1, 1e-8, 10000);
	std::cout << "Testing this x = " << poly.eval(x) << std::endl;

	// Get the complex relaxation of this TODO
	std::cout << "Attempting to find a relaxed root..." << std::endl;
	Polynomial<double> relaxed = poly.getComplexRelaxation();
	int origVars = poly.numVars;
	relaxed.numVars += 1;
	std::vector<double> x2 = relaxed.findRoot(relaxed.numVars-1);
	std::vector<std::complex<double>> xComp(origVars);
	for (int i=0; i<origVars; i++) {
		xComp[i] = std::complex<double>(x2[i], x2[i+origVars]);
		std::cout << xComp[i] << " " << xComp[i]*xComp[i] << std::endl;
	}
	std::cout << "Testing this x = " << relaxed.eval(x2) << std::endl;

	// Maybe a Groebner basis? TODO

	return 0;

}
	
