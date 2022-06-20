#include "poly.h"

// Some useful definitions
using namespace std::complex_literals;

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
	int numVars = 2*n*d*d;
	int conjDelta = numVarsNonConj;
	double rt2 = 1.0/std::sqrt(2.0);

	// Known ideal
	std::vector<double> idealX(numVars);
	if (d == 2 && n == 2) {
		idealX = {1, 0, 0, 1, rt2, rt2, rt2, -rt2, 0, 0, 0, 0, 0,   0,   0,    0};
	} else if (d == 2 && n == 3) {
		idealX = {1, 0, 0, 1, rt2, rt2, rt2, -rt2,  rt2, 0,   rt2, 0,
									  0, 0, 0, 0, 0,   0,   0,    0,    0,   rt2, 0,   -rt2};
	}

	// The list of equations to fill
	std::vector<Polynomial> eqns;

	// Generate equations
	std::cout << "Generating equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			for (int k=0; k<d; k++) {
				for (int l=0; l<d; l++) {

					// Get the equation before squaring
					Polynomial eqnPreSquare(numVars);
					for (int m=0; m<d; m++) {
						int var1 = i*d*d + k*d + m;
						int var2 = j*d*d + l*d + m;
						int var3 = var1 + conjDelta;
						int var4 = var2 + conjDelta;
						eqnPreSquare.addTerm(1, {var1, var2});
						eqnPreSquare.addTerm(1, {var3, var4});
						eqnPreSquare.addTerm(1i, {var1, var4});
						eqnPreSquare.addTerm(-1i, {var2, var3});
					}

					// For the normalisation equations
					if (i == j && k == l) {
						eqnPreSquare.addTerm(-1.0, {});
					}

					// Square this equation
					Polynomial eqn = eqnPreSquare.conjugate() * eqnPreSquare;

					// For the MUB-ness equations
					if (i != j) {
						eqn.addTerm(-1.0/d, {});
						eqn = eqn.conjugate() * eqn;
					}

					// Add it to the list
					eqns.push_back(eqn);

				}
			}
		}
	}

	// Combine these to create a single polynomial
	std::cout << "Creating single polynomial..." << std::endl;
	Polynomial poly(numVars);
	for (int i=0; i<eqns.size(); i++) {
		poly += eqns[i];
	}

	// Integrate and then find a local minimum of this
	std::vector<double> x = poly.integrate(0).findLocalMinimum(0.1);

	return 0;

}
	
