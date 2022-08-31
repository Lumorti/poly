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

	// The list of equations to fill
	PolynomialSystem<double> eqns;

	// Generate equations
	std::cout << "Generating equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			for (int k=0; k<d; k++) {
				for (int l=0; l<d; l++) {

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

					// Both the real and imaginary parts should be 0
					eqns.addPolynomial(Polynomial<double>(eqn.real()));
					eqns.addPolynomial(Polynomial<double>(eqn.imag()));

				}
			}
		}
	}

	// Combine these to create a single polynomial
	std::cout << "Creating single polynomial..." << std::endl;
	Polynomial<double> poly(numVars);
	for (int i=0; i<eqns.size(); i++) {
		poly += eqns[i]*eqns[i];
	}

	// Find a root of this, the first value will go to zero but that's okay
	// (could also use an auxiliary variable, but here not needed)
	std::cout << "Attempting to find a root..." << std::endl;
	std::vector<double> x = poly.findRoot(0, 0.1, 1e-10, 10000000);
	std::cout << "Testing this x = " << poly.eval(x) << std::endl;

	return 0;

}
	
