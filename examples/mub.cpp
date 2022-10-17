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
	int numVars = 2*numVarsNonConj + 1;
	int conjDelta = numVarsNonConj;
	double rt2 = 1.0/std::sqrt(2.0);

	// The list of equations to fill
	std::vector<Polynomial<double>> eqns;

	// Generate equations
	std::cout << "Generating equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=i; j<n; j++) {
			for (int k=0; k<4; k++) {
				for (int l=0; l<4; l++) {

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

					// Only allow one element in the final basis
					//if ((i == n-1 && k > 0) || (j == n-1 && l > 0)) {
						//continue;
					//}
					std::cout << "between " << i << "-" << k << " and " << j << "-" << l << std::endl;

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

	// Combine these to create a single polynomial
	std::cout << "Creating single polynomial..." << std::endl;
	Polynomial<double> poly(numVars);
	for (int i=0; i<eqns.size(); i++) {
		poly += eqns[i]*eqns[i];
	}

	// Minimize number of equations needed TODO
	// for d2n4   10 (only one vector each basis)
	// for d3n5   57 (one less vector per and only one vec in last basis)
	// for d4n6   225 (one less vector per and only one vec in last basis)
	// for d5n7   433 (one less vector per and only one vec in last basis)
	// for d6n4   226 (one less vector per and only one vec in last basis)
	
	// Find a lower bound TODO
	//std::cout << "Finding lower bound..." << std::endl;
	//PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, {});
	//std::cout << prob << std::endl;
	//prob.lowerBound();
	//return 0;

	// Find a root of this using an auxilary variable
	std::cout << "Number of equations: " << eqns.size() << std::endl;
	std::cout << "Attempting to find a root..." << std::endl;
	std::vector<double> x = poly.findRoot(0, 0.5, 1e-10, 100000, 4, false);
	std::cout << "Testing this x = " << poly.eval(x) << std::endl;
	return 0;

}
	
