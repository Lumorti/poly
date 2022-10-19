#include "../poly.h"

bool spaceSize(std::vector<int> dims1, std::vector<int> dims2) {
	int space1 = 0;
	for (int i=0; i<dims1.size(); i++) {
		space1 += dims1[i];
	}
	int space2 = 0;
	for (int i=0; i<dims2.size(); i++) {
		space2 += dims2[i];
	}
	return space1 < space2;
}

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

	// TODO remember we can set the first basis to computational basis

	std::vector<std::vector<int>> dLimits = {
		//{1, 1, 1, 1},
		//{2, 1, 1, 1},
		//{2, 2, 1, 1},
		//{2, 2, 2, 1},
		//{2, 2, 2, 2},
		//{3, 1, 1, 1},
		//{3, 2, 1, 1},
		//{3, 2, 2, 1},
		//{3, 2, 2, 2},
		//{3, 3, 1, 1},
		//{3, 3, 2, 1},
		//{3, 3, 2, 2},
		//{3, 3, 3, 1},
		//{3, 3, 3, 2},
		//{3, 3, 3, 3},
		//{4, 1, 1, 1},
		//{4, 2, 1, 1},
		//{4, 2, 2, 1},
		//{4, 2, 2, 2},
		//{4, 3, 1, 1},
		//{4, 3, 2, 1},
		//{4, 3, 2, 2},
		//{4, 3, 3, 1},
		//{4, 3, 3, 2},
		//{4, 3, 3, 3},
		//{4, 4, 1, 1},
		//{4, 4, 2, 1},
		//{4, 4, 2, 2},
		//{4, 4, 3, 1},
		//{4, 4, 3, 2},
		{4, 4, 3, 3},
		//{4, 4, 4, 1},
		{4, 4, 4, 2},
		{4, 4, 4, 3},
		{4, 4, 4, 4},
		//{5, 1, 1, 1},
		//{5, 2, 1, 1},
		//{5, 2, 2, 1},
		//{5, 2, 2, 2},
		//{5, 3, 1, 1},
		//{5, 3, 2, 1},
		//{5, 3, 2, 2},
		//{5, 3, 3, 1},
		//{5, 3, 3, 2},
		{5, 3, 3, 3},
		//{5, 4, 1, 1},
		//{5, 4, 2, 1},
		//{5, 4, 2, 2},
		//{5, 4, 3, 1},
		{5, 4, 3, 2},
		{5, 4, 3, 3},
		{5, 4, 4, 1},
		{5, 4, 4, 2},
		{5, 4, 4, 3},
		{5, 4, 4, 4},
		//{5, 5, 1, 1},
		//{5, 5, 2, 1},
		//{5, 5, 2, 2},
		{5, 5, 3, 1},
		{5, 5, 3, 2},
		{5, 5, 3, 3},
		{5, 5, 4, 1},
		{5, 5, 4, 2},
		{5, 5, 4, 3},
		{5, 5, 4, 4},
		{5, 5, 5, 1},
		{5, 5, 5, 2},
		{5, 5, 5, 3},
		{5, 5, 5, 4},
		{5, 5, 5, 5},
		//{6, 1, 1, 1},
		//{6, 2, 1, 1},
		//{6, 2, 2, 1},
		//{6, 2, 2, 2},
		//{6, 3, 1, 1},
		//{6, 3, 2, 1},
		//{6, 3, 2, 2},
		//{6, 3, 3, 1},
		//{6, 3, 3, 2},
		{6, 3, 3, 3},
		//{6, 4, 1, 1},
		//{6, 4, 2, 1},
		//{6, 4, 2, 2},
		//{6, 4, 3, 1},
		{6, 4, 3, 2},
		{6, 4, 3, 3},
		{6, 4, 4, 1},
		{6, 4, 4, 2},
		{6, 4, 4, 3},
		{6, 4, 4, 4},
		//{6, 5, 1, 1},
		//{6, 5, 2, 1},
		{6, 5, 2, 2},
		{6, 5, 3, 1},
		{6, 5, 3, 2},
		{6, 5, 3, 3},
		{6, 5, 4, 1},
		{6, 5, 4, 2},
		{6, 5, 4, 3},
		{6, 5, 4, 4},
		{6, 5, 5, 1},
		{6, 5, 5, 2},
		{6, 5, 5, 3},
		{6, 5, 5, 4},
		{6, 5, 5, 5},
		//{6, 6, 1, 1},
		{6, 6, 2, 1},
		{6, 6, 2, 2},
		{6, 6, 3, 1},
		{6, 6, 3, 2},
		{6, 6, 3, 3},
		{6, 6, 4, 1},
		{6, 6, 4, 2},
		{6, 6, 4, 3},
		{6, 6, 4, 4},
		{6, 6, 5, 1},
		{6, 6, 5, 2},
		{6, 6, 5, 3},
		{6, 6, 5, 4},
		{6, 6, 5, 5},
		{6, 6, 6, 1}
	};
	std::sort(dLimits.begin(), dLimits.end(), spaceSize);

	for (int i2=0; i2<dLimits.size(); i2++) {

		// The list of equations to fill
		std::vector<Polynomial<double>> eqns;

		// Generate equations
		//std::cout << "Generating equations..." << std::endl;
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

						// Only allow one element in the final basis
						//if ((i == n-1 && k > 0) || (j == n-1 && l > 0)) {
							//continue;
						//}
						//std::cout << "between " << i << "-" << k << " and " << j << "-" << l << std::endl;

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
		//std::cout << "Creating single polynomial..." << std::endl;
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
		int numVarsMin = 0;
		for (int k=0; k<dLimits[i2].size(); k++) {
			numVarsMin += dLimits[i2][k];
		}
		numVarsMin *= 2*d;
		std::cout << dLimits[i2] << "    " << eqns.size() << "   " << numVarsMin << std::endl;
		
		// Find a lower bound TODO
		//std::cout << "Finding lower bound..." << std::endl;
		//PolynomialProblem<double> prob(Polynomial<double>(numVars), eqns, {});
		//std::cout << prob << std::endl;
		//prob.lowerBound();
		//return 0;

		// Find a root of this using an auxilary variable
		//std::cout << "Number of equations: " << eqns.size() << std::endl;
		//std::cout << "Attempting to find a root..." << std::endl;
		std::vector<double> x = poly.findRoot(2, 0.5, 1e-10, 5000, 16, false);
		//std::cout << "   Testing this x = " << poly.eval(x) << std::endl;

	}

	return 0;

}
	
