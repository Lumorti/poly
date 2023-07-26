#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Create a randomised Ising spin system
	int numSpins = 10;
	Polynomial<double> testObj(numSpins);
	for (int i=0; i<numSpins; i++) {
		for (int j=i+1; j<numSpins; j++) {
			testObj.addTerm(2*(double(rand())/RAND_MAX)-1, {i,j});
		}
	}

	// Create the binary problem
	PolynomialProblem<double> testProb(testObj, {}, {});
	testProb.isBinary = true;

	// Brute force to find the optimum for comparison
	auto res = testProb.bruteForce();
	std::cout << "from brute forcing:" << std::endl;
	std::cout << res.first << " " << res.second << std::endl;
	std::cout << std::endl;

	return 0;

}
