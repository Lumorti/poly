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
	PolynomialBinaryProblem<double> testProb(testObj, {}, {});

	// Brute force to find the optimum for comparison
	auto res = testProb.bruteForce();
	std::cout << "from brute forcing:" << std::endl;
	std::cout << res.first << " " << res.second << std::endl;
	std::cout << std::endl;

	// Try to find a lower bound for this valid
	auto resLower = testProb.lowerBound(10000, 2, false, false, 100);
	std::cout << std::endl;
	std::cout << "from lower bounding:" << std::endl;
	std::cout << resLower << std::endl;
	std::cout << std::endl;

	return 0;

}
