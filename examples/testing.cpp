#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Some params
	int n = 5;
	if (argc > 1) {
		n = std::stoi(argv[1]);
	}
	std::srand(time(NULL));
	std::string seed = std::to_string(std::rand());
	if (argc > 2) {
		seed = argv[2];
	}
	std::cout << "seed = " << seed << std::endl;
	std::srand(std::hash<std::string>{}(seed));
	Polynomial<double> obj(n);

	// All to all random -1 to 1
	//double maxVal = 100;
	//for (int i=0; i<n; i++) {
		//for (int j=i+1; j<n; j++) {
			//double r = 2*maxVal*(float(rand()) / RAND_MAX) - maxVal;
			//obj.addTerm(r, {i,j});
		//}
	//}

	// Factorization of 143
	obj.addTerm(-12, {0});
	obj.addTerm(-50, {1});
	obj.addTerm(-25, {2});
	obj.addTerm(12, {3});
	obj.addTerm(-24, {4});
	obj.addTerm(34, {0,1});
	obj.addTerm(-4, {0,2});
	obj.addTerm(-4, {0,3});
	obj.addTerm(8, {0,4});
	obj.addTerm(-8, {1,2});
	obj.addTerm(-8, {1,3});
	obj.addTerm(-8, {1,4});
	obj.addTerm(17, {2,3});
	obj.addTerm(-24, {2,4});
	obj.addTerm(-24, {3,4});

	// Limited connectivity all -1
	// 20 1 has min at -43, upper bound doesn't match
	// {-1, -1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1}
	// seed = 1941826014
	//for (int i=0; i<n; i++) {
		//for (int j=0; j<5; j++) {
			//int k = n*(float(rand()) / RAND_MAX);
			//if (k != i) {
				//obj.addTerm(1, {i,k});
			//}
		//}
	//}

	// Create the problem
	PolynomialBinaryProblem<double> prob(obj, {}, {});
	//std::cout << obj << std::endl;
	if (argc > 3) {
		auto sol = prob.bruteForce();
		std::cout << "brute force = " << sol.first << " " << sol.second << std::endl;
	}
	prob.lowerBound2(10000);

}
