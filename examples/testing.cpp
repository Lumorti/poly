#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Some params
	int n = 7;
	if (argc > 1) {
		n = std::stoi(argv[1]);
	}
	if (argc > 2) {
		std::srand(std::hash<std::string>{}(std::string(argv[2])));
	} else {
		std::srand(time(0));
	}
	Polynomial<double> obj(n);

	// All to all random -1 to 1
	//for (int i=0; i<n; i++) {
		//for (int j=i+1; j<n; j++) {
			//double r = 2.0*(float(rand()) / RAND_MAX) - 1.0;
			//obj.addTerm(r, {i,j});
		//}
	//}

	// Limited connectivity all -1
	// 20 1 has min at -43, upper bound doesn't match
	// {-1, -1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1}
	for (int i=0; i<n; i++) {
		for (int j=0; j<5; j++) {
			int k = n*(float(rand()) / RAND_MAX);
			if (k != i) {
				obj.addTerm(1, {i,k});
			}
		}
	}

	// Create the problem
	PolynomialBinaryProblem<double> prob(obj, {}, {});
	//std::cout << obj << std::endl;
	if (argc > 3) {
		auto sol = prob.bruteForce();
		std::cout << "brute force = " << sol.first << " " << sol.second << std::endl;
	}
	prob.lowerBound2(10000);

}
