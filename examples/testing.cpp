#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Generate a random Ising
	int n = 5;
	//std::srand(time(0));
	std::srand(7);
	Polynomial<double> obj(n);
	for (int i=0; i<n; i++) {
		for (int j=i+1; j<n; j++) {
			double r = 2.0*(float(rand()) / RAND_MAX) - 1.0;
			obj.addTerm(r, {i,j});
		}
	}
	std::cout << "obj = " << obj << std::endl;

	// Ising model system
	//int n = 5;
	//Polynomial<double> obj(n);
	//obj.addTerm(1, {0,1});
	//obj.addTerm(1, {0,2});
	//obj.addTerm(1, {1,2});
	//obj.addTerm(1, {2,3});

	// Constraints
	std::vector<Polynomial<double>> consPositive;
	std::vector<Polynomial<double>> consZero;
	//{
		//Polynomial<double> newCon(n);
		//newCon.addTerm(-0.5, {});
		//newCon.addTerm(1, {0});
		//newCon.addTerm(1, {1});
		//newCon.addTerm(1, {2});
		//consZero.push_back(newCon);
	//}

	// Create the problem
	PolynomialBinaryProblem<double> prob(obj, consZero, consPositive);
	auto sol = prob.bruteForce();
	std::cout << "brute force = " << sol.first << " " << sol.second << std::endl;
	prob.lowerBound2(10);

}
