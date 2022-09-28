#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Generate a random objective
	int n = 100;
	std::srand(time(0));
	Polynomial<double> obj(n);
	for (int i=0; i<n; i++) {
		double r = 10.0*(float(rand()) / RAND_MAX) - 5.0;
		obj.addTerm(r, {i});
	}

	// Constraints
	std::vector<Polynomial<double>> consPositive;
	std::vector<Polynomial<double>> consZero;
	{
		Polynomial<double> newCon(n);
		newCon.addTerm(-0.5, {});
		newCon.addTerm(1, {0});
		newCon.addTerm(1, {1});
		newCon.addTerm(1, {2});
		consZero.push_back(newCon);
	}

	// Create the problem
	PolynomialBinaryProblem<double> prob(obj, consZero, consPositive);
	prob.lowerBoundNew(30);

}
