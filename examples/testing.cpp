#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// obj = x_1 - 2*x_2 
	// min is -3 at {-1,1}
	Polynomial<double> obj(2);
	obj.addTerm(1, {0});
	obj.addTerm(2, {1});

	// x_1 + x_2 >= 1
	//Polynomial<double> con(2);
	//obj.addTerm(1, {0});
	//obj.addTerm(1, {1});
	//obj.addTerm(-1, {});

	// Create the problem
	PolynomialBinaryProblem<double> prob(obj, {}, {});
	//PolynomialBinaryProblem<double> prob(obj, {}, {con});


}
