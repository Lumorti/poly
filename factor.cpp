#include "poly.h"

// Some useful definitions
using namespace std::complex_literals;

// Pretty print a vector
template <typename type>
void pretty(std::vector<type> a) {
	std::cout << "{ ";
	for (int i=0; i<a.size(); i++) {
		std::cout << a[i];
		if (i < a.size()-1) {
			std::cout << ", ";
		}
	}
	std::cout << "}" << std::endl;
}

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// 947*751=711197
	// Set up the problem
	//int a = 711197;
	int a = 27;
	if (argc > 1) {
		a = std::stoi(argv[1]);
	}
	int bitsPerFactor = std::ceil(std::log2(a));
	//int bitsPerFactor = std::log2(a);
	int numVars = 2*bitsPerFactor;

	// The list of equations to fill
	std::vector<Polynomial> eqns;

	// Generate equations x^2-x=0 to guarentee integers
	std::cout << "Generating equations..." << std::endl;
	for (int i=0; i<numVars; i++) {
		Polynomial eqnPreSquare(numVars);
		eqnPreSquare.addTerm(1, {i,i});
		eqnPreSquare.addTerm(-1, {i});
		eqns.push_back(eqnPreSquare*eqnPreSquare);
	}

	// Generate equation such that a=(x_1 + 2*x_2 + 4*x_3)*(x_4 + 2*x_5 + 4*x_6)
	Polynomial factor1(numVars);
	Polynomial factor2(numVars);
	for (int i=0; i<bitsPerFactor; i++) {
		factor1.addTerm(std::pow(2,i), {i});
		factor2.addTerm(std::pow(2,i), {i+bitsPerFactor});
	}
	Polynomial mainObj = factor1*factor2;
	mainObj.addTerm(-a, {});
	eqns.push_back(mainObj*mainObj);

	// Combine these to create a single polynomial
	Polynomial poly(numVars);
	for (int i=0; i<eqns.size(); i++) {
		poly += eqns[i];
	}

	// Integrate and then find a local minimum
	//std::vector<double> x = poly.integrate(0).findLocalMinimum(0.01, 1e-3);
	std::vector<double> x = { 0, 0, 0, 1, 1, 1, 0, 0, 0, 0};
	std::cout << poly.eval(x) << std::endl;
	std::cout << poly.integrate(0).differentiate(0).eval(x) << std::endl;

	// Link back to the original problem
	int val1 = 0;
	int val2 = 0;
	for (int i=0; i<bitsPerFactor; i++) {
		val1 += std::pow(2, i)*std::round(x[i]);
		val2 += std::pow(2, i)*std::round(x[i+bitsPerFactor]);
	}
	std::cout << a << " = " << val1 << " * " << val2 << std::endl;

	return 0;

}
	
