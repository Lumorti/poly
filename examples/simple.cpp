#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Create a polynomial with a maximum of two variables
	Polynomial<float> poly(2);

	// In this library variables are represented by a zero-based index
	// Here we will say x = variable 0 and y = variable 1

	// This adds 5x^2 to the polynomial
	poly.addTerm(5, {0,0});

	// This adds 3xy to the polynomial
	poly.addTerm(3, {0,1});

	// This adds -2 to the polynomial
	poly.addTerm(-2, {});

	// Output the polynomial
	std::cout << poly << std::endl;

	// Output the polynomial squared
	std::cout << poly*poly << std::endl;

}

