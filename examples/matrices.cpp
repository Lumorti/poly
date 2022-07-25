#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Create a 2x2 matrix of polynomials with a max of 5 variables
	// (  x^2    2y  )
	// (   5     xy  )
	PolynomialMatrix<float> mat(5, 2, 2);
	mat[0][0].addTerm(1, {0,0});
	mat[1][0].addTerm(5, {});
	mat[0][1].addTerm(2, {1});
	mat[1][1].addTerm(1, {0,1});

	// Create a 2-vector of polynomials with a max of 5 variables
	// (  x  )
	// (  y  )
	PolynomialMatrix<float> vec(5, 2);
	vec[0][0].addTerm(1, {0});
	vec[1][0].addTerm(1, {1});

	// Output both
	std::cout << mat << std::endl;
	std::cout << vec << std::endl;

	// Output the product
	std::cout << mat*vec << std::endl;

}

