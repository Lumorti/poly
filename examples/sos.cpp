#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Create a polynomial known to be SOS
	// f(x,y) = (x-3)^2 + (y-x)^2 + (yyy)^2
	//        = 2xx - 6x + 9 + yy - 2yx + yyyyyy
	Polynomial<float> polySOS(3);
	polySOS.addTerm(2, {0,0});
	polySOS.addTerm(-6, {0});
	polySOS.addTerm(9, {});
	polySOS.addTerm(1, {1,1});
	polySOS.addTerm(-2, {0,1});
	polySOS.addTerm(1, {1,1,1,1,1,1});

	// Create a polynomial known to not be SOS
	// f(x,y) = 1 + xx yyyy + xxxxyy - 3xxyy
	Polynomial<float> polyNotSOS(3);
	polyNotSOS.addTerm(1, {});
	polyNotSOS.addTerm(1, {0,0,1,1,1,1});
	polyNotSOS.addTerm(1, {0,0,0,0,1,1});
	polyNotSOS.addTerm(-3, {0,0,1,1});

	// Check to make sure we're right
	std::cout << polySOS << std::endl;
	std::cout << "should be 1: " << polySOS.isSOS(0) << std::endl;
	std::cout << "should be 0: " << polyNotSOS.isSOS(0) << std::endl;

}

