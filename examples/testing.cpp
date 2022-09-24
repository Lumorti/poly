#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// min x_1 + 2x_2 - 3x_3
	// s.t. -1 <= x_i <= 1
	// min = -6 at -1,-1,1
	Polynomial<double> obj(3);
	obj.addTerm(1, {0});
	obj.addTerm(2, {1});
	obj.addTerm(-3, {2});

	std::vector<Polynomial<double>> cons;
	//{
		//Polynomial<double> con;
		//con.addTerm(1, {});
		//con.addTerm(-1, {0});
		//cons.push_back(con);
	//}

	// Create the problem
	PolynomialBinaryProblem<double> prob(obj, {}, cons);
	prob.lowerBoundNew(20);

}
