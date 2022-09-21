#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// min x_1 + x_2 
	Polynomial<double> obj(2);
	obj.addTerm(1, {0});
	obj.addTerm(1, {1});

	// x_1 + 1 >= 0
	// x_2 + 1 >= 0
	std::vector<Polynomial<double>> cons;
	{
		Polynomial<double> con;
		con.addTerm(1, {0});
		con.addTerm(1, {});
		cons.push_back(con);
	}
	{
		Polynomial<double> con;
		con.addTerm(1, {1});
		con.addTerm(1, {});
		cons.push_back(con);
	}

	// Create the problem
	PolynomialBinaryProblem<double> prob(obj, {}, cons);
	prob.lowerBoundNew(100);

}
