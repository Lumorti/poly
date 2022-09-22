#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// min -x_1 -x_2 
	// s.t. -1 <= x_i <= 1
	// min = -2 at 1,1
	Polynomial<double> obj(2);
	obj.addTerm(-1, {0});
	obj.addTerm(-1, {1});

	// min -(2x_1-1) -(2x_2-1)
	// s.t. 0 <= x <= 1
	// min = -2 at 1,1
	
	// c.x + d
	// c.(2x-1) + d = 2cx-c1+d

	std::vector<Polynomial<double>> cons;

	// x <= 1
	// 1 - x >= 0
	//{
		//Polynomial<double> con;
		//con.addTerm(1, {});
		//con.addTerm(-1, {0});
		//cons.push_back(con);
	//}
	//{
		//Polynomial<double> con;
		//con.addTerm(1, {});
		//con.addTerm(-1, {1});
		//cons.push_back(con);
	//}

	// Create the problem
	PolynomialBinaryProblem<double> prob(obj, {}, cons);
	prob.lowerBoundNew(100);

}
