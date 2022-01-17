#include <fstream>

// MOSEK
#include "fusion.h"

// Standard cpp entry point
int main (int argc, char ** argv) {

	// Load the data from file
	std::cout << "Loading file..." << std::endl;
	std::ifstream infile("d2n2.csv");
	std::vector<int> xLocs;
	std::vector<int> yLocs;
	std::vector<double> coeffs;
	double val;
	int x, y;
	int maxI = 0;
	int maxJ = 0;
	while (infile >> x >> y >> val) {
		xLocs.push_back(x);
		yLocs.push_back(y);
		coeffs.push_back(val/2.0);
		xLocs.push_back(y);
		yLocs.push_back(x);
		coeffs.push_back(val/2.0);
	}
	int arraySize = *max_element(xLocs.begin(), xLocs.end())+1;

	// Minimize this array TODO
	int sizeBefore = arraySize;
	int i=0;
	int colsRemoved=0;
	int rowsRemoved=0;
	while (i < sizeBefore) {
		if (std::find(xLocs.begin(), xLocs.end(), i) == xLocs.end()) {
			for (int j=0; j<xLocs.size(); j++) {
				if (xLocs[j] > i) {
					xLocs[j] -= 1;
				}
			}
			colsRemoved++;
		}
		i++;
	}
	i=0;
	while (i < sizeBefore) {
		if (std::find(yLocs.begin(), yLocs.end(), i) == yLocs.end()) {
			for (int j=0; j<yLocs.size(); j++) {
				if (yLocs[j] > i) {
					yLocs[j] -= 1;
				}
			}
			rowsRemoved++;
		}
		i++;
	}
	arraySize = *max_element(xLocs.begin(), xLocs.end())+1;
	std::cout << "array size before = " << sizeBefore << " x " << sizeBefore << std::endl;
	std::cout << "array size after = " << arraySize << " x " << arraySize << std::endl;
	std::cout << "cols removed = " << colsRemoved << std::endl;
	std::cout << "rows removed = " << rowsRemoved << std::endl;

	// Turn it into MOSEK form
	std::cout << "Creating model..." << std::endl;
	auto C = mosek::fusion::Matrix::sparse(arraySize, arraySize, monty::new_array_ptr(xLocs), monty::new_array_ptr(yLocs), monty::new_array_ptr(coeffs));

	// Create the MOSEK model 
	mosek::fusion::Model::t model = new mosek::fusion::Model(); 
	model->setLogHandler([=](const std::string & msg) {std::cout << msg << std::flush;});

	// The matrices to optimise
	auto dimRef = monty::new_array_ptr(std::vector<int>({arraySize, arraySize}));
	mosek::fusion::Variable::t X = model->variable(mosek::fusion::Domain::inPSDCone(arraySize));

	// Objective
	model->objective(mosek::fusion::ObjectiveSense::Minimize, mosek::fusion::Expr::dot(C, X));

	// Top left element of X should be 1
	model->constraint(X->index(0,0), mosek::fusion::Domain::equalsTo(1.0));

	// Solve the SDP
	std::cout << "Solving model..." << std::endl;
	model->solve();

	// Output the result
	double res = model->primalObjValue();
	std::cout << res << std::endl;

}
