#include "poly.h"
#include <time.h>

// Generalised Pauli X matrix
Eigen::MatrixXcd generalPauliX(int d) {
	Eigen::MatrixXcd generalX = Eigen::MatrixXcd::Zero(d, d);
	for (int i=0; i<d; i++) {
		generalX((i+1)%d, i) = 1;
	}
	return generalX;
}

// Generalised Pauli Z matrix
Eigen::MatrixXcd generalPauliZ(int d) {
	Eigen::MatrixXcd generalZ = Eigen::MatrixXcd::Zero(d, d);
	double s = (d-1) / 2.0;
	for (int i=0; i<d; i++) {
		generalZ(i, i) = (s-i);
	}
	return generalZ;
}

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Sometimes we want this, sometimes we don't
	std::srand(std::time(NULL));

	// Get the problem from the args
	int d = 2;
	int n = 2;
	int level = 2;
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}
	int knownBases = n-1;
	if (argc > 3) {
		knownBases = std::stoi(argv[3]);
	}
	if (argc > 4) {
		level = std::stoi(argv[4]);
	}

	// Useful quantities
	int numVarsNonConj = n*d*d;
	int numVars = 2*numVarsNonConj;
	int conjDelta = numVarsNonConj;
	double rt2 = 1.0/std::sqrt(2.0);

	// The list of equations to fill
	PolynomialSystem<double> eqns;

	// Generate equations
	std::cout << "Generating equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			for (int k=0; k<d; k++) {
				for (int l=0; l<d; l++) {

					// (a+ib)*(c+id)
					Polynomial<std::complex<double>> eqn(numVars);
					for (int m=0; m<d; m++) {
						int var1 = i*d*d + k*d + m;
						int var2 = j*d*d + l*d + m;
						int var3 = var1 + conjDelta;
						int var4 = var2 + conjDelta;
						eqn.addTerm(1, {var1, var2});
						eqn.addTerm(1, {var3, var4});
						eqn.addTerm(1i, {var1, var4});
						eqn.addTerm(-1i, {var2, var3});
					}

					// For the normalisation equations
					if (i == j && k == l) {
						eqn.addTerm(-1.0, {});
					}

					// For the MUB-ness equations
					if (i != j) {
						eqn = eqn.conjugate() * eqn;
						eqn.addTerm(-1.0/d, {});
	 				}

					// Both the real and imag parts should be 0
					eqns.addPolynomial(Polynomial<double>(eqn.real()));
					eqns.addPolynomial(Polynomial<double>(eqn.imag()));

				}
			}
		}
	}

	// For replacing variables with known ideal values
	std::vector<std::vector<double>> toAdd = {};

	// General known bases (eigenbasis of Z, X, XZ) TODO
	std::vector<double> valsToReplaceReal = {};
	std::vector<double> valsToReplaceImag = {};
	std::vector<Eigen::VectorXcd> set1;
	std::vector<Eigen::VectorXcd> set2;
	if (knownBases >= 1) {
		Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver;
		eigensolver.compute(generalPauliZ(d));
		Eigen::MatrixXcd vecs = eigensolver.eigenvectors();
		std::cout << generalPauliZ(d) << std::endl << std::endl;
		std::cout << vecs << std::endl << std::endl;
		for (int j=0; j<vecs.cols(); j++) {
			set1.push_back(vecs.col(j));
			for (int i=0; i<vecs.rows(); i++) {
				valsToReplaceReal.push_back(vecs(i,j).real());
				valsToReplaceImag.push_back(vecs(i,j).imag());
			}
		}
	} 
	if (knownBases >= 2) {
		Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver;
		eigensolver.compute(generalPauliX(d));
		Eigen::MatrixXcd vecs = eigensolver.eigenvectors();
		std::cout << generalPauliX(d) << std::endl << std::endl;
		std::cout << vecs << std::endl << std::endl;
		for (int j=0; j<vecs.cols(); j++) {
			set2.push_back(vecs.col(j));
			for (int i=0; i<vecs.rows(); i++) {
				valsToReplaceReal.push_back(vecs(i,j).real());
				valsToReplaceImag.push_back(vecs(i,j).imag());
			}
		}
	} 
	if (knownBases >= 3) {
		Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver;
		eigensolver.compute(generalPauliX(d)*generalPauliZ(d));
		Eigen::MatrixXcd vecs = eigensolver.eigenvectors();
		std::cout << vecs << std::endl << std::endl;
		for (int j=0; j<vecs.cols(); j++) {
			for (int i=0; i<vecs.rows(); i++) {
				valsToReplaceReal.push_back(vecs(i,j).real());
				valsToReplaceImag.push_back(vecs(i,j).imag());
			}
		}
	}

	// Reals at the start, imaginaries at the end
	std::vector<double> valsToReplace= {};
	valsToReplace.insert(valsToReplace.end(), valsToReplaceReal.begin(), valsToReplaceReal.end());
	valsToReplace.insert(valsToReplace.end(), valsToReplaceImag.begin(), valsToReplaceImag.end());

	// d == 2
	//if (d == 2) {
		//if (knownBases >= 1) {
			//toAdd.push_back({1, 0, 0, 1, 0, 0, 0, 0});
		//} 
		//if (knownBases >= 2) {
			//toAdd.push_back({rt2, rt2, rt2, -rt2, 0, 0, 0, 0});
		//} 
		//if (knownBases >= 3) {
			//toAdd.push_back({rt2, 0, rt2, 0, 0, rt2, 0, -rt2});
		//}
	//}

	// For d == 4
	//if (d == 4) {
		//if (knownBases >= 1) {
			//toAdd.push_back({1, 0, 0, 0,    0, 1, 0, 0,    0, 0, 1, 0,   0, 0, 0, 1,      
							 //0, 0, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,  0, 0, 0, 0,});
		//} 
		//if (knownBases >= 2) {
			//toAdd.push_back({0.5, 0.5, 0.5, 0.5,    0.5, 0.5, -0.5, -0.5,    0.5, -0.5, -0.5, 0.5,   0.5, -0.5, 0.5, -0.5,      
							 //0, 0, 0, 0,    0, 0, 0, 0,   0, 0, 0, 0,  0, 0, 0, 0,});
		//} 
		//if (knownBases >= 3) {
			//toAdd.push_back({0.5, -0.5, 0, 0,    0.5, -0.5, 0, 0,    0.5, 0.5, 0, 0,   	0.5, 0.5, 0, 0,      
							 //0, 0, -0.5, -0.5,     0, 0, 0.5, 0.5,    0, 0, 0.5, -0.5,   0, 0, -0.5, 0.5,});
		//}
		//if (knownBases >= 4) {
			//toAdd.push_back({0.5, 0, 0, -0.5,    0.5, 0, 0, 0.5,    0.5, 0, 0, -0.5,   0.5, 0, 0, 0.5,    
							   //0, -0.5, -0.5, 0,   0, -0.5, 0.5, 0,   0, 0.5, 0.5, 0,   0, 0.5, -0.5, 0,});
		//}
		//if (knownBases >= 5) {
			//toAdd.push_back({0.5, 0, -0.5, 0,    0.5, 0, 0.5, 0,    0.5, 0, -0.5, 0,   0.5, 0, 0.5, 0,    
							   //0, -0.5, 0, -0.5,   0, -0.5, 0, 0.5,    0, 0.5, 0, 0.5,   0, 0.5, 0, -0.5,});
		//}
	//}

	// For d == 6
	//if (d == 6) {
		//if (knownBases >= 1) {
			//toAdd.push_back({1,0,0,0,0,0,  0,1,0,0,0,0,  0,0,1,0,0,0,  0,0,0,1,0,0,  0,0,0,0,1,0,  0,0,0,0,0,1,  0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0});
		//} 
		//if (knownBases >= 2) {
			////toAdd.push_back({});
		//} 
		//if (knownBases >= 3) {
			////toAdd.push_back({});
		//}
	//}

	// Combine the real and imag vecs
	//std::vector<double> valsToReplace = {};
	//std::vector<double> valsToReplaceReal = {};
	//std::vector<double> valsToReplaceImag = {};
	//for (int i=0; i<toAdd.size(); i++) {
		//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd[i].begin(), toAdd[i].begin()+d*d);
		//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd[i].begin()+d*d, toAdd[i].end());
	//}
	//valsToReplace.insert(valsToReplace.end(), valsToReplaceReal.begin(), valsToReplaceReal.end());
	//valsToReplace.insert(valsToReplace.end(), valsToReplaceImag.begin(), valsToReplaceImag.end());

	// The indices to change
	std::vector<int> indsToReplace = {};
	for (int i=0; i<valsToReplaceReal.size(); i++) {
		indsToReplace.push_back(i);
	}
	for (int i=0; i<valsToReplaceImag.size(); i++) {
		indsToReplace.push_back(i+numVarsNonConj);
	}

	// Replace everything
	PolynomialSystem<double> simpEqns = eqns.substitute(indsToReplace, valsToReplace).simplify();

	// Only keep equations with the first basis element (i.e. d*2 elements)
	std::vector<int> vars = simpEqns.getVariables();
	int conjDeltaTemp = vars.size() / 2;
	std::vector<int> toKeep;
	std::vector<int> toRemove;
	for (int i=0; i<vars.size(); i++) {
		if ((vars[i] < d) || (vars[i] >= conjDeltaTemp && vars[i] < conjDeltaTemp+d)) {
			toKeep.push_back(i);
		}
	}
	//PolynomialSystem<double> reducedEqns = simpEqns.withVariables(toKeep).removeDuplicates().simplify().prune();
	//PolynomialSystem<double> reducedEqns = simpEqns.withOnlyVariables(toKeep).removeDuplicates().simplify().prune();
	PolynomialSystem<double> reducedEqns = simpEqns.removeDuplicates().simplify().prune();

	// Output the system
	std::cout << "has solutions if MUBs exist:" << std::endl;
	std::cout << reducedEqns << std::endl;
	if (reducedEqns.size() == 0) {
		return 0;
	}

	// Try to prove this is infeasible with Hillbert's Nullstellenstatz TODO
	reducedEqns.proveInfeasible(level);

	// Combine these to create a single polynomial
	//std::cout << "Final poly has " << mapTo.size() << " variables" << std::endl;
	//std::cout << "Creating single polynomial..." << std::endl;
	//Polynomial<double> poly(numVars);
	//for (int i=0; i<newEqns.size(); i++) {
		//poly += newEqns[i]*newEqns[i];
	//}

	// Find a root of this
	//std::cout << "Attempting to find a root..." << std::endl;
	//poly.numVars += 1;
	//std::vector<double> x = poly.findRoot(poly.numVars-1, 0.5, 1e-10, 10000000);
	//std::vector<double> x = poly.findRoot(2, 0.01, 1e-10, 10000000);
	//std::cout << "Testing this x = " << poly.eval(x) << std::endl;
	//for (int i=0; i<newEqns.size(); i++) {
		//std::cout << "Testing eqn = " << newEqns[i].eval(x) << std::endl;
	//}

	// Get the complex relaxation of this
	//std::cout << "Attempting to find a relaxed root..." << std::endl;
	//Polynomial<double> relaxed = poly.getComplexRelaxation();
	//int origVars = poly.numVars;
	//relaxed.numVars += 1;
	//std::vector<double> x2 = relaxed.findRoot(relaxed.numVars-1);
	//std::vector<std::complex<double>> xComp(x2.size());
	//for (int i=0; i<origVars; i++) {
		//xComp[i] = std::complex<double>(x2[i], x2[i+origVars]);
		//std::cout << xComp[i] << std::endl;
	//}
	//std::cout << "Testing this x = " << relaxed.eval(x2) << std::endl;
	//std::cout << "Testing this x = " << Polynomial<std::complex<double>>(poly).eval(xComp) << std::endl;
	//for (int i=0; i<newEqns.size(); i++) {
		//std::cout << "Testing eqn = " << Polynomial<std::complex<double>>(newEqns[i]).eval(xComp) << std::endl;
	//}

	return 0;

}
	
