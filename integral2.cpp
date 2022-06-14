#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <chrono>
#include <math.h>
#include <unordered_map>

// Use Eigen for matrix/vector ops
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

// Some useful definitions
using namespace std::complex_literals;
typedef std::unordered_map<std::string, std::complex<double>> eqn;

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

// Remove zero terms from an equation
eqn prune(eqn a) {
	eqn aPruned;
    for (auto const &pair: a) {
		if (std::abs(pair.second) > 1e-8) {
			aPruned[pair.first] = pair.second;
		}
    }
	return aPruned;
}

// Check equivalence of two equations
bool equiv(eqn a, eqn b) {
    for (auto const &pair: a) {
		if (b[pair.first] != pair.second) {
			return false;
		}
	}
	return true;
}

// Integrate a polynomial with respect to a variable index
eqn integrate(eqn a, int ind) {

	// Process each term
	eqn result;
    for (auto const &pair: a) {

		// Count the number of times the variable appears here
		int degree = std::count(pair.first.begin(), pair.first.end(), ind);

		// x**3 -> 0.25*x**4
		std::vector<int> newInd = pair.first;
		newInd.push_back(ind);

		// Add to the result
		if (result.find(newInd) != result.end()) {
			result[newInd] += pair.second / (degree+1);
		} else {
			result[newInd] = pair.second / (degree+1);
		}

	}

	return result;

}

// Differentiate a polynomial with respect to a variable index
eqn differentiate(eqn a, int ind) {

	// Process each term
	eqn result;
    for (auto const &pair: a) {

		// Count the number of times the variable appears here
		int degree = std::count(pair.first.begin(), pair.first.end(), ind);

		// x**3 -> 3*x**2
		std::vector<int> newInd = pair.first;
		int indFound = -1;
		for (int j=newInd.size()-1; j>=0; j--) {
			if (newInd[j] == ind) {
				indFound = j;
				break;
			}
		}
		if (indFound >= 0) {
			newInd.erase(newInd.begin()+indFound);
		}

		// Add to the result
		if (result.find(newInd) != result.end()) {
			result[newInd] += pair.second * degree;
		} else {
			result[newInd] = pair.second * degree;
		}

	}

	return prune(result);

}

// Eval an equation
std::complex<double> eval(eqn a, std::vector<double> x) {
	std::complex<double> soFar = 0;
    for (auto const &pair: a) {
		std::complex<double> sub = 1;
		for (int j=0; j<pair.first.size(); j++) {
			sub *= x[pair.first[j]];
		}
		soFar += sub*pair.second;
	}
	return soFar;
}

// Eval an equation
std::complex<double> eval(eqn a, Eigen::VectorXd x) {
	std::complex<double> soFar = 0;
    for (auto const &pair: a) {
		std::complex<double> sub = 1;
		for (int j=0; j<pair.first.size(); j++) {
			sub *= x[pair.first[j]];
		}
		soFar += sub*pair.second;
	}
	return soFar;
}

// Conjugate an equation
eqn conj(eqn a) {
	eqn aConj;
    for (auto const &pair: a) {
		aConj[pair.first] = std::conj(pair.second);
	}
	return aConj;
}

// Add two equations
eqn add(eqn a, eqn b) {

	// Start with one equation
	eqn result = a;

	// For each term of the other
    for (auto const &pair: b) {

		// If it's new add it, otherwise combine with the existing
		if (result.find(pair.first) != result.end()) {
			result[newInd] += pair.second;
		} else {
			result[newInd] = pair.second;
		}

	}

	// Remove any zeros and return
	return prune(result);

}

// Multiply two equations
eqn multiply(eqn a, eqn b) {

	// For each term of both equations
	eqn result;
    for (auto const &pair1: a) {
		for (auto const &pair2: b) {

			// Combine the term list
			std::vector<int> combined;
			combined.insert(combined.end(), pair1.second.begin(), pair1.second.end());
			combined.insert(combined.end(), pair2.second.begin(), pair2.second.end());

			// Sort the term list
			std::sort(combined.begin(), combined.end());

			// If it's new add it, otherwise combine with the existing
			if (result.find(combined) != result.end()) {
				result[combined] += pair1.second*pair2.second;
			} else {
				result[combined] = pair1.second*pair2.second;
			}

		}
	}

	// Remove any zeros and return
	return prune(result);

}

// Print an equation
void pretty(eqn a) {
	for (auto const &pair: a) {
		std::cout << pair.second << "*{";
		for (int j=0; j<pair.first.size(); j++) {
			std::cout << pair.first[j];
			if (j < pair.first.size()-1) {
				std::cout << ", ";
			}
		}
		std::cout << "} + ";
    }
	std::cout << " = 0" << std::endl;
}

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Get the problem from the args
	int d = 2;
	int n = 2;
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}

	// Useful quantities
	int numVarsNonConj = n*d*d;
	int numVars = 2*n*d*d;
	int conjDelta = numVarsNonConj;
	double rt2 = 1.0/std::sqrt(2.0);

	// Known ideal TODO
	std::vector<double> idealX(numVars);
	if (d == 2 && n == 2) {
		idealX = {1, 0, 0, 1, rt2, rt2, rt2, -rt2, 0, 0, 0, 0, 0,   0,   0,    0};
	} else if (d == 2 && n == 3) {
		idealX = {1, 0, 0, 1, rt2, rt2, rt2, -rt2,  rt2, 0,   rt2, 0,
									  0, 0, 0, 0, 0,   0,   0,    0,    0,   rt2, 0,   -rt2};
	}

	// The list of equations to fill
	std::vector<eqn> eqns;

	// Generate orthogonality equations
	std::cout << "Generating equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			for (int k=0; k<d; k++) {
				for (int l=0; l<d; l++) {

					// Get the equation before squaring
					eqn eqnPreSquare = {};
					for (int m=0; m<d; m++) {
						int var1 = i*d*d + k*d + m;
						int var2 = j*d*d + l*d + m;
						int var3 = var1 + conjDelta;
						int var4 = var2 + conjDelta;
						eqnPreSquare.push_back(term(1, {var1, var2}));
						eqnPreSquare.push_back(term(1, {var3, var4}));
						eqnPreSquare.push_back(term(1i, {var1, var4}));
						eqnPreSquare.push_back(term(-1i, {var2, var3}));
					}

					// For the normalisation equations
					if (i == j && k == l) {
						eqnPreSquare.push_back(term(-1.0, {}));
					}

					// Square this equation
					eqn eqn = multiply(conj(eqnPreSquare), eqnPreSquare);

					// For the MUB-ness equations
					if (i != j) {
						eqn.push_back(term(-1.0/d, {}));
						eqn = multiply(conj(eqn), eqn);
					}

					// Add it to the list
					//pretty(eqn);
					//std::cout << "eval = " << eval(eqn, idealX) << std::endl;
					//std::cout << std::endl;
					eqns.push_back(eqn);

				}
			}
		}
	}

	//// Combine these to create a single polynomial
	std::cout << "Creating single polynomial..." << std::endl;
	eqn poly;
	for (int i=0; i<eqns.size(); i++) {
		poly = add(poly, eqns[i]);
	}
	//pretty(poly);
	std::cout << "poly eval = " << eval(poly, idealX) << std::endl;
	std::cout << std::endl;

	// Integrate by a variable
	std::cout << "Integrating..." << std::endl;
	eqn integrated = integrate(poly, 0);
	//pretty(integrated);
	//std::cout << std::endl;

	// Get the derivates by each vars 
	std::cout << "Differentiating..." << std::endl;
	std::vector<eqn> gradient(numVars);
	std::vector<std::vector<eqn>> hessian(numVars, std::vector<eqn>(numVars));
	for (int i=0; i<numVars; i++) {
		gradient[i] = differentiate(integrated, i);
	}
	for (int i=0; i<numVars; i++) {
		for (int j=0; j<numVars; j++) {
			hessian[i][j] = differentiate(gradient[i], j);
		}
	}

	// Ensure the first element is the same as the orig
	if (!equiv(gradient[0], poly)) {
		std::cout << "ERROR - differential of integral is not equivalent!" << std::endl;
		return 0;
	}

	std::cout << "Optimizing..." << std::endl;

	// Random starting x
	//std::srand(unsigned(std::time(nullptr)));
	std::srand(0);
	Eigen::VectorXd x = Eigen::VectorXd::Random(numVars);
	//x(0) = 1e-29;
	//for (int i=0; i<x.size(); i++) {
		//x(i) = double(std::rand()) / RAND_MAX;
	//}
	std::cout << "init x = " << std::endl;
	for (int i=0; i<x.size(); i++) {
		std::cout << x(i) << ", ";
	}
	std::cout << std::endl;

	// Initial Hessian calc
	Eigen::MatrixXd H(numVars, numVars);
	//for (int i=0; i<numVars; i++) {
		//for (int j=0; j<numVars; j++) {
			//H(i,j) = std::real(eval(hessian[i][j], x));
		//}
	//}

	// Initial gradient calc
	Eigen::VectorXd g(numVars);
	//Eigen::VectorXd gNew(numVars);
	//for (int i=0; i<numVars; i++) {
		//g(i) = std::real(eval(gradient[i], x));
	//}

	// Perform gradient descent using this info
	double alpha = 0.1;
	if (argc > 3) {
		alpha = std::stod(argv[3]);
	}
	Eigen::MatrixXd inv(numVars, numVars);
	Eigen::VectorXd p(numVars);
	Eigen::VectorXd q(numVars);
	Eigen::VectorXd s(numVars);
	Eigen::VectorXd y(numVars);
	double maxX = 0;
	for (int iter=0; iter<10000000; iter++) {

		// Calculate the gradient
		for (int i=0; i<numVars; i++) {
			g(i) = std::real(eval(gradient[i], x));
		}

		// Calculate the Hessian
		for (int i=0; i<numVars; i++) {
			for (int j=0; j<numVars; j++) {
				H(i,j) = std::real(eval(hessian[i][j], x));
			}
		}

		// Calculate the inverse
		inv = H.inverse();

		// Calculate the direction
		p = -inv*g;
		//p = -H.colPivHouseholderQr().solve(g);
		//Eigen::LeastSquaresConjugateGradient<Eigen::MatrixXd> lscg;
		//lscg.compute(H);
		//p = -lscg.solve(g);

		// Perform a line search
		//alpha = 0.1;
		//double c = 0.0;
		//double tau = 0.9;
		//double m = g.dot(p);
		//double t = -c*m;
		//for (int i=0; i<10; i++) {

			//// Keep going until the Wolfe conditions hold TODO
			//if (std::real(eval(integrated, x)) - std::real(eval(integrated, x+alpha*p)) >= t*alpha) {
				//break;
			//}

			//// Reduce alpha
			//alpha *= tau;

		//}

		// Perform the update
		x += alpha*p;

		// Calculate the gradient at the new point
		//for (int i=0; i<numVars; i++) {
			//gNew(i) = std::real(eval(gradient[i], x));
		//}

		// BFGS-style Hessian update TODO
		//s = alpha*p;
		//y = gNew - g;
		//g = gNew;
		//H += (y*y.transpose()) / (y.dot(s)) - ((H*s*s.transpose()*H.transpose()) / (s.dot(H*s)));

		// Per-iteration output
		std::cout << iter << " " << std::real(eval(integrated, x)) << " " << g.norm() << " " << alpha << std::endl;

		// Convergence criteria
		if (std::abs(g(0)) < 1e-10) {
			break;
		}

	}

	// Output the final result
	std::cout << "final x = " << std::endl;
	for (int i=0; i<x.size(); i++) {
		std::cout << x(i) << ", ";
	}
	std::cout << std::endl;

	return 0;

}
	
