#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <chrono>
#include <math.h>

// Use Eigen for matrix/vector ops
#include <Eigen/Dense>

// Some useful definitions
using namespace std::complex_literals;
typedef std::pair<std::complex<double>, std::vector<int>> term;
typedef std::vector<term> eqn;

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
	for (int i=0; i<a.size(); i++) {
		if (std::abs(a[i].first) > 1e-8) {
			aPruned.push_back(a[i]);
		}
	}
	return aPruned;
}

// Check equivalence of two equations
bool equiv(eqn a, eqn b) {

	// If there aren't the same number of non-zero term, they aren't equal
	if (a.size() != b.size()) {
		std::cout << "wrong size" << std::endl;
		return false;
	}

	// For now check that they're the same order
	for (int i=0; i<a.size(); i++) {
		if (a[i].first != b[i].first) {
			std::cout << "diff coeff = " << i << std::endl;
			return false;
		}
		if (a[i].second != b[i].second) {
			std::cout << "diff indices = " << i << std::endl;
			return false;
		}
	}

	return true;

}

// Integrate a polynomial with respect to a variable index
eqn integrate(eqn a, int ind) {

	// Process each term
	eqn result;
	for (int i=0; i<a.size(); i++) {
		term newTerm = a[i];

		// Count the number of times the variable appears here
		int degree = std::count(a[i].second.begin(), a[i].second.end(), ind);

		// x**3 -> 0.25*x**4
		newTerm.first *= 1.0 / (degree+1);
		newTerm.second.push_back(ind);

		// Add to the result
		result.push_back(newTerm);

	}

	return result;

}

// Differentiate a polynomial with respect to a variable index
eqn differentiate(eqn a, int ind) {

	// Process each term
	eqn result;
	for (int i=0; i<a.size(); i++) {
		term newTerm = a[i];

		// Count the number of times the variable appears here
		int degree = std::count(a[i].second.begin(), a[i].second.end(), ind);

		// x**3 -> 3*x**2
		newTerm.first *= degree;
		int indFound = -1;
		for (int j=newTerm.second.size()-1; j>=0; j--) {
			if (newTerm.second[j] == ind) {
				indFound = j;
				break;
			}
		}
		if (indFound >= 0) {
			newTerm.second.erase(newTerm.second.begin()+indFound);
		}

		// Add to the result
		result.push_back(newTerm);

	}

	return prune(result);

}


// Eval an equation
std::complex<double> eval(eqn a, std::vector<double> x) {
	std::complex<double> soFar = 0;
    for (int i=0; i<a.size(); i++) {
		std::complex<double> sub = 1;
		for (int j=0; j<a[i].second.size(); j++) {
			sub *= x[a[i].second[j]];
		}
		soFar += sub*a[i].first;
	}
	return soFar;
}

// Eval an equation
std::complex<double> eval(eqn a, Eigen::VectorXd x) {
	std::complex<double> soFar = 0;
    for (int i=0; i<a.size(); i++) {
		std::complex<double> sub = 1;
		for (int j=0; j<a[i].second.size(); j++) {
			sub *= x[a[i].second[j]];
		}
		soFar += sub*a[i].first;
	}
	return soFar;
}

// Conjugate an equation
eqn conj(eqn a) {
	eqn aConj;
	for (int i=0; i<a.size(); i++) {
		aConj.push_back(a[i]);
		aConj[i].first = std::conj(aConj[i].first);
	}
	return aConj;
}

// Add two equations
eqn add(eqn a, eqn b) {

	// Start with one equation
	eqn result = a;

	// For each term of the other
	for (int i=0; i<b.size(); i++) {

		// See if it's already in the equation
		int indFound = -1;
		for (int j=0; j<result.size(); j++) {
			if (result[j].second == b[i].second) {
				indFound = j;
				break;
			}
		}

		// If it's new add it, otherwise combine with the existing
		if (indFound == -1) {
			result.push_back(b[i]);
		} else {
			result[indFound].first += b[i].first;
		}

	}

	// Remove any zeros and return
	return prune(result);

}

// Multiply two equations
eqn multiply(eqn a, eqn b) {

	// For each term of both equations
	eqn result;
	for (int t1=0; t1<a.size(); t1++) {
		for (int t2=0; t2<b.size(); t2++) {

			// Combine the term list
			std::vector<int> combined;
			combined.insert(combined.end(), a[t1].second.begin(), a[t1].second.end());
			combined.insert(combined.end(), b[t2].second.begin(), b[t2].second.end());

			// Sort the term list
			std::sort(combined.begin(), combined.end());

			// See if it's already in the equation
			int indFound = -1;
			for (int i=0; i<result.size(); i++) {
				if (result[i].second == combined) {
					indFound = i;
					break;
				}
			}

			// If it's new add it, otherwise combine with the existing
			if (indFound == -1) {
				result.push_back(term(a[t1].first*b[t2].first, combined));
			} else {
				result[indFound].first += a[t1].first*b[t2].first;
			}

		}
	}

	// Remove any zeros and return
	return prune(result);

}

// Print an equation
void pretty(eqn a) {
    for (int i=0; i<a.size(); i++) {
		std::cout << a[i].first << "*{";
		for (int j=0; j<a[i].second.size(); j++) {
			std::cout << a[i].second[j];
			if (j < a[i].second.size()-1) {
				std::cout << ", ";
			}
		}
		std::cout << "}";
		if (i < a.size()-1) {
			std::cout << " + ";
		}
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
	pretty(poly);
	std::cout << eval(poly, idealX) << std::endl;
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
			hessian[i][j] = differentiate(differentiate(integrated, i), j);
		}
	}

	// Ensure the first element is the same as the orig
	if (!equiv(gradient[0], poly)) {
		std::cout << "ERROR - differential of integral is not equivalent!" << std::endl;
		return 0;
	}

	// Random starting x
	Eigen::VectorXd x(numVars);
	//std::srand(unsigned(std::time(nullptr)));
	std::srand(0);
	for (int i=0; i<x.size(); i++) {
		x(i) = double(std::rand()) / RAND_MAX;
	}
	for (int i=0; i<x.size(); i++) {
		std::cout << x(i) << ", ";
	}
	std::cout << std::endl;

	// Perform gradient descent using this info
	Eigen::VectorXd g(numVars);
	Eigen::MatrixXd H(numVars, numVars);
	Eigen::MatrixXd inv(numVars, numVars);
	Eigen::VectorXd p(numVars);
	Eigen::VectorXd q(numVars);
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
		//p = H.householderQr().solve(-g);
		//p = -g;

		// Perform a line search
		double alpha = 0.1;
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

		// TODO bounded, maybe with a barrier function?
		bool shouldBreak = false;
		for (int j=0; j<numVars; j++) {
			if (x(j) > maxX) {
				maxX = x(j);
			}
			//if (x(j) > 10.0 || x(j) < -10.0) {
				//shouldBreak = true;
				//break;
			//}
		}

		// Per-iteration output
		std::cout << iter << " " << std::real(eval(integrated, x)) << " " << g.norm() << " " << alpha << " " << maxX << std::endl;

		// Convergence criteria
		if (shouldBreak || g(0) < 1e-10) {
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
	
