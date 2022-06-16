#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <chrono>
#include <math.h>
#include <unordered_map>

// Use Eigen for matrix/vector ops
#include <Eigen/Dense>

// Some useful definitions
using namespace std::complex_literals;

// Class allowing manipulation of polynomials
class Polynomial {
public:

	// General poly properties
	int numVars = 1;
	int digitsPerInd = 1;
	std::unordered_map<std::string, std::complex<double>> coeffs;

	// Used for fast eval
	std::vector<std::vector<int>> inds;
	std::vector<double> vals;

	// Default constructor
	Polynomial(int numVars_) {
		numVars = numVars_;
		digitsPerInd = std::ceil(std::log10(numVars));
	}

	// Add a term given a coefficient and an index list
	void addTerm(std::complex<double> coeff, std::vector<int> list) {

		// Convert the vector to a string
		std::string asString = "";
		for (int i=0; i<list.size(); i++) {
			std::string padded = std::to_string(list[i]);
			padded.insert(0, digitsPerInd-padded.size(), ' ');
			asString += padded;
		}

		// Add it, combining it if the term already exists
		if (coeffs.find(asString) != coeffs.end()) {
			coeffs[asString] += coeff;
			if (std::abs(coeffs[asString]) < 1e-10) {
				coeffs.erase(asString);
			}
		} else {
			coeffs[asString] = coeff;
		}

	}

	// Remove zeros from a polynomial
	Polynomial prune() {
		Polynomial newPoly(numVars);
		for (auto const &pair: coeffs) {
			if (std::abs(pair.second) > 1e-8) {
				newPoly.coeffs[pair.first] = pair.second;
			}
		}
		return newPoly;
	}

	// Differentiate a polynomial with respect to a variable index
	Polynomial differentiate(int ind) {

		// Convert to a string and cache
		std::string indString = std::to_string(ind);
		indString.insert(0, digitsPerInd-indString.size(), ' ');

		// Process each term
		Polynomial result(numVars);
		for (auto const &pair: coeffs) {

			// Count the number of times the variable appears here
			int indFound = -1;
			double degree = 0.0;
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
				if (pair.first.substr(i,digitsPerInd) == indString) {
					degree += 1.0;
					indFound = i;
				}
			}

			// Only do something if there is something with this ind
			if (degree > 0.0) {

				// Remove the latest term of this ind
				std::string newInd = pair.first;
				newInd.erase(newInd.begin()+indFound, newInd.begin()+indFound+digitsPerInd);

				// Add to the result
				if (result.coeffs.find(newInd) != result.coeffs.end()) {
					result.coeffs[newInd] += pair.second * degree;
				} else {
					result.coeffs[newInd] = pair.second * degree;
				}

			}

		}

		// Remove any zero terms and return
		return result.prune();

	}

	// Integrate a polynomial with respect to a variable index
	Polynomial integrate(int ind) {

		// Convert to a string and cache
		std::string indString = std::to_string(ind);
		indString.insert(0, digitsPerInd-indString.size(), ' ');

		// Process each term
		Polynomial result(numVars);
		for (auto const &pair: coeffs) {

			// Count the number of times the variable appears here
			double degree = 0.0;
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
				if (pair.first.substr(i,digitsPerInd) == indString) {
					degree += 1.0;
				}
			}

			// Add an extra term at the right place
			std::string newInd = pair.first;
			bool hasAdded = false;
			for (int i=0; i<newInd.size(); i+=digitsPerInd) {
				if (std::stoi(newInd.substr(i,digitsPerInd)) > ind) {
					newInd.insert(i, indString);
					hasAdded = true;
					break;
				}
			}

			// If the index is larger than everything
			if (!hasAdded) {
				newInd += indString;
			}

			// Add to the result
			if (result.coeffs.find(newInd) != result.coeffs.end()) {
				result.coeffs[newInd] += pair.second / (degree+1);
			} else {
				result.coeffs[newInd] = pair.second / (degree+1);
			}

		}

		// Remove any zero terms and return
		return result.prune();

	}

	// Overload the addition operator
	Polynomial operator+(const Polynomial& other) {

		// Start with one equation
		Polynomial result = other;

		// For each term of the other
		for (auto const &pair: coeffs) {

			// If it's new add it, otherwise combine with the existing
			if (result.coeffs.find(pair.first) != result.coeffs.end()) {
				result.coeffs[pair.first] += pair.second;
				if (std::abs(result.coeffs[pair.first]) < 1e-10) {
					result.coeffs.erase(pair.first);
				}
			} else {
				result.coeffs[pair.first] = pair.second;
			}

		}

		// Remove any zeros and return
		return result;

	}

	// Overload for in-place addition (+=)
	Polynomial& operator+=(const Polynomial& other){

		// For each term of the other
		for (auto const &pair: other.coeffs) {

			// If it's new add it, otherwise combine with the existing
			if (coeffs.find(pair.first) != coeffs.end()) {
				coeffs[pair.first] += pair.second;
				if (std::abs(coeffs[pair.first]) < 1e-10) {
					coeffs.erase(pair.first);
				}
			} else {
				coeffs[pair.first] = pair.second;
			}

		}

		return *this;
	}

	// Overload the multiplication operator
	Polynomial operator*(const Polynomial& other) {

		// For each term of both equations
		Polynomial result(numVars);
		for (auto const &pair1: coeffs) {
			for (auto const &pair2: other.coeffs) {

				// Combine the term list
				std::string combined = "";
				int ind1 = 0;
				int ind2 = 0;
				while (ind1 < pair1.first.size() && ind2 < pair2.first.size()) {
					if (std::stoi(pair1.first.substr(ind1,digitsPerInd)) < std::stoi(pair2.first.substr(ind2,digitsPerInd))) {
						combined += pair1.first.substr(ind1,digitsPerInd);
						ind1 += digitsPerInd;
					} else {
						combined += pair2.first.substr(ind2,digitsPerInd);
						ind2 += digitsPerInd;
					}
				}
				combined += pair1.first.substr(ind1) + pair2.first.substr(ind2);

				// If it's new add it, otherwise combine with the existing
				if (result.coeffs.find(combined) != result.coeffs.end()) {
					result.coeffs[combined] += pair1.second*pair2.second;
				} else {
					result.coeffs[combined] = pair1.second*pair2.second;
				}

			}
		}

		// Remove any zeros and return
		return result.prune();

	}

	// Get the conjugate of the polynomial
	Polynomial conj() {
		Polynomial con(numVars);
		for (auto const &pair: coeffs) {
			con.coeffs[pair.first] = std::conj(pair.second);
		}
		return con;
	}

	// Evaluate a polynomial with x values
	template <typename type>
	std::complex<double> eval(type x) {

		// For each term being added
		std::complex<double> soFar = 0;
		for (auto const &pair: coeffs) {

			// Multiply all the values
			std::complex<double> sub = 1;
			for (int j=0; j<pair.first.size(); j+=digitsPerInd) {
				sub *= x[std::stoi(pair.first.substr(j, digitsPerInd))];
			}

			// Add to the total
			soFar += sub*pair.second;

		}

		return soFar;

	}

	// Prepares this equation for rapid eval by caching things TODO
	void prepareEvalFast() {

		// For each coefficient
		for (auto const &pair: coeffs) {

			// Get the indices as a std::vector
			std::vector<int> indVec;
			for (int j=0; j<pair.first.size(); j+=digitsPerInd) {
				indVec.push_back(std::stoi(pair.first.substr(j, digitsPerInd)));
			}
			inds.push_back(indVec);

			// Convert the coefficient to real
			vals.push_back(std::real(pair.second));

		}

		// Free some space
		coeffs.clear();

	}

	// Evaluate a polynomial with x values
	template <typename type>
	double evalFast(type x) {

		// For each term being added
		double soFar = 0;
		for (int i=0; i<vals.size(); i++) {

			// Multiply all the values
			double sub = 1;
			for (int j=0; j<inds[i].size(); j++) {
				sub *= x[inds[i][j]];
			}

			// Add to the total
			soFar += sub*vals[i];

		}

		return soFar;

	}

	// When doing std::cout << Polynomial
	friend std::ostream &operator<<(std::ostream &output, const Polynomial &other) {

		// For each term
		int numSoFar = 0;
		for (auto const &pair: other.coeffs) {

			// First the coeff
			output << pair.second << "*{";

			// Then the indices
			output << pair.first << "}";

			// Output an addition on everything but the last
			numSoFar += 1;
			if (numSoFar < other.coeffs.size()) {
				output << " + ";
			}

		}

		return output;

	}

	// Checking equality between two polynomials
	bool operator==(const Polynomial &other) {
		for (auto const &pair: other.coeffs) {
			if (coeffs[pair.first] != pair.second) {
				return false;
			}
		}
		return true;
	}

	// Checking inequality between two polynomials
	bool operator!=(const Polynomial &other) {
		for (auto const &pair: other.coeffs) {
			if (coeffs[pair.first] == pair.second) {
				return false;
			}
		}
		return true;
	}

};

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

	// Known ideal
	std::vector<double> idealX(numVars);
	if (d == 2 && n == 2) {
		idealX = {1, 0, 0, 1, rt2, rt2, rt2, -rt2, 0, 0, 0, 0, 0,   0,   0,    0};
	} else if (d == 2 && n == 3) {
		idealX = {1, 0, 0, 1, rt2, rt2, rt2, -rt2,  rt2, 0,   rt2, 0,
									  0, 0, 0, 0, 0,   0,   0,    0,    0,   rt2, 0,   -rt2};
	}

	// The list of equations to fill
	std::vector<Polynomial> eqns;

	// Generate equations
	std::cout << "Generating equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			for (int k=0; k<d; k++) {
				for (int l=0; l<d; l++) {

					// Get the equation before squaring
					Polynomial eqnPreSquare(numVars);
					for (int m=0; m<d; m++) {
						int var1 = i*d*d + k*d + m;
						int var2 = j*d*d + l*d + m;
						int var3 = var1 + conjDelta;
						int var4 = var2 + conjDelta;
						eqnPreSquare.addTerm(1, {var1, var2});
						eqnPreSquare.addTerm(1, {var3, var4});
						eqnPreSquare.addTerm(1i, {var1, var4});
						eqnPreSquare.addTerm(-1i, {var2, var3});
					}

					// For the normalisation equations
					if (i == j && k == l) {
						eqnPreSquare.addTerm(-1.0, {});
					}

					// Square this equation
					Polynomial eqn = eqnPreSquare.conj() * eqnPreSquare;

					// For the MUB-ness equations
					if (i != j) {
						eqn.addTerm(-1.0/d, {});
						eqn = eqn.conj() * eqn;
					}

					// Add it to the list
					eqns.push_back(eqn);

				}
			}
		}
	}

	//// Combine these to create a single polynomial
	std::cout << "Creating single polynomial..." << std::endl;
	Polynomial poly(numVars);
	for (int i=0; i<eqns.size(); i++) {
		poly += eqns[i];
	}

	// Integrate by a variable
	std::cout << "Integrating..." << std::endl;
	Polynomial integrated = poly.integrate(0);

	// Get the gradient TODO openmp
	std::cout << "Differentiating..." << std::endl;
	std::vector<Polynomial> gradient(numVars, Polynomial(numVars));
	for (int i=0; i<numVars; i++) {
		gradient[i] = integrated.differentiate(i);
	}

	// Ensure the first element is the same as the orig
	if (gradient[0] != poly) {
		std::cout << "ERROR - differential of integral is not equivalent!" << std::endl;
		return 0;
	}

	// Get the Hessian TODO openmp
	std::vector<std::vector<Polynomial>> hessian(numVars, std::vector<Polynomial>(numVars, Polynomial(numVars)));
	for (int i=0; i<numVars; i++) {
		for (int j=i; j<numVars; j++) {
			hessian[i][j] = gradient[i].differentiate(j);
		}
	}

	// Pre-cache things to allow for much faster evals
	std::cout << "Optimizing polynomials for fast eval..." << std::endl;
	poly.prepareEvalFast();
	integrated.prepareEvalFast();
	for (int i=0; i<numVars; i++) {
		gradient[i].prepareEvalFast();
		for (int j=i; j<numVars; j++) {
			hessian[i][j].prepareEvalFast();
		}
	}

	// Random starting x
	//std::srand(unsigned(std::time(nullptr)));
	std::srand(0);
	Eigen::VectorXd x = Eigen::VectorXd::Random(numVars);

	// Perform gradient descent using this info
	std::cout << "Minimizing..." << std::endl;
	double alpha = 0.1;
	if (argc > 3) {
		alpha = std::stod(argv[3]);
	}
	Eigen::MatrixXd inv(numVars, numVars);
	Eigen::VectorXd p(numVars);
	Eigen::MatrixXd H(numVars, numVars);
	Eigen::VectorXd g(numVars);
	double maxX = 0;
	for (int iter=0; iter<10000000; iter++) {

		// Calculate the gradient TODO openmp
		for (int i=0; i<numVars; i++) {
			g(i) = gradient[i].evalFast(x);
		}

		// Calculate the Hessian TODO openmp
		for (int i=0; i<numVars; i++) {
			for (int j=i; j<numVars; j++) {
				H(i,j) = hessian[i][j].evalFast(x);
				H(j,i) = H(i,j);
			}
		}

		// Determine the direction
		//p = H.inverse()*g;
		//p = H.partialPivLu().solve(g);
		p = H.householderQr().solve(g);
		
		// Perform a line search
		//alpha = 0.1;
		//double c = 0.0;
		//double tau = 0.9;
		//double m = g.dot(p);
		//double t = -c*m;
		//for (int i=0; i<10; i++) {

			//// Keep going until the Wolfe conditions hold
			//if (std::real(eval(integrated, x)) - std::real(eval(integrated, x+alpha*p)) >= t*alpha) {
				//break;
			//}

			//// Reduce alpha
			//alpha *= tau;

		//}

		// Perform the update
		x -= alpha*p;

		// Per-iteration output
		std::cout << iter << " " << integrated.evalFast(x) << " " << g.norm() << " " << alpha << "          \r" << std::flush;

		// Convergence criteria
		if (std::abs(g(0)) < 1e-10) {
			break;
		}

	}

	// Output the final result
	std::cout << std::endl << std::endl;
	std::cout << "final x = " << std::endl;
	for (int i=0; i<x.size(); i++) {
		std::cout << x(i) << ", ";
	}
	std::cout << std::endl;
	std::cout << "final eval = " << poly.eval(x) << std::endl;

	// Final check for MUB-ness TODO

	return 0;

}
	
