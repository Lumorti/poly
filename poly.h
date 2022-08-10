#ifndef _poly
#define _poly

#include <limits>
#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <unordered_map>
#include <unordered_set>

// Use Eigen for matrix/vector ops
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

// OpenMP for parallelisation
#include <omp.h>

// MOSEK for SOS checking
#include "fusion.h"

// Allow complex literals like 1i
using namespace std::complex_literals;

// Class allowing manipulation of polynomials
template <class polyType>
class Polynomial {

private:

	// Return combinations with repeats
	// choose 2 from {0,1} = 00, 01, 10, 11, 0, 1
	std::vector<std::vector<int>> getAllMonomials(int numVars, int dimension) {

		// Stop when asking for single order monomials
		std::vector<std::vector<int>> toReturn;
		if (dimension == 1) {
			for (int i=0; i<numVars; i++) {
				toReturn.push_back({i});
			}
			return toReturn;
		}

		// For each var, consider this var as the first and then recurse
		for (int i=0; i<numVars; i++) {
			std::vector<std::vector<int>> y = getAllMonomials(numVars, dimension-1);
			for (int j=0; j<y.size(); j++) {
				y[j].insert(y[j].begin(), i);
			}
			toReturn.insert(toReturn.end(), y.begin(), y.end());
		}

		return toReturn;

	}

	// Same as above but converts to Polynomial
	std::vector<Polynomial> getAllMonomialsAsPoly(int numVars, int dimension) {

		// Get the monomials
		std::vector<std::vector<int>> toReturn;
		for (int d=1; d<=dimension; d++) {
			std::vector<std::vector<int>> nextDim = getAllMonomials(numVars, d);
			toReturn.insert(toReturn.end(), nextDim.begin(), nextDim.end());
		}

		// Convert to Polynomial
		std::vector<Polynomial> toReturnPoly(toReturn.size(), Polynomial(numVars));
		for (int i=0; i<toReturn.size(); i++) {
			toReturnPoly[i].addTerm(1, toReturn[i]);
		}

		return toReturnPoly;

	}
	
public:

	// General poly properties
	double zeroTol = 1e-15;
	int maxVariables = 1;
	int digitsPerInd = 1;
	std::unordered_map<std::string, polyType> coeffs;

	// Used for fast eval
	bool fastEvalReady = false;
	std::vector<std::vector<int>> inds;
	std::vector<polyType> vals;

	// Default constructor
	Polynomial(int maxVariables_) {
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));
	}

	// Constructor from another poly
	template <typename type2>
	Polynomial(Polynomial<type2> other) {

		// Copy metadata
		maxVariables = other.maxVariables;
		digitsPerInd = other.digitsPerInd;
		fastEvalReady = other.fastEvalReady;

		// Copy coeffs
		coeffs = std::unordered_map<std::string, polyType>(other.coeffs.size());
		for (auto const &pair: other.coeffs) {
			coeffs[pair.first] = polyType(pair.second);
		}

		// Copy fast eval coeffs
		inds = other.inds;
		vals = std::vector<polyType>(other.vals.size());
		for (int i=0; i<other.vals.size(); i++) {
			vals[i] = polyType(other.vals[i]);
		}

	}

	// Constructor from a complex poly
	template <typename type2>
	Polynomial(Polynomial<std::complex<type2>> other) {

		// Copy metadata
		maxVariables = other.maxVariables;
		digitsPerInd = other.digitsPerInd;
		fastEvalReady = other.fastEvalReady;

		// Copy coeffs
		coeffs = std::unordered_map<std::string, polyType>(other.coeffs.size());
		for (auto const &pair: other.coeffs) {
			coeffs[pair.first] = polyType(std::real(pair.second));
		}

		// Copy fast eval coeffs
		inds = other.inds;
		vals = std::vector<polyType>(other.vals.size());
		for (int i=0; i<other.vals.size(); i++) {
			vals[i] = polyType(std::real(other.vals[i]));
		}

	}

	// Get a list of all the monomials
	std::vector<std::string> getMonomials() {
		std::vector<std::string> monoms;
		for (auto const &pair: coeffs) {
			monoms.push_back(pair.first);
		}
		return monoms;
	}

	// Get a list of all the variable indices
	std::vector<int> getVariables() {
		std::vector<int> listToReturn;

		// For each term
		for (auto const &pair: coeffs) {

			// Check each var
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
				int varInd = std::stoi(pair.first.substr(i, digitsPerInd));

				// Check if this ind is already in the list
				bool found = false;
				for (int j=0; j<listToReturn.size(); j++) {
					if (listToReturn[j] == varInd) {
						found = true;
					}
				}

				// If it's not, add it
				if (!found) {
					listToReturn.push_back(varInd);
				}

			}
		}

		// Sort and return
		return listToReturn;

	}

	// Get the output width of the polynomial
	int getOutputWidth() const {

		// If poly is empty
		if (coeffs.size() == 0) {
			return 4;
		}

		// For each element in this polynomial
		int w = 0;
		for (auto const &pair: coeffs) {

			// The fixed width elements "*", "{" and "}"
			w += 3;

			// The width of the coefficient
			w += 1+std::floor(std::log10(std::abs(pair.second)));

			// An extra if it's negative
			if (pair.second < 0) {
				w += 1;
			}
			
			// The width of the variable list
			w += pair.first.size();

		}

		// We have (size-1) of " + "
		w += (coeffs.size()-1)*3;

		return w;

	}

	// Given a list of input and a list of output, map variables to different indices
	Polynomial changeVariables(std::vector<int> from, std::vector<int> to) {

		// Count the number of unique new indices
		std::unordered_set<int> uniques(to.begin(), to.end());

		// Create a new polynomial with this number of vars
		Polynomial newPoly(uniques.size());

		// Turn the int vecs into string vec
		std::unordered_map<std::string,std::string> mapString;
		for (int i=0; i<from.size(); i++) {
			std::string fromString = std::to_string(from[i]);
			std::string toString = std::to_string(to[i]);
			fromString.insert(0, digitsPerInd-fromString.size(), ' ');
			toString.insert(0, newPoly.digitsPerInd-toString.size(), ' ');
			mapString[fromString] = toString;
		}

		// For each term
		for (auto const &pair: coeffs) {

			std::string newInd = "";
			polyType coeff = pair.second;

			// For each var in the original term
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
				newInd += mapString[pair.first.substr(i, digitsPerInd)];
			}

			// Add this new term if it's non-zero
			if (std::abs(coeff) > zeroTol) {
				newPoly.coeffs[newInd] = coeff;
			}

		}

		return newPoly;

	}

	// Given a list of input and a list of output, map variables to different indices
	Polynomial changeVariables(std::unordered_map<std::string,std::string> mapping) {

		// Create a new polynomial with this number of vars
		Polynomial newPoly(mapping.size());

		// For each term
		for (auto const &pair: coeffs) {

			// Add this new term if it's non-zero
			if (std::abs(pair.second) > zeroTol) {
				newPoly.coeffs[mapping[pair.first]] = pair.second;
			}

		}

		return newPoly;

	}

	// Add a term given a coefficient and an index list
	void addTerm(polyType coeff, std::vector<int> list) {

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
			if (std::abs(coeffs[asString]) < zeroTol) {
				coeffs.erase(asString);
			}
		} else {
			coeffs[asString] = coeff;
		}

	}

	// Add a term given a coefficient and an index string
	void addTermString(polyType coeff, std::string asString) {

		// Add it, combining it if the term already exists
		if (coeffs.find(asString) != coeffs.end()) {
			coeffs[asString] += coeff;
			if (std::abs(coeffs[asString]) < zeroTol) {
				coeffs.erase(asString);
			}
		} else {
			coeffs[asString] = coeff;
		}

	}

	// Return the max power of any monomial
	int getDegree() const {
		int maxDegree = 0;
		for (auto const &pair: coeffs) {
			maxDegree = std::max(maxDegree, int(pair.first.size()) / digitsPerInd);
		}
		return maxDegree;
	}

	// Check if the polynomial can be represented as a sum-of-squares
	bool isSOS() {

		// Get the degree and half it
		int d = getDegree();
		int dHalf = d / 2;

		// Generate the vector of monoms
		std::vector<int> vars = getVariables();
		std::vector<Polynomial> x = getAllMonomialsAsPoly(vars.size(), dHalf);

		// Which matrix elements must sum to each coefficient
		std::unordered_map<std::string,std::vector<std::vector<int>>> cons;
		for (int i=0; i<x.size(); i++) {
			for (int j=0; j<x.size(); j++) {
				cons[(x[i]*x[j]).getMonomials()[0]].push_back({i,j});
			}
		}

		// Create a model
		mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});

		// Create the variable
		mosek::fusion::Variable::t xM = M->variable(mosek::fusion::Domain::inPSDCone(x.size()));

		// For each el + el + el = coeff
		for (auto const &pair: cons) {
			M->constraint(mosek::fusion::Expr::sum(xM->pick(monty::new_array_ptr(pair.second))), mosek::fusion::Domain::equalsTo(coeffs[pair.first]));
		}

		// Solve the problem
		M->solve();

		// If it's infeasible, it's not SOS
		if (M->getProblemStatus() == mosek::fusion::ProblemStatus::PrimalInfeasible) {
			return false;
		} else {
			return true;
		}

	}

	// Substitute a variable for a value
	Polynomial substitute(int ind, polyType toReplace) {

		// Cache the ind to replace as a string
		std::string indString = std::to_string(ind);
		indString.insert(0, digitsPerInd-indString.size(), ' ');

		// For each element in this polynomial
		Polynomial newPoly(maxVariables);
		for (auto const &pair: coeffs) {

			// Remove any instances of this index
			std::string newKey = "";
			polyType newVal = pair.second;
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
				if (pair.first.substr(i, digitsPerInd) == indString) {
					newVal *= toReplace;
				} else {
					newKey += pair.first.substr(i,digitsPerInd); 
				}
			}

			// Then add or create this new term
			if (newPoly.coeffs.find(newKey) != newPoly.coeffs.end()) {
				newPoly.coeffs[newKey] += newVal;
				if (std::abs(newPoly.coeffs[newKey]) < zeroTol) {
					newPoly.coeffs.erase(newKey);
				}
			} else {
				newPoly.coeffs[newKey] = newVal;
			}

		}

		return newPoly;

	}	

	// Substitute several variables for several values 
	Polynomial substitute(std::vector<int> ind, std::vector<polyType> toReplace) {
		if (ind.size() == 0) {
			return *this;
		}
		Polynomial newPoly = substitute(ind[0], toReplace[0]);
		for (int i=1; i<ind.size(); i++) {
			newPoly = newPoly.substitute(ind[i], toReplace[i]);
		}
		return newPoly;
	}

	// Remove zeros from a polynomial
	Polynomial prune() {
		Polynomial newPoly(maxVariables);
		for (auto const &pair: coeffs) {
			if (std::abs(pair.second) > zeroTol) {
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
		Polynomial result(maxVariables);
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

		return result;

	}

	// Integrate a polynomial with respect to a variable index
	Polynomial integrate(int ind) {

		// Convert to a string and cache
		std::string indString = std::to_string(ind);
		indString.insert(0, digitsPerInd-indString.size(), ' ');

		// Process each term
		Polynomial result(maxVariables);
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

		return result;

	}

	// Overload the addition operator with a constant
	template <class otherType>
	Polynomial operator+(const otherType& other) const {

		// Create an equation for the constant
		Polynomial con(maxVariables);
		con.addTerm(polyType(other), {});

		return (*this) + con;

	}

	// Overload the addition operator with a constant
	template <class otherType>
	Polynomial operator-(const otherType& other) const {

		// Create an equation for the constant
		Polynomial con(maxVariables);
		con.addTerm(-polyType(other), {});

		return (*this) + con;

	}

	// Overload the addition operator with another polynomial
	Polynomial operator+(const Polynomial& other) const {

		// Start with one equation
		Polynomial result = other;

		// For each term of the other
		for (auto const &pair: coeffs) {

			// If it's new add it, otherwise combine with the existing
			if (result.coeffs.find(pair.first) != result.coeffs.end()) {
				result.coeffs[pair.first] += pair.second;
				if (std::abs(result.coeffs[pair.first]) < zeroTol) {
					result.coeffs.erase(pair.first);
				}
			} else {
				result.coeffs[pair.first] = pair.second;
			}

		}

		return result;

	}

	// Overload the self-subtraction operator
	Polynomial operator-() const {
		Polynomial negPoly(maxVariables);
		for (auto const &pair: coeffs) {
			negPoly.coeffs[pair.first] = -pair.second;
		}
		return negPoly;
	}

	// Overload the subtraction operator (using the addition)
	Polynomial operator-(const Polynomial& other) const {
		return (*this + (-other));
	}

	// Overload for in-place addition
	Polynomial& operator+=(const Polynomial& other) {
		*this = (*this) + other;
		return *this;
	}

	// Overload for in-place subtraction
	Polynomial& operator-=(const Polynomial& other) {
		*this = (*this) - other;
		return *this;
	}

	// Overload for in-place multiplication 
	Polynomial& operator*=(const Polynomial& other) {
		*this = (*this) * other;
		return *this;
	}

	// Overload the multiplication operator with a constant
	template <class otherType>
	Polynomial operator*(const otherType& other) const {

		// For each term of both equations
		Polynomial result(maxVariables);
		polyType otherConverted = polyType(other);
		for (auto const &pair: coeffs) {

			// Multiply by the constant
			result.coeffs[pair.first] = otherConverted*pair.second;

		}

		return result;

	}
	
	// Overload the multiplication operator with another poly
	Polynomial operator*(const Polynomial& other) const {

		// For each term of both equations
		Polynomial result(std::max(maxVariables, other.maxVariables));
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
				if (ind1 < pair1.first.size()) {
					combined += pair1.first.substr(ind1);
				}
				if (ind2 < pair2.first.size()) {
					combined += pair2.first.substr(ind2);
				}

				// If it's new add it, otherwise combine with the existing
				if (result.coeffs.find(combined) != result.coeffs.end()) {
					result.coeffs[combined] += pair1.second*pair2.second;
				} else {
					result.coeffs[combined] = pair1.second*pair2.second;
				}

			}
		}

		return result;

	}

	// Get the conjugate of the polynomial
	Polynomial conjugate() {
		Polynomial con(maxVariables);
		for (auto const &pair: coeffs) {
			con.coeffs[pair.first] = std::conj(pair.second);
		}
		return con;
	}

	// Get the real part of the polynomial
	Polynomial real() {
		Polynomial con(maxVariables);
		for (auto const &pair: coeffs) {
			con.coeffs[pair.first] = std::real(pair.second);
		}
		return con;
	}

	// Get the imaginary part of the polynomial
	Polynomial imag() {
		Polynomial con(maxVariables);
		for (auto const &pair: coeffs) {
			con.coeffs[pair.first] = std::imag(pair.second);
		}
		return con;
	}

	// Evaluate a polynomial with x values
	polyType eval(std::vector<polyType> x) {

		// For each term being added
		polyType soFar = 0;
		for (auto const &pair: coeffs) {

			// Multiply all the values
			polyType sub = 1;
			for (int j=0; j<pair.first.size(); j+=digitsPerInd) {
				sub *= x[std::stoi(pair.first.substr(j, digitsPerInd))];
			}

			// Add to the total
			soFar += sub*pair.second;

		}

		return soFar;

	}

	// Go from a1*a2 = 0
	//         -> f = (a1+ib1)*(a2+ib2)
	//         -> real(f)**2 + imag(f)**2 = 0
	Polynomial<double> getComplexRelaxation() {

		// First convert this poly to complex
		Polynomial<std::complex<double>> complexPoly(maxVariables*2);

		// For each term in the original
		for (auto const &pair: coeffs) {

			// Start with the coefficient
			Polynomial<std::complex<double>> newTerm(maxVariables*2);
			newTerm.addTerm(pair.second, {});

			// Then for every index
			for (int j=0; j<pair.first.size(); j+=digitsPerInd) {

				// Go from c -> a+ib
				Polynomial<std::complex<double>> thisTerm(maxVariables*2);
				int ind = std::stoi(pair.first.substr(j, digitsPerInd));
				thisTerm.addTerm(1, {ind});
				thisTerm.addTerm(1i, {ind+maxVariables});

				// If there's c1*c2 need (a1+ib1)*(a2+ib2)
				newTerm *= thisTerm;

			}

			// Add this term to the full equation
			complexPoly += newTerm;

		}

		// Square the real and imag parts
		Polynomial<double> newPoly = complexPoly.real()*complexPoly.real() + complexPoly.imag()*complexPoly.imag();
		return newPoly;

	}

	// Prepares this equation for rapid eval by caching things
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

		// Flag that ready for fast eval now
		fastEvalReady = true;

	}

	// Evaluate a polynomial with x values
	template <typename type>
	polyType evalFast(type x) {

		// For each term being added
		polyType soFar = 0;
		for (int i=0; i<vals.size(); i++) {

			// Multiply all the values
			polyType sub = 1;
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

		// If empty
		if (other.coeffs.size() == 0) {
			output << "0*{}";
			return output;
		}

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

	// Return true if the polynomial contains this variable
	bool contains(int varInd) {

		// Cache the ind to find as a string
		std::string indString = std::to_string(varInd);
		indString.insert(0, digitsPerInd-indString.size(), ' ');

		// For each coefficient
		for (auto const &pair: coeffs) {

			// Get the indices as a std::vector
			for (int j=0; j<pair.first.size(); j+=digitsPerInd) {

				// Check if this index is the correct
				if (pair.first.substr(j, digitsPerInd) == indString) {
					return true;
				}

			}

		}

		// Default false
		return false;

	}

	// Try to find a root, with one variable being optimized towards zero
	std::vector<polyType> findRoot(int zeroInd=0, double alpha=0.1, double tolerance=1e-10, int maxIters=10000000) {
		return integrate(zeroInd).findLocalMinimum(alpha, tolerance, maxIters, zeroInd);
	}

	// Size operator (returns number of non-zero monomials)
	long int size() const {
		return coeffs.size();
	}

	// Overload index operator
	std::string& operator[](int index) {
		return coeffs[index];
	}

	// Use the Newton method to find a local minimum
	std::vector<polyType> findLocalMinimum(double alpha=0.1, double tolerance=1e-10, int maxIters=10000000, int zeroInd=0, int threads=4) {

		// Prepare everything for parallel computation
		omp_set_num_threads(threads);
		Eigen::setNbThreads(threads);

		// Get the gradient
		std::vector<Polynomial> gradient(maxVariables, Polynomial(maxVariables));
		#pragma omp parallel for
		for (int i=0; i<maxVariables; i++) {
			gradient[i] = differentiate(i);
		}

		// Get the Hessian
		std::vector<std::vector<Polynomial>> hessian(maxVariables, std::vector<Polynomial>(maxVariables, Polynomial(maxVariables)));
		#pragma omp parallel for
		for (int i=0; i<maxVariables; i++) {
			for (int j=i; j<maxVariables; j++) {
				hessian[i][j] = gradient[i].differentiate(j);
			}
		}

		// Pre-cache things to allow for much faster evals
		#pragma omp parallel for
		for (int i=0; i<maxVariables; i++) {
			gradient[i].prepareEvalFast();
			for (int j=i; j<maxVariables; j++) {
				hessian[i][j].prepareEvalFast();
			}
		}

		// Random starting x
		Eigen::VectorXd x = Eigen::VectorXd::Random(maxVariables);

		// Perform gradient descent using this info
		Eigen::MatrixXd inv(maxVariables, maxVariables);
		Eigen::VectorXd p(maxVariables);
		Eigen::MatrixXd H(maxVariables, maxVariables);
		Eigen::VectorXd g(maxVariables);
		double maxX = 0;
		int iter = 0;
		double norm = 1;
		for (iter=0; iter<maxIters; iter++) {

			// Calculate the gradient
			#pragma omp parallel for
			for (int i=0; i<maxVariables; i++) {
				g(i) = gradient[i].evalFast(x);
			}

			// Calculate the Hessian
			#pragma omp parallel for
			for (int i=0; i<maxVariables; i++) {
				for (int j=i; j<maxVariables; j++) {
					H(i,j) = hessian[i][j].evalFast(x);
					H(j,i) = H(i,j);
				}
			}

			// Determine the direction
			//p = H.householderQr().solve(-g);
			p = H.colPivHouseholderQr().solve(-g);
			//p = (H.transpose() * H).ldlt().solve(H.transpose() * (-g));

			// Perform the update
			x += alpha*p;

			// Per-iteration output
			norm = std::abs(g.norm());
			std::cout << iter << " " << norm << " " << g(zeroInd) << " " << x(zeroInd) << "          \r" << std::flush;

			// Convergence criteria
			if (norm < tolerance) {
				break;
			}

		}
		std::cout << iter << " " << norm << " " << g(zeroInd) << " " << x(zeroInd) << std::endl;

		// Final output
		std::cout.precision(std::numeric_limits<double>::max_digits10);
		std::cout << "Finished in " << iter << " iterations       " << std::endl;
		if (iter == maxIters) {
			std::cout << "WARNING - reached iteration limit" << std::endl;
		}
		std::cout << "Final x = ";
		std::cout << "{ ";
		for (int i=0; i<maxVariables; i++) {
			std::cout << x(i);
			if (i < x.size()-1) {
				std::cout << ", ";
			}
		}
		std::cout << "}" << std::endl;
		std::cout << "Final gradient = " << g(0) << std::endl;

		// Convert the eigen vec into a normal vec
		std::vector<polyType> toReturn(maxVariables);
		for (int i=0; i<maxVariables; i++) {
			toReturn[i] = x(i);
		}
		return toReturn;

	}

};

// For manipulating a whole system of equations at once
template <class polyType>
class PolynomialSystem {
public:

	// Contains the list of polynomials
	std::vector<Polynomial<polyType>> system;

	// Default constructor
	PolynomialSystem() {}

	// Init with a vector of equations
	PolynomialSystem(std::vector<Polynomial<polyType>> system_) {
		system = system_;
	}

	// Add a polynomial
	void addPolynomial(Polynomial<polyType> newPoly) {
		system.push_back(newPoly);
	}

	// Overload index operator
	Polynomial<polyType>& operator[](int index) {
		return system[index];
	}

	// Get a list of all the monomials
	std::vector<std::string> getMonomials() {

		// Use an unordered set to get the unique monomial list
		std::unordered_set<std::string> monoms;
		for (int i=0; i<system.size(); i++) {
			std::vector<std::string> tempList = system[i].getMonomials();
			monoms.insert(tempList.begin(), tempList.end());
		}

		// Turn this into a vector
		std::vector<std::string> monomList;
		for (const std::string& mon: monoms) {
			monomList.push_back(mon);
		}

		return monomList;

	}

	// Return a sorted list of var indices
	std::vector<int> getVariables() {

		// Use an unordered set to get the unique var list
		std::unordered_set<int> varList;
		for (int i=0; i<system.size(); i++) {
			std::vector<int> tempList = system[i].getVariables();
			varList.insert(tempList.begin(), tempList.end());
		}

		// Turn this into a vector
		std::vector<int> varListVec;
		for (const int& var: varList) {
			varListVec.push_back(var);
		}

		// Ensure this keeps the same order
		std::sort(varListVec.begin(), varListVec.end());

		return varListVec;

	}

	// Given a poly system of vars 2,5,6, return poly system of 0,1,2
	PolynomialSystem simplify() {

		// Start with a blank system
		PolynomialSystem newPolySystem;

		// Turn this into a vector
		std::vector<int> varList = getVariables();

		// Map to the simplest indices
		std::vector<int> mapTo;
		for (int i=0; i<varList.size(); i++) {
			mapTo.push_back(i);
		}

		// Replace each var and copy to the new system
		for (int i=0; i<system.size(); i++) {
			newPolySystem.system.push_back(system[i].changeVariables(varList, mapTo));
		}

		return newPolySystem;

	}

	// Remove any equations which are exactly equal
	PolynomialSystem removeDuplicates() {

		// Start with a blank system
		PolynomialSystem newPolySystem;

		// Replace each var and copy to the new system
		for (int i=0; i<system.size(); i++) {

			// See if there exists a copy already
			bool valid = true;
			for (int j=0; j<newPolySystem.system.size(); j++) {
				if (newPolySystem.system[j] == system[i]) {
					valid = false;
					break;
				}
			}

			// If it's not a duplicate, add it
			if (valid) {
				newPolySystem.system.push_back(system[i]);
			}

		}

		return newPolySystem;

	}

	// Given a list of var indices and values, replace everything
	PolynomialSystem substitute(std::vector<int> indsToReplace, std::vector<polyType> valsToReplace) {

		// Start with a blank poly system
		PolynomialSystem newPolySystem;

		// Copy each equation, substituting
		for (int i=0; i<system.size(); i++) {
			newPolySystem.system.push_back(system[i].substitute(indsToReplace, valsToReplace));
		}

		return newPolySystem;

	}

	// Remove zeros and zero equations from the system
	PolynomialSystem prune() {

		// Start with a blank poly system
		PolynomialSystem newPolySystem;

		// Copy each equation if it has at least one term
		for (int i=0; i<system.size(); i++) {
			Polynomial<polyType> pruned = system[i].prune();
			if (pruned.coeffs.size() > 0) {
				newPolySystem.system.push_back(pruned);
			}
		}

		return newPolySystem;

	}

	// Discard equations which contain these variables
	PolynomialSystem withoutVariables(std::vector<int> toRemove) {

		// Start with a blank poly system
		PolynomialSystem newPolySystem;

		// Copy each equation if it doesn't contain this variable
		for (int i=0; i<system.size(); i++) {

			// Get the list of vars for this equation
			std::vector<int> varList = system[i].getVariables();

			// Check if there's any overlap with these two lists
			bool valid = true;
			for (int j=0; j<varList.size(); j++) {
				if (std::find(toRemove.begin(), toRemove.end(), varList[j]) != toRemove.end()) {
					valid = false;
					break;
				}
			}

			// If there isn't, add it
			if (valid) {
				newPolySystem.system.push_back(system[i]);
			}

		}

		return newPolySystem;

	}

	// Only keep equations which contain these variables
	PolynomialSystem withVariables(std::vector<int> toKeep) {

		// Start with a blank poly system
		PolynomialSystem newPolySystem;

		// Copy each equation if it doesn't contain this variable
		for (int i=0; i<system.size(); i++) {

			// Get the list of vars for this equation
			std::vector<int> varList = system[i].getVariables();

			// Check if there's any overlap with these two lists
			bool valid = false;
			for (int j=0; j<varList.size(); j++) {
				if (std::find(toKeep.begin(), toKeep.end(), varList[j]) != toKeep.end()) {
					valid = true;
					break;
				}
			}

			// If there isn't, add it
			if (valid) {
				newPolySystem.system.push_back(system[i]);
			}

		}

		return newPolySystem;

	}

	// Only keep equations which only contain these variables
	PolynomialSystem withOnlyVariables(std::vector<int> toKeep) {

		// Start with a blank poly system
		PolynomialSystem newPolySystem;

		// Copy each equation if it doesn't contain this variable
		for (int i=0; i<system.size(); i++) {

			// Get the list of vars for this equation
			std::vector<int> varList = system[i].getVariables();

			// Check if there's any overlap with these two lists
			bool valid = true;
			for (int j=0; j<varList.size(); j++) {
				if (std::find(toKeep.begin(), toKeep.end(), varList[j]) == toKeep.end()) {
					valid = false;
					break;
				}
			}

			// If there isn't, add it
			if (valid) {
				newPolySystem.system.push_back(system[i]);
			}

		}

		return newPolySystem;

	}
	// Acting like a normal c++ container
	int size() {
		return system.size();
	}

	// When doing stream output
	friend std::ostream &operator<<(std::ostream &output, const PolynomialSystem &other) {

		// For each polynomial
		int numSoFar = 0;
		for (int i=0; i<other.system.size(); i++) {

			// Output this poly plus a space
			output << other.system[i] << std::endl;

		}

		return output;

	}

	// Use Hilbert's Nullstellensatz to attempt to prove infeasibility
	double proveInfeasible(int level=2, int threads=4) {

		// Start with the minimum amount of monomials
		std::vector<int> varList = getVariables();

		// Start with our system of equations
		PolynomialSystem newSys;
		newSys.system = system;

		// For each level, figure out the terms to multiply by
		std::vector<Polynomial<polyType>> toMultiply;
		if (level >= 1) {
			for (int i=0; i<varList.size(); i++) {
				Polynomial<polyType> termToMult(system[0].maxVariables);
				termToMult.addTerm(1, {varList[i]});
				toMultiply.push_back(termToMult);
			}
		}
		if (level >= 2) {
			for (int i=0; i<varList.size(); i++) {
				for (int j=0; j<varList.size(); j++) {
					Polynomial<polyType> termToMult(system[0].maxVariables);
					termToMult.addTerm(1, {varList[i], varList[j]});
					toMultiply.push_back(termToMult);
				}
			}
		}
		if (level >= 3) {
			for (int i=0; i<varList.size(); i++) {
				for (int j=0; j<varList.size(); j++) {
					for (int k=0; k<varList.size(); k++) {
						Polynomial<polyType> termToMult(system[0].maxVariables);
						termToMult.addTerm(1, {varList[i], varList[j], varList[k]});
						toMultiply.push_back(termToMult);
					}
				}
			}
		}
		if (level >= 4) {
			for (int i=0; i<varList.size(); i++) {
				for (int j=0; j<varList.size(); j++) {
					for (int k=0; k<varList.size(); k++) {
						for (int l=0; l<varList.size(); l++) {
							Polynomial<polyType> termToMult(system[0].maxVariables);
							termToMult.addTerm(1, {varList[i], varList[j], varList[k], varList[l]});
							toMultiply.push_back(termToMult);
						}
					}
				}
			}
		}

		// Multiply everything
		for (int i=0; i<toMultiply.size(); i++) {
			for (int j=0; j<system.size(); j++) {
				newSys.addPolynomial(system[j]*toMultiply[i]);
			}
		}

		// Get mapping of monomials
		std::vector<std::string> monomialList = newSys.getMonomials();
		std::unordered_map<std::string,int> monomialMap;
		for (int i=0; i<monomialList.size(); i++) {
			monomialMap[monomialList[i]] = i;
		}

		// Convert everything into a linear system
		std::vector<Eigen::Triplet<double>> tripsletsA;
		for (int i=0; i<newSys.system.size(); i++) {
			for (auto const &pair: newSys.system[i].coeffs) {
				tripsletsA.push_back(Eigen::Triplet<double>(monomialMap[pair.first], i, pair.second));
			}
		}
		Eigen::SparseMatrix<double> A(monomialList.size(), newSys.size());
		A.setFromTriplets(tripsletsA.begin(), tripsletsA.end());

		// We want this to be able to reach 1
		Eigen::SparseVector<double> b(monomialList.size());
		b.insert(monomialMap[""]) = 1.0;

		// Prepare everything for parallel computation
		omp_set_num_threads(threads);
		Eigen::setNbThreads(threads);

		// Try to solve the linear system
		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver(A);
		Eigen::VectorXd cert = solver.solve(b);

		// Calculate the residuals
		double res = (A*cert-b).norm();
		return res;

	}

};


// For manipulating matrices of polynomials
template <class polyType>
class PolynomialMatrix {
public:

	// Contains the list of polynomials
	std::vector<std::vector<Polynomial<polyType>>> system;

	// Matrix info
	int columns = 1;
	int rows = 1;
	int maxVariables = 1;

	// Default constructor
	PolynomialMatrix(int maxVariables_=1, int rows_=1, int columns_=1) {
		rows = rows_;
		columns = columns_;
		maxVariables = maxVariables_;
		system = std::vector<std::vector<Polynomial<polyType>>>(rows, std::vector<Polynomial<polyType>>(columns, Polynomial<polyType>(maxVariables)));
	}

	// Constructor from brace-enclosed list
	PolynomialMatrix(std::vector<double> vec) {
		rows = vec.size();
		columns = 1;
		maxVariables = 1;
		system = std::vector<std::vector<Polynomial<polyType>>>(rows, std::vector<Polynomial<polyType>>(columns, Polynomial<polyType>(maxVariables)));
		for (int i=0; i<rows; i++) {
			system[i][0].addTerm(vec[i], {});
		}
	}

	// Constructor from an Eigen matrix
	PolynomialMatrix(Eigen::MatrixXd other) {
		rows = other.rows();
		columns = other.cols();
		maxVariables = 1;
		system = std::vector<std::vector<Polynomial<polyType>>>(rows, std::vector<Polynomial<polyType>>(columns, Polynomial<polyType>(maxVariables)));
		for (int i=0; i<rows; i++) {
			for (int j=0; j<rows; j++) {
				if (other(i,j) > 0) {
					system[i][j].addTerm(other(i,j), {});
				}
			}
		}
	}

	// Constructor from a PolynomialMatrix
	PolynomialMatrix(const PolynomialMatrix &other) {
		rows = other.rows;
		columns = other.columns;
		maxVariables = other.maxVariables;
		system = std::vector<std::vector<Polynomial<polyType>>>(rows, std::vector<Polynomial<polyType>>(columns, Polynomial<polyType>(maxVariables)));
		for (int i=0; i<rows; i++) {
			for (int j=0; j<columns; j++) {
				system[i][j] = other.system[i][j];
			}
		}
	}

	// Assignment from a PolynomialMatrix
	PolynomialMatrix& operator=(const PolynomialMatrix& other) {
		rows = other.rows;
		columns = other.columns;
		maxVariables = other.maxVariables;
		system = std::vector<std::vector<Polynomial<polyType>>>(rows, std::vector<Polynomial<polyType>>(columns, Polynomial<polyType>(maxVariables)));
		for (int i=0; i<rows; i++) {
			for (int j=0; j<columns; j++) {
				system[i][j] = other.system[i][j];
			}
		}
		return *this;
	}

	// Overload index operator
	std::vector<Polynomial<polyType>>& operator[](int index) {
		return system[index];
	}

	

	// Get the max output width of the polynomial matrix
	std::vector<int> getOutputWidths() const {

		// For each column, figure out the output widths
		std::vector<int> columnWidths(columns, 0);
		for (int j=0; j<columns; j++) {
			for (int i=0; i<rows; i++) {

				// Each column should be the max width of any element
				columnWidths[j] = std::max(columnWidths[j], system[i][j].getOutputWidth()+2);

			}
		}

		return columnWidths;

	}

	// When doing std::cout << PolynomialMatrix
	friend std::ostream &operator<<(std::ostream &output, const PolynomialMatrix &other) {

		// Get the widths of each column
		std::vector<int> columnWidths = other.getOutputWidths();

		// For each row
		for (int i=0; i<other.rows; i++) {

			// Start the row with a symbol
			output << "(  ";

			// Output each column
			for (int j=0; j<other.columns; j++) {
				output << other.system[i][j];

				// Add padding if needed
				int paddingNeeded = columnWidths[j]-other.system[i][j].getOutputWidth();
				output << std::string(paddingNeeded, ' ');

			}

			// Then a newline
			output << ")" << std::endl;

		}

		return output;

	}

	// Overload the multiplication operator with another matrix
	PolynomialMatrix operator*(const PolynomialMatrix& other) {

		// (n x m) * (m x l) = (n x l)
		PolynomialMatrix result(maxVariables, rows, other.columns);

		// For each element of the new matrix
		for (int i=0; i<result.rows; i++) {
			for (int j=0; j<result.columns; j++) {

				// For each element in the multiplication
				for (int k=0; k<columns; k++) {
					result[i][j] += (system[i][k]*other.system[k][j]);
				}

			}
		}

		return result;

	}

	// Overload the addition operator
	PolynomialMatrix operator+(const PolynomialMatrix& other) {

		// Same size as the original
		PolynomialMatrix result(std::max(maxVariables, other.maxVariables), rows, columns);

		// For each element of the new matrix
		for (int i=0; i<rows; i++) {
			for (int j=0; j<columns; j++) {

				// Perform the operation
				result[i][j] = system[i][j] + other.system[i][j];

			}
		}

		return result;

	}

	// Overload the addition operator
	PolynomialMatrix operator-(const PolynomialMatrix& other) {

		// Same size as the original
		PolynomialMatrix result(std::max(maxVariables, other.maxVariables), rows, columns);

		// For each element of the new matrix
		for (int i=0; i<rows; i++) {
			for (int j=0; j<columns; j++) {

				// Perform the operation
				result[i][j] = system[i][j] - other.system[i][j];

			}
		}

		return result;

	}

	// Return the identity as a PolynomialMatrix
	static PolynomialMatrix identity(int matSize) {
		PolynomialMatrix iden(1, matSize, matSize);
		for (int i=0; i<matSize; i++) {
			iden[i][i].addTerm(1, {});
		}
		return iden;
	}

	// Get a list of all the monomials
	std::vector<std::string> getMonomials() {

		// Use an unordered set to get the unique monomial list
		std::unordered_set<std::string> monoms;
		for (int i=0; i<rows; i++) {
			for (int j=0; j<columns; j++) {
				std::vector<std::string> tempList = system[i][j].getMonomials();
				monoms.insert(tempList.begin(), tempList.end());
			}
		}

		// Turn this into a vector
		std::vector<std::string> monomList;
		for (const std::string& mon: monoms) {
			monomList.push_back(mon);
		}

		return monomList;

	}

	// Set each to a new variable
	PolynomialMatrix changeVariables(std::unordered_map<std::string,std::string> mapping) {
		
		// The new matrix may have more variables
		PolynomialMatrix newMat(mapping.size(), rows, columns);

		// For each element
		for (int i=0; i<rows; i++) {
			for (int j=0; j<columns; j++) {

				// Apply the mapping
				newMat[i][j] = system[i][j].changeVariables(mapping);

			}
		}

		return newMat;

	}

};

// Overload the multiplication operator with a polynomial
template <class polyType>
PolynomialMatrix<polyType> operator*(const PolynomialMatrix<polyType>& mat, const Polynomial<polyType>& poly) {

	// Matrix doesn't change shape
	PolynomialMatrix<polyType> result(mat.maxVariables, mat.rows, mat.columns);

	// Needed because otherwise const Poly * const Poly apparently isn't valid
	Polynomial<polyType> polyCopy = poly;

	// For each element
	for (int i=0; i<mat.rows; i++) {
		for (int j=0; j<mat.columns; j++) {
			result[i][j] = mat[i][j] * polyCopy;
		}
	}

	return result;

}

// Overload the multiplication operator with a polynomial
template <class polyType>
PolynomialMatrix<polyType> operator*(const Polynomial<polyType>& poly, const PolynomialMatrix<polyType>& mat) {

	// Matrix doesn't change shape
	PolynomialMatrix<polyType> result(mat.maxVariables, mat.rows, mat.columns);

	// Needed because otherwise const Poly * const Poly apparently isn't valid
	Polynomial<polyType> polyCopy = poly;

	// For each element
	for (int i=0; i<mat.rows; i++) {
		for (int j=0; j<mat.columns; j++) {
			result[i][j] = polyCopy*mat.system[i][j];
		}
	}

	return result;

}

// Switch the order
template <class polyType, class otherType>
Polynomial<polyType> operator*(const otherType& other, const Polynomial<polyType>& poly) {
	return poly*other;
}

// Switch the order
template <class polyType, class otherType>
Polynomial<polyType> operator+(const otherType& other, const Polynomial<polyType>& poly) {
	return poly+other;
}

// Switch the order
template <class polyType, class otherType>
Polynomial<polyType> operator-(const otherType& other, const Polynomial<polyType>& poly) {
	return poly-other;
}

// Generic overload for outputting vector
template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "{";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii) {
        os << *ii;
		if (ii != v.end()) {
			os << ", ";
		}
    }
    os << "}";
    return os;
}

#endif
	
