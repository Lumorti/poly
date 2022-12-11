#ifndef _poly
#define _poly

#include <limits>
#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <random>
#include <algorithm>
#include <iomanip>

// Use Eigen for matrix/vector ops
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

// OpenMP for parallelisation
#include <omp.h>

// MOSEK for SOS checking
#include "fusion.h"

// Allow complex literals like 1i
using namespace std::complex_literals;

// Generic overload for outputting vector of strings
std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& v) {
    os << "{";
    for (typename std::vector<std::string>::const_iterator ii = v.begin(); ii != v.end(); ++ii) {
        os << "\"" << *ii << "\"";
		if (ii+1 != v.end()) {
			os << ", ";
		}
    }
    os << "}";
    return os;
}

// Generic overload for outputting a pair
template <class type1, class type2>
std::ostream& operator<<(std::ostream& os, const std::pair<type1,type2>& v) {
	os << "{" << v.first << ", " << v.second << "}";
    return os;
}

// Generic overload for outputting a tuple
template <class type1, class type2, class type3>
std::ostream& operator<<(std::ostream& os, const std::tuple<type1,type2,type3>& v) {
	os << "{" << std::get<0>(v) << ", " << std::get<1>(v) << ", " << std::get<2>(v) << "}";
    return os;
}

// Generic overload for outputting map
template <class T, class T2>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<T,T2>& v) {
    os << "{";
	int vSize = v.size();
	int i = 0;
    for (typename std::unordered_map<T,T2>::const_iterator ii = v.begin(); ii != v.end(); ++ii) {
        os << (*ii).first << ": " << (*ii).second;
		if (i < vSize) {
			os << ", ";
		}
		i++;
    }
    os << "}";
    return os;
}

// Generic overload for outputting vector
template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "{";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii) {
        os << *ii;
		if (ii+1 != v.end()) {
			os << ", ";
		}
    }
    os << "}";
    return os;
}

// Generic overload for outputting vector of vector
template <typename type> 
std::ostream &operator<<(std::ostream &output, const std::vector<std::vector<type>> &arr) {

	// Used fixed precision
	int precision = 7;
	output << std::fixed << std::setprecision(precision);

	// Loop over the array
	std::string rowText;
	for (int y=0; y<arr.size(); y++) {

		// For the first line, add the pre text
		if (y == 0) {
			rowText = " {";

		// Otherwise pad accordingly
		} else {
			rowText = "  ";
		}

		// Spacing
		output << rowText << " { ";

		// For the x values, combine them all on one line
		for (int x=0; x<arr[y].size(); x++) {
			output << std::setw(precision+1) << arr[y][x];
			if (x < arr[y].size()-1) {
				output << ", ";
			}
		}

		// Output the row
		if (y < arr.size() - 1) {
			output << "}, " << std::endl;
		}

	}

	// Output the final closing braces
	output << "} } ";
	
	return output;

}

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
	bool isNaN = false;
	std::unordered_map<std::string, polyType> coeffs;

	// Used for fast eval
	bool fastEvalReady = false;
	std::vector<std::vector<int>> inds;
	std::vector<polyType> vals;

	// Default constructor
	Polynomial() {
		isNaN = true;
	}

	// Simplest constructor
	Polynomial(int maxVariables_) {
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));
	}

	// Constructor with a single term
	Polynomial(int maxVariables_, polyType coeff, std::vector<int> inds) {
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));
		addTerm(coeff, inds);
	}

	// Constructor with a single term
	Polynomial(int maxVariables_, polyType coeff, std::initializer_list<int> inds) {
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));
		addTerm(coeff, inds);
	}

	// Constructor with a single constant term
	Polynomial(int maxVariables_, polyType coeff) {
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));
		addTerm(coeff);
	}

	// Constructor with a single term from string index
	Polynomial(int maxVariables_, polyType coeff, std::string inds) {
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));
		addTerm(coeff, inds);
	}

	// Constructor from the string form
	Polynomial(int maxVariables_, std::string asString) {
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));

		// Loop over the string form
		std::string currentCoeff = "";
		std::string currentInd = "";
		int lookingFor = 0; // 0 = coeff, 1 = ind
		for (int i=0; i<asString.size(); i++) {

			// When we find this, we know the next thing will be an index list
			if (asString[i] == '*') {
				lookingFor = 1;
				i++;

			// When we find this, we have a coeff + ind
			} else if (asString[i] == '}') {
				addTerm(std::stod(currentCoeff), currentInd);
				currentInd = "";
				currentCoeff = "";
				lookingFor = 0;

			// Otherwise add to the correct variable
			} else if (lookingFor == 0) {
				currentCoeff += asString[i];
			} else if (lookingFor == 1) {
				currentInd += asString[i];
			}

		}


	}


	// Change the digits per index
	Polynomial changeMaxVariables(int newMaxVariables) {

		// If no change, do nothing
		if (newMaxVariables == maxVariables) {
			return *this;
		}

		// For each coefficient
		Polynomial newPoly(newMaxVariables);
		for (auto const &pair: coeffs) {

			// Convert to the new size
			std::string newIndString = "";
			for (int j=0; j<pair.first.size(); j+=digitsPerInd) {
				int ind = std::stoi(pair.first.substr(j, digitsPerInd));
				std::string padded = std::to_string(ind);
				padded.insert(0, newPoly.digitsPerInd-padded.size(), ' ');
				newIndString += padded;
			}

			// Update this coeff
			newPoly.coeffs[newIndString] = pair.second;

		}

		return newPoly;

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

	// Get a list of all the monomials
	std::vector<std::string> getMonomials() const {
		std::vector<std::string> monoms;
		for (auto const &pair: coeffs) {
			monoms.push_back(pair.first);
		}
		return monoms;
	}

	// Get a list of all the monomials as Polynomials
	std::vector<Polynomial<polyType>> getMonomialsAsPolynomials() const {
		std::vector<Polynomial<polyType>> monoms;
		for (auto const &pair: coeffs) {
			monoms.push_back(Polynomial<polyType>(maxVariables, 1, pair.first));
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

		// Just return (no sorting)
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

	// Add a term given just a coefficient
	void addTerm(polyType coeff) {

		// Convert the nothing to a string
		std::string asString = "";

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

	// Add a term given a coefficient and an index list
	void addTerm(polyType coeff, std::vector<int> list) {

		// Sort the list
		std::sort(list.begin(), list.end());

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
	
	// Add a term given a coefficient and an index list
	void addTerm(polyType coeff, std::initializer_list<int> initList) {

		// Sort the list
		std::vector<int> list(initList.size());
		std::copy(initList.begin(), initList.end(), list.begin());
		std::sort(list.begin(), list.end());

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
	void addTerm(polyType coeff, std::string asString) {

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

	// Return the max power of any variable in any monomial
	int getDegreePerVariable() const {
		int maxDegree = 0;
		for (auto const &pair: coeffs) {
			int currentDegree = 0;
			std::string currentVar = "";
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
				if (pair.first.substr(i, digitsPerInd) == currentVar) {
					currentDegree += 1;
				} else {
					currentDegree = 1;
					currentVar = pair.first.substr(i, digitsPerInd);
				}
				maxDegree = std::max(maxDegree, currentDegree);
			}
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

	// Remove doubles from indices, such that given 0012 -> 12, 233 -> 2
	Polynomial removeDuplicates() {

		// For each element in this polynomial
		Polynomial newPoly(maxVariables);
		for (auto const &pair: coeffs) {

			// Remove any instances of this index
			std::string newKey = "";
			polyType newVal = pair.second;
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
				if (i < pair.first.size() - digitsPerInd && pair.first.substr(i, digitsPerInd) == pair.first.substr(i+digitsPerInd, digitsPerInd)) {
					i += digitsPerInd;
				} else {
					newKey += pair.first.substr(i, digitsPerInd); 
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

	// Apply an int index to int index mapping
	Polynomial replaceWithVariable(std::unordered_map<int, int> mapping) {

		// Create a new polynomial with this number of vars
		Polynomial newPoly(maxVariables);

		// Turn the int vecs into string vec
		std::unordered_map<std::string, std::string> mapString;
		for (auto const &pair: mapping) {
			std::string fromString = std::to_string(pair.first);
			std::string toString = std::to_string(pair.second);
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
				if (mapString.find(pair.first.substr(i, digitsPerInd)) != mapString.end()) {
					newInd += mapString[pair.first.substr(i, digitsPerInd)];
				} else {
					newInd += pair.first.substr(i, digitsPerInd);
				}
			}

			// Add this new term if it's non-zero
			if (std::abs(coeff) > zeroTol) {

				// Add it, combining it if the term already exists
				if (newPoly.coeffs.find(newInd) != newPoly.coeffs.end()) {
					newPoly.coeffs[newInd] += coeff;
					if (std::abs(newPoly.coeffs[newInd]) < zeroTol) {
						newPoly.coeffs.erase(newInd);
					}
				} else {
					newPoly.coeffs[newInd] = coeff;
				}

			}

		}

		return newPoly;

	}

	// Apply a string index to string index mapping
	Polynomial replaceWithVariable(std::unordered_map<std::string, std::string> mapping) {

		// Create a new polynomial with this number of vars
		Polynomial newPoly(mapping.size());

		// For each term
		for (auto const &pair: coeffs) {

			// Add this new term if it's non-zero
			if (std::abs(pair.second) > zeroTol) {

				// Apply the mapping to get the new index
				std::string newInd = mapping[pair.first];
				polyType coeff = pair.second;

				// Add it, combining it if the term already exists
				if (newPoly.coeffs.find(newInd) != newPoly.coeffs.end()) {
					newPoly.coeffs[newInd] += coeff;
					if (std::abs(newPoly.coeffs[newInd]) < zeroTol) {
						newPoly.coeffs.erase(newInd);
					}
				} else {
					newPoly.coeffs[newInd] = coeff;
				}

			}

		}

		return newPoly;

	}

	// Substitute a variable for a polynomial
	Polynomial replaceWithPoly(std::string indString, Polynomial<polyType> toReplace) {

		// For each element in this polynomial
		Polynomial newPoly(toReplace.maxVariables);
		for (auto const &pair: coeffs) {

			// Find any instance of this index
			Polynomial toMultiply2(toReplace.maxVariables, 1, {});
			std::string remainingIndex = "";
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
				if (pair.first.substr(i, digitsPerInd) == indString) {
					toMultiply2 *= toReplace;
				} else {
					remainingIndex += pair.first.substr(i, digitsPerInd);
				}
			}

			// Get this monomial as a seperate poly
			Polynomial toMultiply1(toReplace.maxVariables);
			toMultiply1.addTerm(pair.second, remainingIndex);

			// Add to the main poly
			newPoly += toMultiply1*toMultiply2;

		}

		return newPoly;

	}

	// Substitute a variable for a polynomial
	Polynomial replaceWithPoly(int ind, Polynomial<polyType> toReplace) {

		// Cache the ind to replace as a string
		std::string indString = std::to_string(ind);
		indString.insert(0, digitsPerInd-indString.size(), ' ');

		// For each element in this polynomial
		Polynomial newPoly(toReplace.maxVariables);
		for (auto const &pair: coeffs) {

			// Find any instance of this index
			Polynomial toMultiply2(toReplace.maxVariables, 1, {});
			std::string remainingIndex = "";
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
				if (pair.first.substr(i, digitsPerInd) == indString) {
					toMultiply2 *= toReplace;
				} else {
					remainingIndex += pair.first.substr(i, digitsPerInd);
				}
			}

			// Get this monomial as a seperate poly
			Polynomial toMultiply1(toReplace.maxVariables);
			toMultiply1.addTerm(pair.second, remainingIndex);

			// Add to the main poly
			newPoly += toMultiply1*toMultiply2;

		}

		return newPoly;

	}

	// Substitute several variables for several polynomials
	Polynomial replaceWithPoly(std::vector<int> inds, std::vector<Polynomial<polyType>> toReplace) {

		// Cache the inds to replace as a string
		std::unordered_map<std::string,int> stringToLoc;
		for (int i=0; i<inds.size(); i++) {
			std::string indString = std::to_string(inds[i]);
			indString.insert(0, digitsPerInd-indString.size(), ' ');
			stringToLoc[indString] = i;
		}

		// For each element in this polynomial
		Polynomial newPoly(toReplace[0].maxVariables);
		for (auto const &pair: coeffs) {

			// Find any instance of this index
			Polynomial toMultiply2(toReplace[0].maxVariables, 1, {});
			std::string remainingIndex = "";
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
				if (stringToLoc.find(pair.first.substr(i, digitsPerInd)) != stringToLoc.end()) {
					toMultiply2 *= toReplace[stringToLoc[pair.first.substr(i, digitsPerInd)]];
				} else {
					remainingIndex += pair.first.substr(i, digitsPerInd);
				}
			}
			
			// Get whatever's left as a seperate poly
			Polynomial toMultiply1(toReplace[0].maxVariables);
			toMultiply1.addTerm(pair.second, remainingIndex);
			toMultiply2 *= toMultiply1;

			// Add to the main poly
			newPoly += toMultiply2;

		}

		return newPoly;

	}

	// Substitute a variable for a value
	Polynomial replaceWithValue(int ind, polyType toReplace) {

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

	// Substitute a string variable for a value 
	Polynomial replaceWithValue(std::string indString, polyType toReplace) {

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
	Polynomial replaceWithValue(std::vector<int> ind, std::vector<polyType> toReplace) {
		if (ind.size() == 0) {
			return *this;
		}
		Polynomial newPoly = replaceWithValue(ind[0], toReplace[0]);
		for (int i=1; i<ind.size(); i++) {
			newPoly = newPoly.replaceWithValue(ind[i], toReplace[i]);
		}
		return newPoly;
	}

	// Substitute several variables for several values
	Polynomial replaceWithValue(std::vector<polyType> toReplace) {
		if (toReplace.size() == 0) {
			return *this;
		}
		Polynomial newPoly = replaceWithValue(0, toReplace[0]);
		for (int i=1; i<toReplace.size(); i++) {
			newPoly = newPoly.replaceWithValue(i, toReplace[i]);
		}
		return newPoly;
	}

	// Substitute several variables for several values (with strings)
	Polynomial replaceWithValue(std::vector<std::string> ind, std::vector<polyType> toReplace) {
		if (ind.size() == 0) {
			return *this;
		}
		Polynomial newPoly = replaceWithValue(ind[0], toReplace[0]);
		for (int i=1; i<ind.size(); i++) {
			newPoly = newPoly.replaceWithValue(ind[i], toReplace[i]);
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
		Polynomial con(maxVariables, polyType(other), {});
		return (*this) + con;
	}

	// Overload the addition operator with a constant
	template <class otherType>
	Polynomial operator-(const otherType& other) const {
		Polynomial con(maxVariables, -polyType(other), {});
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
	template <class otherType>
	Polynomial& operator+=(const otherType& other) {
		*this = (*this) + other;
		return *this;
	}

	// Overload for in-place subtraction
	template <class otherType>
	Polynomial& operator-=(const otherType& other) {
		*this = (*this) - other;
		return *this;
	}

	// Overload for in-place multiplication
	template <class otherType>
	Polynomial& operator*=(const otherType& other) {
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

	// Overload the division operator with another poly
	Polynomial operator/(const Polynomial& other) const {

		// If dividing by a more complex poly, return a null poly
		if (other.size() != 1 || isNaN || other.isNaN) {
			return Polynomial();
		}

		// Copy the values from the other poly
		std::string otherMonom;
		polyType otherVal;
		for (auto const &pair: other.coeffs) {
			otherMonom = pair.first;
			otherVal = pair.second;
		}

		// For each term of this equation
		Polynomial result(std::max(maxVariables, other.maxVariables));
		for (auto const &pair: coeffs) {

			// For each variable in the dividing term
			std::string thisMonom = pair.first;
			for (int i=0; i<otherMonom.size(); i+=digitsPerInd) {

				// If there's in the numerator left, we have null
				if (thisMonom.size() <= 0) {
					return Polynomial();
				}

				// Otherwise try to remove the corresponding term
				bool removed = false;
				for (int j=0; j<thisMonom.size(); j+=digitsPerInd) {
					if (otherMonom.substr(i, digitsPerInd) == thisMonom.substr(j, digitsPerInd)) {
						thisMonom.erase(j, digitsPerInd);
						removed = true;
						break;
					}
				}

				// If left with an resolved division, it's NaN
				if (!removed) {
					return Polynomial();
				}

			}

			// If it's new add it, otherwise combine with the existing
			if (result.coeffs.find(thisMonom) != result.coeffs.end()) {
				result.coeffs[thisMonom] += pair.second/otherVal;
			} else {
				result.coeffs[thisMonom] = pair.second/otherVal;
			}

		}

		return result;

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

	// Prepares this equation for rapid eval by caching things
	void prepareEvalMixed() {

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

		// Flag that ready for fast eval now
		fastEvalReady = true;

	}

	// Prepares this equation for ONLY fast eval (saves memory)
	void prepareEvalFast() {
		
		// Prepare for fast eval
		prepareEvalMixed();

		// Free some space
		coeffs.clear();

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

	// Output in a form suitable for use elsewhere
	std::string asMathematica() const {

		// For each term
		std::string toReturn = "";
		int numSoFar = 0;
		for (auto const &pair: coeffs) {

			// If it contains at least one variable
			if (pair.first != "") {

				// First the coeff
				if (pair.second != 1 && std::abs(pair.second) > zeroTol) {
					toReturn += std::to_string(pair.second);
					toReturn += "*";
				}

				// Then the indices
				for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
					toReturn += "x[";
					toReturn += pair.first.substr(i, digitsPerInd);
					toReturn += "]";
					if (i < pair.first.size()-1) {
						toReturn += "*";
					}
				}

			// For the constant, don't need any variables
			} else {
				toReturn += std::to_string(pair.second);
			}

			// Output an addition on everything but the last
			numSoFar += 1;
			if (numSoFar < coeffs.size()) {
				toReturn += " + ";
			}

		}

		return toReturn;

	}
	
	// Output in a form suitable for use elsewhere
	std::string asLaTeX() const {

		// For each term
		std::string toReturn = "";
		int numSoFar = 0;
		for (auto const &pair: coeffs) {

			// If it contains at least one variable
			if (pair.first != "") {

				// First the coeff
				if (pair.second != 1 && std::abs(pair.second) > zeroTol) {
					toReturn += std::to_string(pair.second);
				}

				// Then the indices
				for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
					toReturn += "x_";
					toReturn += pair.first.substr(i, digitsPerInd);
				}

			// For the constant, don't need any variables
			} else {
				toReturn += std::to_string(pair.second);
			}

			// Output an addition on everything but the last
			numSoFar += 1;
			if (numSoFar < coeffs.size()) {
				toReturn += " + ";
			}

		}

		return toReturn;

	}

	// Output in a form suitable for use elsewhere
	std::string asString() const {

		// If null poly
		if (isNaN) {
			return "NaN";
		}

		// If empty
		if (coeffs.size() == 0) {
			return "0*{}";
		}

		// For each term
		std::string toReturn = "";
		int numSoFar = 0;
		for (auto const &pair: coeffs) {

			// First the coeff
			toReturn += std::to_string(pair.second) + "*{";

			// Then the indices
			toReturn += pair.first + "}";

			// Output an addition on everything but the last
			numSoFar += 1;
			if (numSoFar < coeffs.size()) {
				toReturn += " + ";
			}

		}

		return toReturn;

	}

	// When doing std::cout << Polynomial
	friend std::ostream &operator<<(std::ostream &output, const Polynomial &other) {
		output << other.asString();
		return output;
	}

	// Assigning with a constant
	Polynomial& operator=(const polyType& other) {
		coeffs.clear();
		coeffs[""] = other;
		return *this;
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

	// Return true if the polynomial contains this variable
	bool contains(std::string indString) {

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

	// Size operator (returns number of non-zero monomials)
	long int size() const {
		return coeffs.size();
	}

	// Overload index operator with int index
	polyType& operator[](std::vector<int> indexAsInts) {
		std::string index = "";
		for (int i=0; i<indexAsInts.size(); i++) {
			std::string indString = std::to_string(indexAsInts[i]);
			indString.insert(0, digitsPerInd-indString.size(), ' ');
			index += indString;
		}
		return coeffs[index];
	}

	// Overload index operator with int index
	polyType& operator[](std::initializer_list<int> initList) {
		std::vector<int> indexAsInts(initList.size());
		std::copy(initList.begin(), initList.end(), indexAsInts);
		std::string index = "";
		for (int i=0; i<indexAsInts.size(); i++) {
			std::string indString = std::to_string(indexAsInts[i]);
			indString.insert(0, digitsPerInd-indString.size(), ' ');
			index += indString;
		}
		return coeffs[index];
	}

	// Overload index operator with string index
	polyType& operator[](std::string index) {
		return coeffs[index];
	}

	// Try to find a root, with one variable being optimized towards zero
	std::vector<polyType> findRoot(int zeroInd=0, double alpha=0.1, double tolerance=1e-10, int maxIters=10000000, int threads=4, bool verbose=false) {
		return integrate(zeroInd).findLocalMinimum(zeroInd, alpha, tolerance, maxIters, threads, verbose);
	}

	// Use the Newton method to find a local minimum
	std::vector<polyType> findLocalMinimum(int zeroInd=0, double alpha=0.1, double tolerance=1e-10, int maxIters=10000000, int threads=4, bool verbose=false) {

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
		double minVal = 1e10;
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
			p = H.colPivHouseholderQr().solve(-g);

			// If this is zero, follow the gradient
			if (p.norm() < 1e-10) {
				p = 100*Eigen::VectorXd::Random(maxVariables);
			}

			// Perform the update
			x += alpha*p;

			// Per-iteration output
			norm = std::abs(g.norm());
			minVal = std::min(norm, minVal);
			std::cout << iter << " " << norm << " " << g(zeroInd) << " " << x(zeroInd) << "   " << minVal << "  "  << p.norm() << "          \r" << std::flush;

			// Convergence criteria
			if (norm < tolerance) {
				break;
			}

		}
		std::cout << iter << " " << norm << " " << g(zeroInd) << " " << x(zeroInd) << std::endl;

		// Final output
		if (verbose) {
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
		}

		// Convert the eigen vec into a normal vec
		std::vector<polyType> toReturn(maxVariables);
		for (int i=0; i<maxVariables; i++) {
			toReturn[i] = x(i);
		}
		return toReturn;

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
	PolynomialMatrix(std::vector<polyType> vec) {
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
	PolynomialMatrix replaceWithVariable(std::unordered_map<std::string,std::string> mapping) {
		
		// The new matrix may have more variables
		PolynomialMatrix newMat(mapping.size(), rows, columns);

		// For each element
		for (int i=0; i<rows; i++) {
			for (int j=0; j<columns; j++) {

				// Apply the mapping
				newMat[i][j] = system[i][j].replaceWithVariable(mapping);

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
template <typename polyType, typename otherType>
Polynomial<polyType> operator*(const otherType& other, const Polynomial<polyType>& poly) {
	return poly*other;
}

// Switch the order
template <typename polyType, typename otherType>
Polynomial<polyType> operator+(const otherType& other, const Polynomial<polyType>& poly) {
	return poly+other;
}

// Switch the order
template <typename polyType, typename otherType>
Polynomial<polyType> operator-(const otherType& other, const Polynomial<polyType>& poly) {
	return poly-other;
}

// Overload std::real
template <typename polyType2, typename polyType>
Polynomial<polyType2> std::real(const Polynomial<polyType>& poly) {
	Polynomial<polyType2> newPoly(poly.maxVariables);
	for (auto const &pair: poly.coeffs) {
		newPoly.coeffs[pair.first] = std::real(pair.second);
	}
	return newPoly;
}

// Overload std::imag
template <typename polyType2, typename polyType>
Polynomial<polyType2> std::imag(const Polynomial<polyType>& poly) {
	Polynomial<polyType2> newPoly(poly.maxVariables);
	for (auto const &pair: poly.coeffs) {
		newPoly.coeffs[pair.first] = std::imag(pair.second);
	}
	return newPoly;
}

// Overload the std::conj
template <typename polyType>
Polynomial<polyType> std::conj(const Polynomial<polyType>& poly) {
	Polynomial<polyType> newPoly(poly.maxVariables);
	for (auto const &pair: poly.coeffs) {
		newPoly.coeffs[pair.first] = std::conj(pair.second);
	}
	return newPoly;
}

// For minimizing a polynomial of binary vars subject to constraints
template <typename polyType>
class PolynomialBinaryProblem {
public:

	// The things definiting the problem
	int maxVariables = 1;
	int digitsPerInd = 1;
	Polynomial<polyType> obj;
	std::vector<Polynomial<polyType>> conZero;
	std::vector<Polynomial<polyType>> conPositive;

	// Constructor with everything
	PolynomialBinaryProblem(Polynomial<polyType> obj_, std::vector<Polynomial<polyType>> conZero_, std::vector<Polynomial<polyType>> conPositive_) {
		obj = obj_;
		conZero = conZero_;
		conPositive = conPositive_;
		maxVariables = obj.maxVariables;
		digitsPerInd = obj.digitsPerInd;
	}

	// Get a list of all the monomials
	std::vector<std::string> getMonomials() {

		// Use an unordered set to get the unique monomial list
		std::unordered_set<std::string> monoms;
		std::vector<std::string> tempList = obj.getMonomials();
		monoms.insert(tempList.begin(), tempList.end());
		for (int i=0; i<conZero.size(); i++) {
			tempList = conZero[i].getMonomials();
			monoms.insert(tempList.begin(), tempList.end());
		}
		for (int i=0; i<conPositive.size(); i++) {
			tempList = conPositive[i].getMonomials();
			monoms.insert(tempList.begin(), tempList.end());
		}

		// Turn this into a vector
		std::vector<std::string> monomList;
		for (const std::string& mon: monoms) {
			monomList.push_back(mon);
		}

		return monomList;

	}

	// Given a list of var indices and values, replace everything
	template <typename otherType>
	PolynomialBinaryProblem replaceWithValue(std::vector<int> indsToReplace, std::vector<otherType> valsToReplace) {

		// Convert the list to the correct type
		std::vector<polyType> convertedList(valsToReplace.size());
		for (int i=0; i<valsToReplace.size(); i++) {
			convertedList[i] = polyType(valsToReplace[i]);
		}

		// Start with a blank poly system
		PolynomialBinaryProblem newPolyProblem({}, {}, {});

		// Copy each equation, substituting
		newPolyProblem.obj = obj.replaceWithValue(indsToReplace, convertedList);
		for (int i=0; i<conZero.size(); i++) {
			newPolyProblem.conZero.push_back(conZero[i].replaceWithValue(indsToReplace, convertedList));
		}
		for (int i=0; i<conPositive.size(); i++) {
			newPolyProblem.conPositive.push_back(conPositive[i].replaceWithValue(indsToReplace, convertedList));
		}

		return newPolyProblem;

	}

	// When doing std::cout << PolynomialBinaryProblem
	friend std::ostream &operator<<(std::ostream &output, const PolynomialBinaryProblem &other) {

		// Output the objective
		output << "Minimize: " << std::endl << std::endl;;
		output << other.obj << std::endl << std::endl;

		// Output each constraint
		int numSoFar = 0;
		if (other.conZero.size() + other.conPositive.size() > 0) {
			output << "Subject to: " << std::endl;
			for (int i=0; i<other.conZero.size(); i++) {
				output << std::endl << other.conZero[i] << " = 0 " << std::endl;;
			}
			for (int i=0; i<other.conPositive.size(); i++) {
				output << std::endl << other.conPositive[i] << " > 0 " << std::endl;
			}
		}

		return output;

	}

	// Returns true if all of the constraints are satisfied
	bool isSatisfied() {

		// Check each constraint
		for (int i=0; i<conZero.size(); i++) {
			if (abs(conZero[i][""]) > conZero[i].zeroTol) {
				return false;
			}
		}
		for (int i=0; i<conPositive.size(); i++) {
			if (conZero[i][""] > conZero[i].zeroTol) {
				return false;
			}
		}

		// If we got to the end, it's satisfied
		return true;

	}

	// Find the exact solution through brute force
	std::pair<polyType,std::vector<int>> bruteForce() {

		// Create a vector listing all the variables
		int numVars = obj.getVariables().size();
		std::vector<int> inds(numVars);
		for (int i=0; i<numVars; i++) {
			inds[i] = i;
		}

		// For each possible set of variables
		int numCombs = std::pow(2, obj.maxVariables);
		polyType bestVal = 10000000;
		std::vector<int> bestSol(numVars);
		for (int k=0; k<numCombs; k++) {

			// Convert to -1/1
			std::vector<int> sol(numVars, -1);
			for (int i=0; i<numVars; i++) {
				if ((k >> i) & 1) {
					sol[numVars-i-1] = 1;
				}
			}

			// Substitute all the values
			PolynomialBinaryProblem testing = replaceWithValue(inds, sol);

			// If it's a valid interior point
			if (testing.isSatisfied()) {

				// See if this objective is better
				if (testing.obj[""] < bestVal) {
					bestVal = testing.obj[""];
					bestSol = sol;
				}

			}

		}

		return {bestVal, bestSol};
	
	}

	// Solve a linear program given an objective and zero/positive constraints
	std::pair<polyType,std::vector<polyType>> solveLinear(Polynomial<polyType>& objLinear, std::vector<Polynomial<polyType>>& conZeroLinear, std::vector<Polynomial<polyType>> conPositiveLinear, std::vector<std::string>& monoms, std::vector<std::vector<int>>& monomPairs) {

		// Set some vars
		int oneIndex = 0;
		int varsTotal = monoms.size();

		// Create the PSD matrices from this list
		for (int j=0; j<monomPairs.size(); j++) {

			// If it's a second-order
			if (monomPairs[j].size() == 3) {

				// Extract the vals and construct the matrix
				int x = monomPairs[j][0];
				int y = monomPairs[j][1];
				int xy = monomPairs[j][2];

				// x-y-z+1
				Polynomial<polyType> poly1(varsTotal);
				poly1.addTerm(1, {x});
				poly1.addTerm(-1, {y});
				poly1.addTerm(-1, {xy});
				poly1.addTerm(1, {oneIndex});
				conPositiveLinear.push_back(poly1);

				// -x+y-z+1
				Polynomial<polyType> poly2(varsTotal);
				poly2.addTerm(-1, {x});
				poly2.addTerm(1, {y});
				poly2.addTerm(-1, {xy});
				poly2.addTerm(1, {oneIndex});
				conPositiveLinear.push_back(poly2);

				// -x-y+z+1
				Polynomial<polyType> poly3(varsTotal);
				poly3.addTerm(-1, {x});
				poly3.addTerm(-1, {y});
				poly3.addTerm(1, {xy});
				poly3.addTerm(1, {oneIndex});
				conPositiveLinear.push_back(poly3);

				// x+y+z+1
				Polynomial<polyType> poly4(varsTotal);
				poly4.addTerm(1, {x});
				poly4.addTerm(1, {y});
				poly4.addTerm(1, {xy});
				poly4.addTerm(1, {oneIndex});
				conPositiveLinear.push_back(poly4);

			// If it's a third-order 
			} else if (monomPairs[j].size() == 7) {

				// Extract the vals and construct the matrix
				int x = monomPairs[j][0];
				int y = monomPairs[j][1];
				int z = monomPairs[j][2];
				int xy = monomPairs[j][3];
				int xz = monomPairs[j][4];
				int yz = monomPairs[j][5];
				int xyz = monomPairs[j][6];

				//1+x+xy+xyz+xz+y+yz+z
				{
					Polynomial<double> newCon(varsTotal);
					newCon.addTerm(1);
					newCon.addTerm(1, {x});
					newCon.addTerm(1, {y});
					newCon.addTerm(1, {z});
					newCon.addTerm(1, {xy});
					newCon.addTerm(1, {xz});
					newCon.addTerm(1, {yz});
					newCon.addTerm(1, {xyz});
					conPositiveLinear.push_back(newCon);
				}

				//1+x+xy-xyz-xz+y-yz-z
				{
					Polynomial<double> newCon(varsTotal);
					newCon.addTerm(1);
					newCon.addTerm(1, {x});
					newCon.addTerm(1, {y});
					newCon.addTerm(-1, {z});
					newCon.addTerm(1, {xy});
					newCon.addTerm(-1, {xz});
					newCon.addTerm(-1, {yz});
					newCon.addTerm(-1, {xyz});
					conPositiveLinear.push_back(newCon);
				}

				//1-x-xy-xyz-xz+y+yz+z
				{
					Polynomial<double> newCon(varsTotal);
					newCon.addTerm(1);
					newCon.addTerm(-1, {x});
					newCon.addTerm(1, {y});
					newCon.addTerm(1, {z});
					newCon.addTerm(-1, {xy});
					newCon.addTerm(-1, {xz});
					newCon.addTerm(1, {yz});
					newCon.addTerm(-1, {xyz});
					conPositiveLinear.push_back(newCon);
				}

				//1+x-xy-xyz+xz-y-yz+z
				{
					Polynomial<double> newCon(varsTotal);
					newCon.addTerm(1);
					newCon.addTerm(1, {x});
					newCon.addTerm(1, {z});
					newCon.addTerm(-1, {y});
					newCon.addTerm(-1, {xy});
					newCon.addTerm(1, {xz});
					newCon.addTerm(-1, {yz});
					newCon.addTerm(-1, {xyz});
					conPositiveLinear.push_back(newCon);
				}

				//1-x+xy+xyz-xz-y-yz+z
				{
					Polynomial<double> newCon(varsTotal);
					newCon.addTerm(1);
					newCon.addTerm(-1, {x});
					newCon.addTerm(-1, {y});
					newCon.addTerm(1, {z});
					newCon.addTerm(1, {xy});
					newCon.addTerm(-1, {xz});
					newCon.addTerm(-1, {yz});
					newCon.addTerm(1, {xyz});
					conPositiveLinear.push_back(newCon);
				}

				//1-x+xy-xyz+xz-y+yz-z
				{
					Polynomial<double> newCon(varsTotal);
					newCon.addTerm(1);
					newCon.addTerm(-1, {x});
					newCon.addTerm(-1, {y});
					newCon.addTerm(-1, {z});
					newCon.addTerm(1, {xy});
					newCon.addTerm(1, {xz});
					newCon.addTerm(1, {yz});
					newCon.addTerm(-1, {xyz});
					conPositiveLinear.push_back(newCon);
				}

				//1+x-xy+xyz-xz-y+yz-z
				{
					Polynomial<double> newCon(varsTotal);
					newCon.addTerm(1);
					newCon.addTerm(1, {x});
					newCon.addTerm(-1, {y});
					newCon.addTerm(-1, {z});
					newCon.addTerm(-1, {xy});
					newCon.addTerm(-1, {xz});
					newCon.addTerm(1, {yz});
					newCon.addTerm(1, {xyz});
					conPositiveLinear.push_back(newCon);
				}

				//1-x-xy+xyz+xz+y-yz-z
				{
					Polynomial<double> newCon(varsTotal);
					newCon.addTerm(1);
					newCon.addTerm(-1, {x});
					newCon.addTerm(1, {y});
					newCon.addTerm(-1, {z});
					newCon.addTerm(-1, {xy});
					newCon.addTerm(1, {xz});
					newCon.addTerm(-1, {yz});
					newCon.addTerm(1, {xyz});
					conPositiveLinear.push_back(newCon);
				}

			}

		}

		// Convert the objective to MOSEK form
		std::vector<polyType> c(varsTotal);
		for (auto const &pair: objLinear.coeffs) {
			int ind = oneIndex;
			if (pair.first != "") {
				ind = std::stoi(pair.first);
			}
			c[ind] = pair.second;
		}
		auto cM = monty::new_array_ptr<polyType>(c);

		// Convert the linear equality constraints to MOSEK form
		std::vector<int> ARows;
		std::vector<int> ACols;
		std::vector<polyType> AVals;
		for (int i=0; i<conZeroLinear.size(); i++) {
			for (auto const &pair: conZeroLinear[i].coeffs) {
				ARows.push_back(i);
				if (pair.first == "") {
					ACols.push_back(oneIndex);
				} else {
					ACols.push_back(std::stoi(pair.first));
				}
				AVals.push_back(pair.second);
			}
		}
		auto AM = mosek::fusion::Matrix::sparse(conZeroLinear.size(), varsTotal, monty::new_array_ptr<int>(ARows), monty::new_array_ptr<int>(ACols), monty::new_array_ptr<polyType>(AVals));

		// Convert the linear positivity constraints to MOSEK form
		std::vector<int> BRows;
		std::vector<int> BCols;
		std::vector<polyType> BVals;
		for (int i=0; i<conPositiveLinear.size(); i++) {
			for (auto const &pair: conPositiveLinear[i].coeffs) {
				BRows.push_back(i);
				if (pair.first == "") {
					BCols.push_back(oneIndex);
				} else {
					BCols.push_back(std::stoi(pair.first));
				}
				BVals.push_back(pair.second);
			}
		}
		auto BM = mosek::fusion::Matrix::sparse(conPositiveLinear.size(), varsTotal, monty::new_array_ptr<int>(BRows), monty::new_array_ptr<int>(BCols), monty::new_array_ptr<polyType>(BVals));

		// Create a model
		mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});

		// DEBUG
		//M->setLogHandler([=](const std::string & msg){std::cout << msg << std::flush;});

		// Create the variable
		mosek::fusion::Variable::t xM = M->variable(varsTotal, mosek::fusion::Domain::inRange(-1, 1));

		// The first element of the vector should be one
		M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));

		// Linear equality constraints
		M->constraint(mosek::fusion::Expr::mul(AM, xM), mosek::fusion::Domain::equalsTo(0.0));

		// Linear positivity constraints
		M->constraint(mosek::fusion::Expr::mul(BM, xM), mosek::fusion::Domain::greaterThan(0));

		// Objective is to minimize the sum of the original linear terms
		M->objective(mosek::fusion::ObjectiveSense::Minimize, mosek::fusion::Expr::dot(cM, xM));

		// Solve the problem
		M->solve();

		// Get the solution values
		auto sol = *(xM->level());
		polyType outer = M->primalObjValue();

		// Output the relevent moments
		std::vector<polyType> solVec(xM->getSize());
		for (int i=0; i<solVec.size(); i++) {
			solVec[i] = sol[i];
		}

		return std::pair<polyType,std::vector<polyType>>(outer, solVec);

	}

	// Given the data from the previous run, guess the next best combo of monoms
	std::vector<int> getBestPair(std::vector<std::string>& monoms, std::vector<Polynomial<polyType>>& monomsAsPolys, std::vector<std::vector<int>>& monomPairs, std::pair<polyType,std::vector<polyType>> prevRes, int matLevel) {

		// Get the list of linear monoms and their values
		std::vector<polyType> linVals(maxVariables);
		int numDone = 0;
		for (int i=0; i<monoms.size(); i++) {
			if (monoms[i].size() == digitsPerInd) {
				if (prevRes.second[i] >= 0) {
					linVals[std::stoi(monoms[i])] = 1;
				} else {
					linVals[std::stoi(monoms[i])] = -1;
				}
				numDone++;
				if (numDone > maxVariables) {
					break;
				}
			}
		}

		// Get the probabilities based on the error for each monomial
		std::vector<polyType> monomResults(monoms.size());
		std::vector<double> monomProbs(monoms.size());
		double totalProb = 0;
		double totalError = 0;
		for (int i=0; i<monoms.size(); i++) {
			if (i < prevRes.second.size()) {
				monomResults[i] = std::abs(prevRes.second[i] - monomsAsPolys[i].evalFast(linVals));
			} else {
				monomResults[i] = 0;
			}
			monomProbs[i] = std::pow(1+monomResults[i],2); // DEBUG
			totalProb += monomProbs[i];
			totalError += monomResults[i];
		}

		// Readjust the probability distribution DEBUG
		for (int i=0; i<monoms.size(); i++) {
			monomProbs[i] = monomProbs[i] / totalProb;
			//std::cout << monoms[i] << "     " << monomResults[i] << "     " << monomProbs[i] << std::endl;
		}

		// Stop if we've found a valid solution
		if (totalError < 1e-5) {
			return {-2,-2,-2};
		}

		// Get the index permutation based on the size
		std::vector<int> orderedIndices(monoms.size(), 0);
		for (int i=0; i<orderedIndices.size(); i++) {
			orderedIndices[i] = i;
		}
		std::sort(orderedIndices.begin(), orderedIndices.end(),
			[&](const int& a, const int& b) {
				return (monoms[a].size() > monoms[b].size());
			}
		);

		// If we're using the second level
		int maxTries = 200;
		if (matLevel == 2) {

			// Now keep searching until we find a good mat that's new
			for (int i2=0; i2<maxTries; i2++) {

				// Pick a random number and then add probs until we reach that number
				double probToReach = (double(rand())/(RAND_MAX));
				int i = -1;
				double probSoFar = 0;
				for (int k=0; k<monoms.size(); k++) {
					probSoFar += monomProbs[k];
					if (probSoFar > probToReach) {
						i = k;
						break;
					}
				}

				// Find the highest order monom which divides this
				for (int j2=0; j2<monoms.size(); j2++) {
					int j = orderedIndices[j2];

					// If it's somewhat appropriate
					if (monoms[j].size() > 0 && monoms[j].size() < monoms[i].size()) {

						// Perform the division
						Polynomial<polyType> otherPoly = monomsAsPolys[i] / monomsAsPolys[j];

						// If it's a nice division
						if (!otherPoly.isNaN) {

							// Get the monoms to add
							std::string firstString = monoms[j];
							std::string secondString = otherPoly.getMonomials()[0];
							std::string combinedString = monoms[i];

							// We know where the first and combined are
							int firstLoc = j;
							int combinedLoc = i;

							// Check if the second exists
							auto secondIter = std::find(monoms.begin(), monoms.end(), secondString);
							int secondLoc = -1;
							bool found = false;

							// If it does exist, see if this combo is new
							if (secondIter != monoms.end()) {
								secondLoc = secondIter - monoms.begin();
								for (int k=0; k<monomPairs.size(); k++) {
									if ((monomPairs[k][0] == firstLoc && monomPairs[k][1] == secondLoc) || (monomPairs[k][0] == secondLoc && monomPairs[k][1] == firstLoc)) {
										found = true;
									}
								}

							// If it doesn't exist, add it
							} else {
								secondLoc = monoms.size();
								monoms.push_back(secondString);
								monomsAsPolys.push_back(Polynomial<polyType>(maxVariables, 1, secondString));
								monomsAsPolys[monomsAsPolys.size()-1].prepareEvalMixed();
							}

							// If it's new, add it
							if (!found) {
								return {firstLoc, secondLoc, combinedLoc};
							}

						}

					}

				}

			}

		// If we're using the third level
		} else if (matLevel == 3) {

			// Now keep searching until we find a good mat that's new
			for (int i2=0; i2<maxTries; i2++) {

				// Pick a random number and then add probs until we reach that number
				double probToReach = (double(rand())/(RAND_MAX));
				int i = -1;
				double probSoFar = 0;
				for (int k=0; k<monoms.size(); k++) {
					probSoFar += monomProbs[k];
					if (probSoFar >= probToReach) {
						i = k;
						break;
					}
				}

				// Perform the division
				Polynomial<polyType> xyzPoly = monomsAsPolys[i];

				// Find the highest order monom which divides this
				for (int j2=0; j2<monoms.size(); j2++) {
					int j = orderedIndices[j2];

					// If it's somewhat appropriate
					if (monoms[j].size() > 0 && monoms[j].size() < monoms[i].size()) {

						// Perform the division
						Polynomial<polyType> yzPoly = xyzPoly / monomsAsPolys[j];

						// If it's a nice division
						if (!yzPoly.isNaN) {

							// Find the highest order monom which divides this
							for (int k2=0; k2<monoms.size(); k2++) {
								int k = orderedIndices[k2];

								// If it's somewhat appropriate
								if (monoms[k].size() > 0 && monoms[k].size() < yzPoly.getMonomials()[0].size()) {

									// Perform the division
									Polynomial<polyType> zPoly = yzPoly / monomsAsPolys[k];

									// If it's a nice division
									if (!zPoly.isNaN) {

										// Calculate the base monomials
										Polynomial<polyType> xPoly = xyzPoly / yzPoly;
										Polynomial<polyType> yPoly = yzPoly / zPoly;

										// Determine all the different multiplications as strings
										std::vector<std::string> monStrings;
										monStrings.push_back(xPoly.getMonomials()[0]);
										monStrings.push_back(yPoly.getMonomials()[0]);
										monStrings.push_back(zPoly.getMonomials()[0]);
										monStrings.push_back((xPoly*yPoly).getMonomials()[0]);
										monStrings.push_back((xPoly*zPoly).getMonomials()[0]);
										monStrings.push_back(yzPoly.getMonomials()[0]);
										monStrings.push_back(xyzPoly.getMonomials()[0]);

										// Find all of these
										std::vector<int> monLocs(monStrings.size());
										for (int k=0; k<monStrings.size(); k++) {
											auto loc = std::find(monoms.begin(), monoms.end(), monStrings[k]);
											if (loc == monoms.end()) {
												monoms.push_back(monStrings[k]);
												monomsAsPolys.push_back(Polynomial<polyType>(maxVariables, 1, monStrings[k]));
												monomsAsPolys[monomsAsPolys.size()-1].prepareEvalMixed();
												monLocs[k] = monoms.size()-1;
											} else {
												monLocs[k] = loc - monoms.begin();
											}
										}

										// Check if it's already been used
										bool isNew = true;
										std::vector<int> list1 = monLocs;
										std::sort(list1.begin(), list1.end());
										for (int k=0; k<monomPairs.size(); k++) {
											std::vector<int> list2 = monomPairs[k];
											std::sort(list2.begin(), list2.end());
											if (list1 == list2) {
												isNew = false;
												break;
											}
										}
										
										// If it's new, that works
										if (isNew) {
											return monLocs;
										}

									}

								}

							}

							// If we didn't find a 3, use the 2
							Polynomial<polyType> xPoly = xyzPoly / yzPoly;

							// Determine all the different multiplications as strings
							std::vector<std::string> monStrings;
							monStrings.push_back(xPoly.getMonomials()[0]);
							monStrings.push_back(yzPoly.getMonomials()[0]);
							monStrings.push_back(xyzPoly.getMonomials()[0]);

							// Find all of these
							std::vector<int> monLocs(monStrings.size());
							for (int k=0; k<monStrings.size(); k++) {
								auto loc = std::find(monoms.begin(), monoms.end(), monStrings[k]);
								if (loc == monoms.end()) {
									monoms.push_back(monStrings[k]);
									monomsAsPolys.push_back(Polynomial<polyType>(maxVariables, 1, monStrings[k]));
									monomsAsPolys[monomsAsPolys.size()-1].prepareEvalMixed();
									monLocs[k] = monoms.size()-1;
								} else {
									monLocs[k] = loc - monoms.begin();
								}
							}

							// Check if it's already been used
							bool isNew = true;
							std::vector<int> list1 = monLocs;
							std::sort(list1.begin(), list1.end());
							for (int k=0; k<monomPairs.size(); k++) {
								std::vector<int> list2 = monomPairs[k];
								std::sort(list2.begin(), list2.end());
								if (list1 == list2) {
									isNew = false;
									break;
								}
							}
							
							// If it's new, that works
							if (isNew) {
								return monLocs;
							}

						}
					}
				}
			}
		}

		// This should never happen
		return {-1, -1, -1};

	}
	
	// Get a lower bound
	polyType lowerBound(int maxIters=100000000, int matLevel=2, bool verbose=false, bool elimAtEnd=false, int matsPerIter=20) {

		// Get the monomial list and sort it
		std::vector<std::string> monoms = getMonomials();
		std::sort(monoms.begin(), monoms.end(), [](const std::string& first, const std::string& second){return first.size() < second.size();});

		// List of semdefinite matrices
		std::vector<std::vector<int>> monomPairs;

		// Random seed
		std::srand(time(0));

		// Single-order monomials should always appear
		for (int i=0; i<maxVariables; i++) {
			std::string newMonom = std::to_string(i);
			newMonom.insert(0, digitsPerInd-newMonom.size(), ' ');
			if (std::find(monoms.begin(), monoms.end(), newMonom) == monoms.end()) {
				monoms.push_back(newMonom);
			}
		}

		// Add the highest order monom
		std::string bigBoi = "";
		for (int i=0; i<maxVariables; i++) {
			std::string newMonom = std::to_string(i);
			newMonom.insert(0, digitsPerInd-newMonom.size(), ' ');
			bigBoi += newMonom;
		}
		if (std::find(monoms.begin(), monoms.end(), bigBoi) == monoms.end()) {
			monoms.push_back(bigBoi);
		}

		// First monom should always be 1
		auto loc = std::find(monoms.begin(), monoms.end(), "");
		if (loc != monoms.end()) {
			monoms.erase(loc);
		}
		monoms.insert(monoms.begin(), "");
		
		// Also get the monomials as polynomials and prepare for fast eval
		std::vector<Polynomial<polyType>> monomsAsPolys(monoms.size());
		for (int i=0; i<monoms.size(); i++) {
			monomsAsPolys[i] = Polynomial<polyType>(maxVariables, 1, monoms[i]);
			monomsAsPolys[i].prepareEvalMixed();
		}

		// Create the mapping from monomials to indices (to linearize)
		std::unordered_map<std::string,std::string> mapping;
		int digitsPerIndAfterLinear = std::ceil(std::log10(monoms.size()+1));
		for (int i=1; i<monoms.size(); i++) {
			std::string newInd = std::to_string(i);
			newInd.insert(0, digitsPerIndAfterLinear-newInd.size(), ' ');
			mapping[monoms[i]] = newInd;
		}

		// Linearize the problem
		Polynomial<polyType> objLinear = obj.replaceWithVariable(mapping);
		std::vector<Polynomial<polyType>> conZeroLinear(conZero.size());
		for (int i=0; i<conZero.size(); i++) {
			conZeroLinear[i] = conZero[i].replaceWithVariable(mapping);
		}
		std::vector<Polynomial<polyType>> conPositiveLinear(conPositive.size());
		for (int i=0; i<conPositive.size(); i++) {
			conPositiveLinear[i] = conPositive[i].replaceWithVariable(mapping);
		}
		
		// Initial solve
		auto prevRes = solveLinear(objLinear, conZeroLinear, conPositiveLinear, monoms, monomPairs);

		// Keep iterating
		std::pair<polyType,std::vector<polyType>> bestVal = {-100000000, {}};
		int highestOrder = 0;
		bool toStop = false;
		int numTheSame = 0;
		for (int i=0; i<maxIters; i++) {

			// Add a monom pair to the list
			if (i >= 1) {

				// Add several at once
				for (int k=0; k<matsPerIter; k++) {

					// Get the best guess
					std::vector<int> monomPairToTry = getBestPair(monoms, monomsAsPolys, monomPairs, prevRes, matLevel);

					// Check for error or convergence
					if (monomPairToTry[0] == -1) {
						continue;
					} else if (monomPairToTry[0] == -2) {
						std::cout << "converged to optimum" << std::endl;
						toStop = true;
						break;
					}

					// Add to various lists
					monomPairs.push_back(monomPairToTry);
					highestOrder = std::max(highestOrder, int(monoms[monomPairToTry[monomPairToTry.size()-1]].size())/digitsPerInd);

				}

				// Break the second loop too
				if (toStop) {
					break;
				}

			}

			// Solve this SDP
			auto res = solveLinear(objLinear, conZeroLinear, conPositiveLinear, monoms, monomPairs);
			prevRes = res;

			// Output
			int totalSize = 0;
			for (int k=0; k<monomPairs.size(); k++) {
				totalSize += 1+monomPairs[k].size();
			}
			std::cout << res.first << "   " << totalSize << "   " << monoms.size() << "   " << highestOrder << "   " << std::endl;

			// If this is similar to the best
			if (std::abs((res.first - bestVal.first) / bestVal.first) < 1e-3) {
				numTheSame++;
			} else {
				numTheSame = 0;
			}

			// Update the best val
			if (res.first > bestVal.first) {
				bestVal = res;
			}

			// If we're stagnating, increase the level
			if (numTheSame > 3 && matLevel < 3) {
				matLevel++;
				std::cout << "increasing matrix level" << std::endl;
			}

		}

		// See how many can be removed without affecting the value
		if (elimAtEnd) {
			for (int k=0; k<monomPairs.size(); k++) {
				auto monomPairsCopy = monomPairs;
				monomPairsCopy.erase(monomPairsCopy.begin()+k);
				auto tempRes = solveLinear(objLinear, conZeroLinear, conPositiveLinear, monoms, monomPairsCopy);
				std::vector<int> tempPair = getBestPair(monoms, monomsAsPolys, monomPairsCopy, tempRes, matLevel);
				std::cout << k << " / " << monomPairs.size() << " = " << tempRes.first << std::endl;
				if (tempPair[0] == -2) {
					monomPairs.erase(monomPairs.begin()+k);
					k--;
				}
			}
			auto finalRes = solveLinear(objLinear, conZeroLinear, conPositiveLinear, monoms, monomPairs);
			int totalSize = 0;
			for (int k=0; k<monomPairs.size(); k++) {
				totalSize += 1+monomPairs[k].size();
			}
			std::cout << "final with " << totalSize << ": " << finalRes.first << std::endl;
		}

		// Output final moment list
		if (verbose) {
			std::cout << std::endl << "final moments:" << std::endl;
			std::cout << monoms << std::endl;
		}

		// Output final matrix list
		if (verbose) {
			std::cout << std::endl << "final mats:" << std::endl;
			std::cout << "{";
			for (int i=0; i<monomPairs.size(); i++) {
				std::cout << "{" << monomPairs[i][0] << "," << monomPairs[i][1] << "," << monomPairs[i][2] << "}";
				if (i < monomPairs.size()-1) {
					std::cout << ",";
				}
			}
			std::cout << "}" << std::endl;
		}

		// Output final matrix list expanded as monoms
		if (verbose) {
			std::cout <<std::endl << "final mats as monoms:" << std::endl;
			std::cout << "{";
			for (int i=0; i<monomPairs.size(); i++) {
				std::cout << "{" << monoms[monomPairs[i][0]] << "," << monoms[monomPairs[i][1]] << "," << monoms[monomPairs[i][2]] << "}";
				if (i < monomPairs.size()-1) {
					std::cout << ",";
				}
			}
			std::cout << "}" << std::endl;
			std::cout << std::endl;
		}

		return bestVal.first;

	}

	// Get a constraint from a monom list and a list of monom polys
	//void addCons(std::vector<std::tuple<std::vector<int>, std::vector<double>, std::vector<std::string>>>& newCons, std::vector<std::string> monoms, std::vector<Polynomial<polyType>> asPolys) {

		//// Calculate the level
		//int L = asPolys.size()-1;

		//// For now we do everything specific to levels
		//if (L == 2) {

			//// Find these polynomials
			//std::vector<std::string> newMonoms;
			//std::vector<int> monomInds(asPolys.size(), -1000000);
			//for (int i=0; i<asPolys.size(); i++) {
				//std::string mon = asPolys[i].getMonomials()[0];
				//auto loc = std::find(monoms.begin(), monoms.end(), mon);
				//if (loc != monoms.end()) {
					//monomInds[i] = loc - monoms.begin();
				//} else {
					//auto locNew = std::find(newMonoms.begin(), newMonoms.end(), mon);
					//if (locNew != newMonoms.end()) {
						//monomInds[i] = -2 -(locNew - newMonoms.begin());
					//} else {
						//newMonoms.push_back(mon);
						//monomInds[i] = -(newMonoms.size()+1);
					//}
				//}
			//}

			//// In this case the -1 index means the constant 1
			//std::vector<int> newInds;
			//newInds.push_back(-1);
			//newInds.push_back(monomInds[1]);
			//newInds.push_back(monomInds[2]);
			//newInds.push_back(monomInds[0]);

			//// Add each con
			//for (int i=0; i<2; i++) {
				//for (int j=0; j<2; j++) {

					//// con(xy,y,x) = 2x(-1)^i + 2y(-1)^j + 2xy(-1)^(i+j) + (1-sumofcoeffs)
					//std::vector<double> newVals;
					//newVals.push_back(1-std::pow(-1, i)-std::pow(-1, j)-std::pow(-1, i+j));
					//newVals.push_back(2*std::pow(-1, i));
					//newVals.push_back(2*std::pow(-1, j));
					//newVals.push_back(2*std::pow(-1, i+j));
					//newCons.push_back(std::make_tuple(newInds, newVals, newMonoms));

				//}
			//}

		//// For level 3
		//} else if (L == 3) {

			//// Multiply to get the mixed
			//asPolys.push_back(asPolys[1]*asPolys[2]);
			//asPolys.push_back(asPolys[1]*asPolys[3]);
			//asPolys.push_back(asPolys[2]*asPolys[3]);

			//// Find these polynomials
			//std::vector<std::string> newMonoms;
			//std::vector<int> monomInds(asPolys.size(), -1000000);
			//for (int i=0; i<asPolys.size(); i++) {
				//std::string mon = asPolys[i].getMonomials()[0];
				//auto loc = std::find(monoms.begin(), monoms.end(), mon);
				//if (loc != monoms.end()) {
					//monomInds[i] = loc - monoms.begin();
				//} else {
					//auto locNew = std::find(newMonoms.begin(), newMonoms.end(), mon);
					//if (locNew != newMonoms.end()) {
						//monomInds[i] = -2 -(locNew - newMonoms.begin());
					//} else {
						//newMonoms.push_back(mon);
						//monomInds[i] = -(newMonoms.size()+1);
					//}
				//}
			//}

			//// In this case the -1 index means the constant 1
			//std::vector<int> newInds;
			//newInds.push_back(-1);
			//newInds.push_back(monomInds[1]);
			//newInds.push_back(monomInds[2]);
			//newInds.push_back(monomInds[3]);
			//newInds.push_back(monomInds[4]);
			//newInds.push_back(monomInds[5]);
			//newInds.push_back(monomInds[6]);
			//newInds.push_back(monomInds[0]);

			//// Add each con
			//for (int i=0; i<2; i++) {
				//for (int j=0; j<2; j++) {
					//for (int k=0; k<2; k++) {

						//// 2x(-1)^i + ... + 2xy(-1)^(i+j) ... + 2xyz(-1)^(i+j+k) + (1-sumofcoeffs)
						//std::vector<double> newVals;
						//newVals.push_back(1-std::pow(-1, i)-std::pow(-1, j)-std::pow(-1, k)-std::pow(-1, i+j)-std::pow(-1, i+k)-std::pow(-1, j+k)-std::pow(-1, i+j+k));
						//newVals.push_back(2*std::pow(-1, i));
						//newVals.push_back(2*std::pow(-1, j));
						//newVals.push_back(2*std::pow(-1, k));
						//newVals.push_back(2*std::pow(-1, i+j));
						//newVals.push_back(2*std::pow(-1, i+k));
						//newVals.push_back(2*std::pow(-1, j+k));
						//newVals.push_back(2*std::pow(-1, i+j+k));
						//newCons.push_back(std::make_tuple(newInds, newVals, newMonoms));

					//}
				//}
			//}

		//}

	//}

	// Get a lower bound with the new method
	//polyType lowerBoundNew(int maxIters=1000000000, int matLevel=2, bool verbose=false, bool elimAtEnd=false, int matsPerIter=20) {

		//// Get the monomial list and sort it
		//std::vector<std::string> monoms = getMonomials();
		//std::sort(monoms.begin(), monoms.end(), [](const std::string& first, const std::string& second){return first.size() < second.size();});

		//// List of semdefinite matrices
		//std::vector<std::vector<int>> monomPairs;

		//// Single-order monomials should always appear
		//for (int i=0; i<maxVariables; i++) {
			//std::string newMonom = std::to_string(i);
			//newMonom.insert(0, digitsPerInd-newMonom.size(), ' ');
			//if (std::find(monoms.begin(), monoms.end(), newMonom) == monoms.end()) {
				//monoms.push_back(newMonom);
			//}
		//}

		//// We don't need a monom for 1
		//auto loc = std::find(monoms.begin(), monoms.end(), "");
		//if (loc != monoms.end()) {
			//monoms.erase(loc);
		//}
		
		//// Also get the monomials as polynomials and prepare for fast eval
		//std::vector<Polynomial<polyType>> monomsAsPolys(monoms.size());
		//for (int i=0; i<monoms.size(); i++) {
			//monomsAsPolys[i] = Polynomial<polyType>(maxVariables, 1, monoms[i]);
			//monomsAsPolys[i].prepareEvalMixed();
		//}

		//// Create the mapping from monomials to indices (to linearize)
		//std::unordered_map<std::string,std::string> mapping;
		//int digitsPerIndAfterLinear = std::ceil(std::log10(monoms.size()+1));
		//mapping[""] = "";
		//for (int i=0; i<monoms.size(); i++) {
			//std::string newInd = std::to_string(i);
			//newInd.insert(0, digitsPerIndAfterLinear-newInd.size(), ' ');
			//mapping[monoms[i]] = newInd;
			////std::cout << monoms[i] << " -> " << newInd << std::endl; // DEBUG
		//}

		//// Increase the number of vars in advance
		//for (int i=0; i<(monoms.size()+conPositive.size())*10; i++) {
			//mapping[std::to_string(monoms.size()+i)+" "] = monoms.size()+i;
		//}

		//// Linearize the problem
		//Polynomial<polyType> objLinear = obj.replaceWithVariable(mapping);
		//std::vector<Polynomial<polyType>> conZeroLinear(conZero.size());
		//for (int i=0; i<conZero.size(); i++) {
			//conZeroLinear[i] = conZero[i].replaceWithVariable(mapping);
		//}
		//std::vector<Polynomial<polyType>> conPositiveLinear(conPositive.size());
		//for (int i=0; i<conPositive.size(); i++) {
			//conPositiveLinear[i] = conPositive[i].replaceWithVariable(mapping);
		//}

		//// All vars should be at most one 
		//for (int i=0; i<monoms.size(); i++) {
			//Polynomial<polyType> newCon(objLinear.maxVariables);
			//newCon.addTerm(-1, {i});
			//newCon.addTerm(1, {});
			//conPositiveLinear.push_back(newCon);
		//}

		//// Add a slack variable to positive cons
		//for (int i=0; i<conPositiveLinear.size(); i++) {
			//Polynomial<polyType> newCon = conPositiveLinear[i];
			//newCon.addTerm(-1, {int(monoms.size())});
			//newCon.addTerm(-1, {});
			//conZeroLinear.push_back(newCon);
			//monoms.push_back("s");
			//monomsAsPolys.push_back(Polynomial<polyType>(objLinear.maxVariables));
		//}

		//// Useful quantities
		//int n = monoms.size();
		//int m = conZeroLinear.size();

		//// Convert cons to a nicer form (Ax = b, -1 <= x <= 1)
		//Eigen::SparseMatrix<double> A(m, n);
		//std::vector<Eigen::Triplet<double>> tripletsA;
		//Eigen::VectorXd b = Eigen::VectorXd::Zero(m);
		//for (int i=0; i<conZeroLinear.size(); i++) {
			//for (auto const &pair: conZeroLinear[i].coeffs) {
				//if (pair.first != "") {
					//tripletsA.push_back(Eigen::Triplet<double>(i, std::stoi(pair.first), pair.second));
				//} else {
					//b(i) = -pair.second;
				//}
			//}
		//}
		//A.setFromTriplets(tripletsA.begin(), tripletsA.end());

		//// Convert obj to a nicer form (c.x+d, -1 <= x <= 1)
		//Eigen::VectorXd c = Eigen::VectorXd::Zero(n);
		//double d = 0;
		//for (auto const &pair: objLinear.coeffs) {
			//if (pair.first != "") {
				//c[std::stoi(pair.first)] = pair.second;
			//} else {
				//d -= pair.second;
			//}
		//}

		//// DEBUG
		////std::cout << "c = {" << c.transpose() << "} + " << d << std::endl;
		////for (int i=0; i<m; i++) {
			////std::cout << "before con " << i << ": {" << A.row(i) << "} = " << b[i] << std::endl;
		////}

		//// Define some matrices
		//Eigen::VectorXd e = Eigen::VectorXd::Ones(n);

		//// Adjust from -1 <= x <= 1
		////          to  0 <= x <= 1
		//b = 0.5*(b + A*e);
		//d = d - c.dot(e);
		//c = 2.0*c;

		//// DEBUG
		////std::cout << std::endl;
		////std::cout << "c = {" << c.transpose() << "} + " << d << std::endl;
		////for (int i=0; i<m; i++) {
			////std::cout << "lin con " << i << ": {" << A.row(i) << "} = " << b[i] << std::endl;
		////}

		//// Initial guess
		//Eigen::VectorXd x = e / 2.0;
		//Eigen::VectorXd s = e / 2.0;
		//Eigen::VectorXd lambda = Eigen::VectorXd::Zero(m);
		////Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> sol;
		////sol.compute(A);
		////Eigen::VectorXd x = sol.solve(b);

		//// Using this:
		//// https://faculty.ksu.edu.sa/sites/default/files/
		//// Interior%20Point%20Methods%20and%20Linear%20Programming.pdf

		//// Keep iterating
		//double mu = 1.0;
		//double rho = 0.1;
		//double primal = 0;
		//bool stepTaken = false;
		//std::cout << std::defaultfloat << std::setprecision(5);
		//std::cout << "----------------------------------------------------------------" << std::endl;
		//std::cout << "| iter |    primal   |    dual     |    alpha    |      mu     |    bonus" << std::endl;
		//std::cout << "----------------------------------------------------------------" << std::endl;
		//for (int i=0; i<maxIters; i++) {

			//// Get the matrix needed for the update
			//Eigen::SparseMatrix<double> matToSolve(n+m+n, n+m+n);
			//std::vector<Eigen::Triplet<double>> tripletsMat;
			//for (int j=0; j<tripletsA.size(); j++) {
				//tripletsMat.push_back(Eigen::Triplet<double>(tripletsA[j].col(), n+tripletsA[j].row(), tripletsA[j].value()));
				//tripletsMat.push_back(Eigen::Triplet<double>(n+tripletsA[j].row(), tripletsA[j].col(), tripletsA[j].value()));
			//}
			//for (int j=0; j<n; j++) {
				//tripletsMat.push_back(Eigen::Triplet<double>(j, n+m+j, 1));
				//tripletsMat.push_back(Eigen::Triplet<double>(n+m+j, j, s[j]));
				//tripletsMat.push_back(Eigen::Triplet<double>(n+m+j, n+m+j, x[j]));
			//}
			//matToSolve.setFromTriplets(tripletsMat.begin(), tripletsMat.end());

			//// Get the vector needed for the update 
			//Eigen::VectorXd vecToSolve = Eigen::VectorXd::Zero(n+m+n);
			//vecToSolve.segment(0, n) = -(A.transpose()*lambda + s - c);
			//vecToSolve.segment(n, m) = -(A*x - b);
			//for (int j=0; j<n; j++) {
				//vecToSolve[n+m+j] = -(x[j]*s[j] - mu);
			//}

			//// Solve this linear system to get the update
			//Eigen::VectorXd delta = Eigen::VectorXd::Zero(n+m+n);
			//Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
			//solver.analyzePattern(matToSolve);
			//solver.factorize(matToSolve);
			//if (solver.info() == Eigen::Success) {
				//delta = solver.solve(vecToSolve);
			//} else {
				//Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
				//solver.compute(matToSolve);
				//delta = solver.solve(vecToSolve);
			//}

			//// DEBUG
			////std::cout << std::endl;
			////std::cout << Eigen::MatrixXd(matToSolve) << std::endl;
			////std::cout << std::endl;
			////std::cout << vecToSolve.transpose() << std::endl;
			////std::cout << std::endl;
			////std::cout << "x = " << x.transpose() << std::endl;
			////std::cout << "s = " << s.transpose() << std::endl;
			////std::cout << "delta = " << delta.transpose() << std::endl;
			////std::cout << std::endl;

			//// Determine the best step-size for x
			//double alphaX = 1.0;
			//bool nonZeroAlphaX = false;
			//for (int j=0; j<1000; j++) {
				//Eigen::VectorXd xMod = x + alphaX*delta.segment(0,n);
				//bool allPositive = true;
				//for (int k=0; k<n; k++) {
					//if (xMod[k] < 0) {
						//allPositive = false;
						//break;
					//}
				//}
				//if (allPositive) {
					//nonZeroAlphaX = true;
					//break;
				//}
				//alphaX *= 0.95;
			//}

			//// Determine the best step-size for s
			//double alphaS = 1.0;
			//bool nonZeroAlphaS = false;
			//for (int j=0; j<1000; j++) {
				//Eigen::VectorXd sMod = s + alphaS*delta.segment(n+m,n);
				//bool allPositive = true;
				//for (int k=0; k<n; k++) {
					//if (sMod[k] < 0) {
						//allPositive = false;
						//break;
					//}
				//}
				//if (allPositive) {
					//nonZeroAlphaS = true;
					//break;
				//}
				//alphaS *= 0.95;
			//}

			//// Per-iter output
			//primal = x.dot(c) + d;
			//double dual = lambda.dot(b) + d;
			//std::cout << "|";
			//std::cout << std::setw(5) << i << " | ";
			//if (stepTaken) {
				//std::cout << std::setw(11) << primal << " | ";
				//std::cout << std::setw(11) << dual << " | ";
			//} else {
				//std::cout << std::setw(11) << "-" << " | ";
				//std::cout << std::setw(11) << "-" << " | ";
			//}
			//if (nonZeroAlphaX) {
				//std::cout << std::setw(11) << alphaX << " | ";
			//} else {
				//std::cout << std::setw(11) << "-" << " | ";
			//}
			//std::cout << std::setw(11) << mu << " | ";
			//std::cout << std::setw(11) << alphaS;
			//std::cout << std::endl;

			//// If we couldn't find a good step-size, try a higher mu
			//if (!nonZeroAlphaX || !nonZeroAlphaS) {
				//mu /= rho;
				//continue;
			//}

			//// Apply the update
			//x += alphaX*delta.segment(0, n);
			//lambda += alphaS*delta.segment(n, m);
			//s += alphaS*delta.segment(n+m, n);
			//stepTaken = true;

			//// Update the barrier parameter if we've done all we can
			//if (vecToSolve.segment(n+m, n).norm() < 100000) {
				//mu *= rho;
			//}

			//// Get the list of rounded values
			//std::vector<polyType> principleValues(obj.maxVariables);
			//int numDone = 0;
			//for (int i=0; i<monoms.size(); i++) {
				//if (monoms[i].size() == digitsPerInd && monoms[i] != "s") {
					//if (x[i] >= 0) {
						//principleValues[std::stoi(monoms[i])] = 1;
					//} else {
						//principleValues[std::stoi(monoms[i])] = -1;
					//}
					//numDone++;
					//if (numDone > maxVariables) {
						//break;
					//}
				//}
			//}

			//// Get the probabilities based on the error for each monomial
			//std::vector<polyType> monomErrors(monoms.size(), 0);
			//std::vector<double> monomProbs(monoms.size(), 0);
			//double totalProb = 0;
			//double totalError = 0;
			//for (int i=0; i<monoms.size(); i++) {
				//if (monoms[i] != "s" && i < x.size()) {
					//monomErrors[i] = std::abs(x[i] - monomsAsPolys[i].evalFast(principleValues));
					//monomProbs[i] = std::pow(monomErrors[i], 2) + std::pow(delta[i], 2); // DEBUG
					//totalProb += monomProbs[i];
					//totalError += monomErrors[i];
				//}
			//}

			//// Readjust the probability distribution
			//for (int i=0; i<monoms.size(); i++) {
				//monomProbs[i] = monomProbs[i] / totalProb;
			//}

			//// Stop if the primal and dual match
			//if (i >= 2 && std::abs(primal-dual) < 1e-4 && vecToSolve.norm() < 1e-5) {
				//break;
			//}

			//// Pick a variable according to this
			//double probToReach = (double(rand())/(RAND_MAX));
			//int monomInd = -1;
			//double probSoFar = 0;
			//for (int k=0; k<monoms.size(); k++) {
				//probSoFar += monomProbs[k];
				//if (probSoFar > probToReach) {
					//monomInd = k;
					//break;
				//}
			//}

			//// L2 constraints
			//std::vector<std::tuple<std::vector<int>, std::vector<double>, std::vector<std::string>>> newCons;
			//for (int k=0; k<monoms.size(); k++) {
				//if (monoms[k] != "s" && monoms[k] != monoms[monomInd]) {
				
					//// Division
					//auto div = monomsAsPolys[monomInd] / monomsAsPolys[k];
					//if (!div.isNaN) {
						////addCons(newCons, monoms, {monomsAsPolys[monomInd], monomsAsPolys[k], div});

						//// L3 constraints
						//for (int l=0; l<monoms.size(); l++) {
							//if (l != k && monoms[l] != "s" && l != monomInd) {

								//// Multiplication
								//auto newBig = monomsAsPolys[monomInd] * monomsAsPolys[l];
								//if (newBig.getDegreePerVariable() == 1) {
									//addCons(newCons, monoms, {newBig, monomsAsPolys[k], div, monomsAsPolys[l]});
								//}

							//}
						//}

					//}

					//// Multiplication
					//auto mul = monomsAsPolys[monomInd] * monomsAsPolys[k];
					//if (mul.getDegreePerVariable() == 1) {
						////addCons(newCons, monoms, {mul, monomsAsPolys[monomInd], monomsAsPolys[k]});

						//// L3 constraints
						//for (int l=0; l<monoms.size(); l++) {
							//if (l != k && monoms[l] != "s" && l != monomInd) {

								//// Multiplication
								//auto newBig = mul * monomsAsPolys[l];
								//if (newBig.getDegreePerVariable() == 1) {
									//addCons(newCons, monoms, {newBig, mul, monomsAsPolys[k], monomsAsPolys[l]});
								//}

							//}
						//}

					//}

				//}
			//}

			//// Test each of these
			//double worstViol = 100000000;
			//int worstInd = -1;
			//for (int k=0; k<newCons.size(); k++) {

				//// Get the con vectors
				//std::vector<int> conInds = std::get<0>(newCons[k]);
				//std::vector<double> conVals = std::get<1>(newCons[k]);

				//// Get the linear sum
				//double viol = 0;
				//for (int m=0; m<conInds.size(); m++) {
					//if (conInds[m] >= 0) {
						//viol += x[conInds[m]]*conVals[m];
					//} else if (conInds[m] == -1) {
						//viol += conVals[m];
					//} else {
						//viol += 2;
					//}
				//}

				//// We want the worst violation
				//if (viol < worstViol) {
					//worstViol = viol;
					//worstInd = k;
				//}

				//// DEBUG
				////std::cout << viol << std::endl;
				////std::cout << conInds << std::endl;
				////std::cout << conVals << std::endl;
				////std::cout << std::get<2>(newCons[k]) << std::endl;
				////std::cout << std::endl;

			//}

			//// DEBUG
			//std::cout << "Ax-b = " << (A*x-b).transpose() << std::endl;
			//std::cout << "x = " << x.transpose() << std::endl;
			//std::cout << "monoms = " << monoms << std::endl;
			//std::cout << "delta = " << delta.segment(0, n).transpose() << std::endl;
			//std::cout << "error = " << monomErrors << std::endl;
			//std::cout << "probs = " << monomProbs << std::endl;
			//std::cout << "chose = " << monoms[monomInd] << std::endl;
			//std::cout << std::endl;

			//if (worstViol > 0) {
				//std::cout << "no new constraints" << std::endl;
				//continue;
			//}

			//// Keep iterating
			//if (i < 10000) {

				//// Get the most constrictive constraint
				//std::vector<int> finalConInds = std::get<0>(newCons[worstInd]);
				//std::vector<double> finalConVals = std::get<1>(newCons[worstInd]);
				//std::vector<std::string> finalConMonoms = std::get<2>(newCons[worstInd]);

				//// DEBUG
				//std::cout << "worst con = " << std::endl;
				//std::cout << worstViol << std::endl;
				//std::cout << finalConInds << std::endl;
				//std::cout << finalConVals << std::endl;
				//std::cout << finalConMonoms << std::endl;

				//// Add the new monoms
				//int oldn = n;
				//int oldm = m;
				//for (int k=finalConMonoms.size()-1; k>=0; k--) {
					//monoms.push_back(finalConMonoms[k]);
					//monomsAsPolys.push_back(Polynomial<polyType>(maxVariables, 1, finalConMonoms[k]));
					//monomsAsPolys[monomsAsPolys.size()-1].prepareEvalMixed();
					//n += 1;
				//}

				//// Add the slack variable
				//monoms.push_back("s");
				//monomsAsPolys.push_back(Polynomial<polyType>(objLinear.maxVariables));
				//n += 1;
				//m += 1;

				//// Expand the various matrices
				//e = Eigen::VectorXd::Ones(n);
				//x.conservativeResize(n);
				//s.conservativeResize(n);
				//lambda.conservativeResize(m);
				//c.conservativeResize(n);
				//b.conservativeResize(m);
				//A.conservativeResize(m, n);

				//// Conservative resize leaves new vals uninitialized
				//x.tail(n-oldn) = Eigen::VectorXd::Ones(n-oldn) / 2.0;
				//x[n-1] = 0.0;
				//s.tail(n-oldn) = Eigen::VectorXd::Ones(n-oldn) / 2.0;
				//lambda.tail(m-oldm) = Eigen::VectorXd::Zero(m-oldm);
				//c.tail(n-oldn) = Eigen::VectorXd::Zero(n-oldn);

				//// Add the relevant values 
				//for (int k=0; k<finalConInds.size(); k++) {
					//if (finalConInds[k] >= 0) {
						//tuple.push_back(Eigen::Triplet<double>(m-1, finalConInds[k], finalConVals[k]));
					//} else if (finalConInds[k] == -1) {
						//b[m-1] = -finalConVals[k];
					//} else {
						//tripletsA.push_back(Eigen::Triplet<double>(m-1, n+finalConInds[k], finalConVals[k]));
					//}
				//}
				//tripletsA.push_back(Eigen::Triplet<double>(m-1, n-1, -1));
				//A.setFromTriplets(tripletsA.begin(), tripletsA.end());

			//}

		//}
		//std::cout << "----------------------------------------------------------------" << std::endl;

		//return primal;

	//}

	// Given the data from the previous run, guess the next best combo of monoms
	std::vector<Polynomial<polyType>> getAllCons(std::vector<std::string>& monoms, std::vector<Polynomial<polyType>>& monomsAsPolys, std::pair<polyType,std::vector<polyType>> prevRes) {
		std::vector<Polynomial<polyType>> allCons;

		// If we're using the second level
		int matLevel = 2;
		double minViol = 100000;
		if (matLevel >= 2) {

			// Consider everything multiplied by everything
			for (int i=1; i<monoms.size(); i++) {
				for (int j=i+1; j<monoms.size(); j++) {

					// Multiply these two polys
					Polynomial<polyType> multi = (monomsAsPolys[i] * monomsAsPolys[j]).removeDuplicates();

					// Find these indices
					std::vector<int> monomIndices = {-1, i, j};
					bool monomNotFound = false;
					auto loc = std::find(monoms.begin(), monoms.end(), multi.getMonomials()[0]);
					if (loc != monoms.end()) {
						monomIndices[0] = loc - monoms.begin();
					} else {
						monomNotFound = true;
					}

					// Don't add new monoms
					if (monomNotFound) {
						continue;
					}

					// The various coefficients of the resulting linear cons
					std::vector<std::vector<int>> levelCoeffs = {
						{1, 1, 1}, 
						{1, -1, -1}, 
						{-1, 1, -1}, 
						{-1, -1, 1}
					};

					// DEBUG3
					//std::cout << monomsAsPolys[i] << " * " << monomsAsPolys[j] << " = " << multi << std::endl;
					//std::cout << prevRes.second[monomIndices[1]] << " * " << prevRes.second[monomIndices[2]] << " - " << prevRes.second[monomIndices[0]] << " = " << prevRes.second[monomIndices[1]]*prevRes.second[monomIndices[2]] - prevRes.second[monomIndices[0]] << std::endl;

					// Check each of these linear cons to see which is violated
					for (int k=0; k<levelCoeffs.size(); k++) {
						double sumValue = 1;
						Polynomial<polyType> sumPoly(maxVariables, 1);
						for (int l=0; l<levelCoeffs[k].size(); l++) {
							if (monomIndices[l] >= 0) { 
								sumValue += levelCoeffs[k][l]*prevRes.second[monomIndices[l]];
								sumPoly += levelCoeffs[k][l]*monomsAsPolys[monomIndices[l]];
							} else {
								sumValue += 0;
								sumPoly += levelCoeffs[k][l]*multi;
							}
						}
						if (sumValue < -1e-2) {
							allCons.push_back(sumPoly);
						}
					}

				}
			}

		}

		// For mat level 3
		if (matLevel >= 3) {

			// Consider everything multiplied by everything
			for (int i=1; i<monoms.size(); i++) {
				for (int j=i+1; j<monoms.size(); j++) {
					for (int k=j+1; k<monoms.size(); k++) {

						// Multiply these polys
						std::vector<Polynomial<polyType>> monomPolys = {
							monomsAsPolys[i],
							monomsAsPolys[j],
							monomsAsPolys[k],
							(monomsAsPolys[i]*monomsAsPolys[j]).removeDuplicates(),
							(monomsAsPolys[i]*monomsAsPolys[k]).removeDuplicates(),
							(monomsAsPolys[j]*monomsAsPolys[k]).removeDuplicates(),
							(monomsAsPolys[i]*monomsAsPolys[j]*monomsAsPolys[k]).removeDuplicates()
						};

						// Find these indices
						std::vector<int> monomIndices = {
							i,
							j,
							k,
							-1,
							-1,
							-1,
							-1
						};
						bool monomNotFound = false;
						for (int l=0; l<monomIndices.size(); l++) {
							if (monomIndices[l] == -1) {
								auto loc = std::find(monoms.begin(), monoms.end(), monomPolys[l].getMonomials()[0]);
								if (loc != monoms.end()) {
									monomIndices[l] = loc - monoms.begin();
								} else {
									monomNotFound = true;
								}
							}
						}

						// Don't add new monoms
						if (monomNotFound) {
							continue;
						}

						// The various coefficients of the resulting linear cons
						std::vector<std::vector<int>> levelCoeffs = {
							{1, 1, 1, 1, 1, 1, 1}, 
							{1, 1, -1, 1, -1, -1, -1}, 
							{1, -1, 1, -1, 1, -1, -1}, 
							{1, -1, -1, -1, -1, 1, 1}, 
							{-1, 1, 1, -1, -1, 1, -1}, 
							{-1, 1, -1, -1, 1, -1, 1}, 
							{-1, -1, 1, 1, -1, -1, 1}, 
							{-1, -1, -1, 1, 1, 1, -1}
						};

						// DEBUG3
						std::cout << monomsAsPolys[i] << " * " << monomsAsPolys[j] << " * " << monomsAsPolys[k] << " = " << monomPolys[6] << std::endl;
						//std::cout << prevRes.second[monomIndices[1]] << " * " << prevRes.second[monomIndices[2]] << " - " << prevRes.second[monomIndices[0]] << " = " << prevRes.second[monomIndices[1]]*prevRes.second[monomIndices[2]] - prevRes.second[monomIndices[0]] << std::endl;

						// Check each of these linear cons to see which is violated
						for (int m=0; m<levelCoeffs.size(); m++) {
							double sumValue = 1;
							Polynomial<polyType> sumPoly(maxVariables, 1);
							for (int l=0; l<levelCoeffs[m].size(); l++) {
								if (monomIndices[l] >= 0) { 
									sumValue += levelCoeffs[m][l]*prevRes.second[monomIndices[l]];
									sumPoly += levelCoeffs[m][l]*monomsAsPolys[monomIndices[l]];
								} else {
									sumValue += 0;
									sumPoly += levelCoeffs[m][l]*monomPolys[l];
								}
							}
							if (sumValue < -1e-2) {
								allCons.push_back(sumPoly);
							}
						}

					}
				}
			}

		}

		return allCons;

	}

	// Given the data from the previous run, guess the next best combo of monoms
	//Polynomial<polyType> getNewCon(std::vector<std::string>& monoms, std::vector<Polynomial<polyType>>& monomsAsPolys, std::pair<polyType,std::vector<polyType>> prevRes) {

		//// If we're using the second level
		//int matLevel = 2;
		//double minViol = 100000;
		//if (matLevel == 2) {

			//// Consider everything multiplied by everything
			//for (int i=1; i<monoms.size(); i++) {
				//for (int j=i+1; j<monoms.size(); j++) {

					//// Multiply these two polys
					//Polynomial<polyType> multi = (monomsAsPolys[i] * monomsAsPolys[j]).removeDuplicates();

					//// Find these indices
					//std::vector<int> monomIndices = {-1, i, j};
					//auto loc = std::find(monoms.begin(), monoms.end(), multi.getMonomials()[0]);
					//if (loc != monoms.end()) {
						//monomIndices[0] = loc - monoms.begin();
					//}

					//if (monomIndices[0] == -1) {
						//continue;
					//}

					//// The various coefficients of the resulting linear cons
					//std::vector<std::vector<int>> levelCoeffs = {
						//{1, 1, 1}, 
						//{1, -1, -1}, 
						//{-1, 1, -1}, 
						//{-1, -1, 1}
					//};

					//// DEBUG3
					////std::cout << monomsAsPolys[i] << " * " << monomsAsPolys[j] << " = " << multi << std::endl;
					////std::cout << prevRes.second[monomIndices[1]] << " * " << prevRes.second[monomIndices[2]] << " - " << prevRes.second[monomIndices[0]] << " = " << prevRes.second[monomIndices[1]]*prevRes.second[monomIndices[2]] - prevRes.second[monomIndices[0]] << std::endl;

					//// Check each of these linear cons to see which is violated
					//for (int k=0; k<levelCoeffs.size(); k++) {
						//double sumValue = 1;
						//Polynomial<polyType> sumPoly(maxVariables, 1);
						//for (int l=0; l<levelCoeffs[k].size(); l++) {
							//if (monomIndices[l] >= 0) { 
								//sumValue += levelCoeffs[k][l]*prevRes.second[monomIndices[l]];
								//sumPoly += levelCoeffs[k][l]*monomsAsPolys[monomIndices[l]];
							//} else {
								//sumValue += 1;
								//sumPoly += levelCoeffs[k][l]*multi;
							//}
						//}
						////std::cout << sumPoly << " = " << sumValue << std::endl; // DEBUG3
						////minViol = std::min(sumValue, minViol); // DEBUG3
						//if (sumValue < -1e-2) {
							////std::cout << "found with " << sumValue << std::endl; // DEBUG3
							//return sumPoly;
						//}
					//}

				//}
			//}

		//// If we're using the third level
		//} else if (matLevel == 3) {

		//}

		//std::cout << "reached end with " << minViol << std::endl;

		//// This should never happen
		//return Polynomial<polyType>();

	//}

	// Simple inefficient sign function
	double sign(double a) {
		if (a >= 0) {
			return 1;
		}
		return -1;
	}
	
	// Get a lower bound
	polyType lowerBound2(int maxIters=100000000, int matLevel=2, bool verbose=false, bool elimAtEnd=false, int matsPerIter=20) {

		// Get the monomial list and sort it
		std::vector<std::string> monoms = getMonomials();
		std::sort(monoms.begin(), monoms.end(), [](const std::string& first, const std::string& second){return first.size() < second.size();});

		// List of semdefinite matrices
		std::vector<std::vector<int>> monomPairs;

		// Random seed
		std::srand(time(0));

		// First monom should always be 1
		auto loc = std::find(monoms.begin(), monoms.end(), "");
		if (loc != monoms.end()) {
			monoms.erase(loc);
		}
		monoms.insert(monoms.begin(), "");

		// Try adding all second order moments
		for (int i=0; i<maxVariables; i++) {
			std::string newMonom1 = std::to_string(i);
			newMonom1.insert(0, digitsPerInd-newMonom1.size(), ' ');
			for (int j=i+1; j<maxVariables; j++) {
				std::string newMonom2 = std::to_string(j);
				newMonom2.insert(0, digitsPerInd-newMonom2.size(), ' ');
				if (std::find(monoms.begin(), monoms.end(), newMonom1+newMonom2) == monoms.end()) {
					monoms.push_back(newMonom1+newMonom2);
				}
			}
		}
		
		// Also get the monomials as polynomials and prepare for fast eval
		std::vector<Polynomial<polyType>> monomsAsPolys(monoms.size());
		for (int i=0; i<monoms.size(); i++) {
			monomsAsPolys[i] = Polynomial<polyType>(maxVariables, 1, monoms[i]);
			monomsAsPolys[i].prepareEvalMixed();
		}

		// Create the mapping from monomials to indices (to linearize)
		std::unordered_map<std::string,std::string> mapping;
		int digitsPerIndAfterLinear = std::ceil(std::log10(monoms.size()+1))*10;
		for (int i=1; i<monoms.size(); i++) {
			std::string newInd = std::to_string(i);
			newInd.insert(0, digitsPerIndAfterLinear-newInd.size(), ' ');
			mapping[monoms[i]] = newInd;
			//std::cout << monoms[i] << " -> " << newInd << std::endl; // DEBUG3
		}

		// Linearize the problem
		Polynomial<polyType> objLinear = obj.replaceWithVariable(mapping);
		std::vector<Polynomial<polyType>> conZeroLinear(conZero.size());
		for (int i=0; i<conZero.size(); i++) {
			conZeroLinear[i] = conZero[i].replaceWithVariable(mapping);
		}
		std::vector<Polynomial<polyType>> conPositiveLinear(conPositive.size());
		for (int i=0; i<conPositive.size(); i++) {
			conPositiveLinear[i] = conPositive[i].replaceWithVariable(mapping);
		}
		
		// Keep iterating
		std::pair<polyType,std::vector<polyType>> bestVal = {-100000000, {}};
		std::pair<polyType,std::vector<polyType>> prevRes = {-100000000, {}};
		for (int i=0; i<maxIters; i++) {

			// Solve this new linear program
			auto res = solveLinear(objLinear, conZeroLinear, conPositiveLinear, monoms, monomPairs);
			prevRes = res;

			// Find an upper bound
			//std::vector<double> roundedVals(monoms.size());
			//for (int k=0; k<monoms.size(); k++) {
				//roundedVals[k] = sign(res.second[k]);
			//}
			//std::vector<double> validSpins(maxVariables, 0);
			//validSpins[0] = 1;
			//for (int j=1; j<maxVariables; j++) {
				//for (int k=1; k<monoms.size(); k++) {
					//int ind1 = std::stoi(monoms[k].substr(0, digitsPerInd));
					//int ind2 = std::stoi(monoms[k].substr(digitsPerInd, digitsPerInd));
					//if (ind1 == j && validSpins[ind2] != 0) {
						//validSpins[ind1] = roundedVals[k]*validSpins[ind2];
					//} else if (ind2 == j && validSpins[ind1] != 0) {
						//validSpins[ind2] = roundedVals[k]*validSpins[ind1];
					//}
				//}
			//}
			//validSpins = {-1, -1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1};
			//bool recreates = true;
			//for (int k=1; k<monoms.size(); k++) {
				//int ind1 = std::stoi(monoms[k].substr(0, digitsPerInd));
				//int ind2 = std::stoi(monoms[k].substr(digitsPerInd, digitsPerInd));
				//if (validSpins[ind1]*validSpins[ind2] != roundedVals[k]) {
					//recreates = false;
				//}
			//}
			//std::cout << monoms << std::endl;
			//std::cout << roundedVals << std::endl;
			//std::cout << validSpins << std::endl;
			//double upperBound = obj.substitute(validSpins)[""];

			// DEBUG3
			double integralityError = 0;
			for (int k=0; k<monoms.size(); k++) {
				integralityError += std::abs(res.second[k]-sign(res.second[k]));
			}

			// Output
			//std::cout << res.first << "   " << upperBound << "   " << monoms.size() << " " << conPositiveLinear.size() << std::endl;
			std::cout << res.first << "   " << integralityError << "   " << monoms.size() << " " << conPositiveLinear.size() << std::endl;

			// Find a constraint that is violated
			//Polynomial<polyType> newCon = getNewCon(monoms, monomsAsPolys, prevRes);
			//std::cout << newCon << std::endl; // DEBUG3

			// Check for error or convergence
			//if (newCon.isNaN) {
				//std::cout << "converged to optimum" << std::endl;
				//break;
			//}

			// Adjust the linear program to remove this variable
			//conPositive.push_back(newCon);
			//conPositiveLinear.push_back(newCon.replaceWithVariable(mapping));

			std::vector<Polynomial<polyType>> allCons = getAllCons(monoms, monomsAsPolys, prevRes);
			for (int j=0; j<allCons.size(); j++) {

				// Add any unknown terms to the mapping
				//std::vector<std::string> newMons = allCons[j].getMonomials();
				//for (int k=0; k<newMons.size(); k++) {
					//auto loc = std::find(monoms.begin(), monoms.end(), newMons[k]);
					//if (loc == monoms.end()) {
						//monoms.push_back(newMons[k]);
						//monomsAsPolys.push_back(Polynomial<polyType>(maxVariables, 1, newMons[k]));
						//monomsAsPolys[monomsAsPolys.size()-1].prepareEvalMixed();
						//std::string newInd = std::to_string(mapping.size());
						//newInd.insert(0, digitsPerIndAfterLinear-newInd.size(), ' ');
						//mapping[newMons[k]] = newInd;
						//std::cout << "added " << newMons[k] << std::endl;
					//}
				//}

				//conPositive.push_back(allCons[j]);
				conPositiveLinear.push_back(allCons[j].replaceWithVariable(mapping));

			}

			if (allCons.size() == 0) {
				break;
			}

		}

		// DEBUG3
		//for (int i=0; i<conZero.size(); i++) {
			//std::cout << conZero[i] << " = " << conZero[i].substitute({-1, -1, -1, -1, 1, -1, 1}) << std::endl;
		//}
		//for (int i=0; i<conPositive.size(); i++) {
			//std::cout << conPositive[i] << " = " << conPositive[i].substitute({-1, -1, 1, 1, 1}) << std::endl;
		//}

		return bestVal.first;

	}

	// Convert from a real problem to a binary problem
	void fromRealProblem(int numBits) {

		// Now we have more variables
		int maxVariablesNew = maxVariables + maxVariables*numBits;

		// Want to replace every x with x1+x2*0.5+x3*0.25-1
		std::vector<int> indsToReplace(maxVariables);
		std::vector<Polynomial<polyType>> polyToReplace(maxVariables);
		for (int i=0; i<maxVariables; i++) {
			indsToReplace[i] = i;
			polyToReplace[i] = Polynomial<polyType>(maxVariablesNew, -1);
			for (int j=0; j<numBits; j++) {
				polyToReplace[i].addTerm(1.0/std::pow(2,j), {maxVariables+numBits*i+j});
			}
		}

		// Perform the replacement
		Polynomial<polyType> objNew = obj.replaceWithPoly(indsToReplace, polyToReplace);
		std::vector<Polynomial<polyType>> conZeroNew(conZero.size());
		for (int i=0; i<conZero.size(); i++) {
			conZero[i] = conZero[i].replaceWithPoly(indsToReplace, polyToReplace);
		}
		std::vector<Polynomial<polyType>> conPositiveNew(conPositive.size());
		for (int i=0; i<conPositive.size(); i++) {
			conPositive[i] = conPositive[i].replaceWithPoly(indsToReplace, polyToReplace);
		}

		// Create the problem
		*this = PolynomialBinaryProblem<polyType>(objNew, conZeroNew, conPositiveNew);

	}

};

// For minimizing a polynomial of vars subject to constraints
template <class polyType>
class PolynomialProblem {
public:

	// The things definiting the problem
	int maxVariables = 1;
	int digitsPerInd = 1;
	Polynomial<polyType> obj;
	std::vector<Polynomial<polyType>> conZero;
	std::vector<Polynomial<polyType>> conPositive;
	std::vector<std::vector<std::pair<int,polyType>>> syms;

	// Constructor with everything but syms
	PolynomialProblem(Polynomial<polyType> obj_, std::vector<Polynomial<polyType>> conZero_, std::vector<Polynomial<polyType>> conPositive_) {
		obj = obj_;
		conZero = conZero_;
		conPositive = conPositive_;
		maxVariables = obj.maxVariables;
		digitsPerInd = obj.digitsPerInd;
	}

	// Constructor with everything 
	PolynomialProblem(Polynomial<polyType> obj_, std::vector<Polynomial<polyType>> conZero_, std::vector<Polynomial<polyType>> conPositive_, std::vector<std::vector<std::pair<int,polyType>>> syms_) {
		obj = obj_;
		conZero = conZero_;
		conPositive = conPositive_;
		syms = syms_;
		maxVariables = obj.maxVariables;
		digitsPerInd = obj.digitsPerInd;
	}

	// Get a list of all the monomials
	std::vector<std::string> getMonomials() {

		// Use an unordered set to get the unique monomial list
		std::unordered_set<std::string> monoms;
		std::vector<std::string> tempList = obj.getMonomials();
		monoms.insert(tempList.begin(), tempList.end());
		for (int i=0; i<conZero.size(); i++) {
			tempList = conZero[i].getMonomials();
			monoms.insert(tempList.begin(), tempList.end());
		}
		for (int i=0; i<conPositive.size(); i++) {
			tempList = conPositive[i].getMonomials();
			monoms.insert(tempList.begin(), tempList.end());
		}

		// Turn this into a vector
		std::vector<std::string> monomList;
		for (const std::string& mon: monoms) {
			monomList.push_back(mon);
		}

		return monomList;

	}

	// Given a list of var indices and values, replace everything
	template <typename otherType>
	PolynomialProblem replaceWithValue(std::vector<int> indsToReplace, std::vector<otherType> valsToReplace) {

		// Convert the list to the correct type
		std::vector<polyType> convertedList(valsToReplace.size());
		for (int i=0; i<valsToReplace.size(); i++) {
			convertedList[i] = polyType(valsToReplace[i]);
		}

		// Start with a blank poly system
		PolynomialProblem newPolyProblem({}, {}, {});

		// Copy each equation, substituting
		newPolyProblem.obj = obj.replaceWithValue(indsToReplace, convertedList);
		for (int i=0; i<conZero.size(); i++) {
			newPolyProblem.conZero.push_back(conZero[i].replaceWithValue(indsToReplace, convertedList));
		}
		for (int i=0; i<conPositive.size(); i++) {
			newPolyProblem.conPositive.push_back(conPositive[i].replaceWithValue(indsToReplace, convertedList));
		}

		return newPolyProblem;

	}

	// When doing std::cout << PolynomialProblem
	friend std::ostream &operator<<(std::ostream &output, const PolynomialProblem &other) {

		// If it's optimisation
		if (other.obj.size() > 0) {
			
			// Output the objective
			output << "Minimize: " << std::endl << std::endl;
			output << other.obj << std::endl << std::endl;

			// Output each constraint
			int numSoFar = 0;
			if (other.conZero.size() + other.conPositive.size() > 0) {
				output << "Subject to: " << std::endl;
				for (int i=0; i<other.conZero.size(); i++) {
					output << std::endl << other.conZero[i] << " = 0 " << std::endl;;
				}
				for (int i=0; i<other.conPositive.size(); i++) {
					output << std::endl << other.conPositive[i] << " > 0 " << std::endl;
				}
			}

		// If it's constraint satisfaction
		} else {

			// Output each constraint
			int numSoFar = 0;
			if (other.conZero.size() + other.conPositive.size() > 0) {
				output << "Find a point satisfying:" << std::endl;
				for (int i=0; i<other.conZero.size(); i++) {
					output << std::endl << other.conZero[i] << " = 0 " << std::endl;;
				}
				for (int i=0; i<other.conPositive.size(); i++) {
					output << std::endl << other.conPositive[i] << " > 0 " << std::endl;
				}
			}

		}

		// Output symmetries
		if (other.syms.size() > 0) {
			output << "With symmetries:" << std::endl;
			for (int i=0; i<other.syms.size(); i++) {
				output << std::endl << other.syms[i] << std::endl;;
			}
		}

		return output;

	}

	// Solve a SDP program given an objective and zero/positive constraints
	std::pair<bool,std::vector<polyType>> solveSDPWithCuts(std::vector<std::pair<double,double>>& varMinMax, std::vector<Polynomial<polyType>>& conZeroLinear, std::vector<Polynomial<polyType>> conPositiveLinear, std::vector<std::string>& monoms, std::vector<std::vector<Polynomial<polyType>>>& monomProducts, std::vector<std::tuple<double,int,int>> qCones={}) {

		// Create the PSD matrices from this list
		std::vector<std::shared_ptr<monty::ndarray<int,1>>> shouldBePSD;
		for (int j=0; j<monomProducts.size(); j++) {

			// Get the list of all monomial locations for the PSD matrix
			std::vector<int> monLocs;
			for (int i=0; i<monomProducts[j].size(); i++) {
				for (int k=0; k<monomProducts[j].size(); k++) {

					// Calculate the product
					std::string monString = (monomProducts[j][i]*monomProducts[j][k]).getMonomials()[0];

					// Find this in the monomial list
					auto loc = std::find(monoms.begin(), monoms.end(), monString);
					if (loc != monoms.end()) {
						monLocs.push_back(loc - monoms.begin());
					} else {
						monLocs.push_back(monoms.size());
						monoms.push_back(monString);
					}

				}
			}

			// This (when reformatted) should be positive-semidefinite
			shouldBePSD.push_back(monty::new_array_ptr<int>(monLocs));

		}

		// Set some vars
		int oneIndex = 0;
		int varsTotal = monoms.size();

		// Get the inds of the first order monomials and their squares
		std::vector<int> firstMonomInds(varMinMax.size());
		std::vector<int> squaredMonomInds(varMinMax.size());
		for (int i=0; i<monoms.size(); i++) {
			if (monoms[i].size() == digitsPerInd) {
				firstMonomInds[std::stoi(monoms[i])] = i;
			} else if (monoms[i].size() == 2*digitsPerInd && monoms[i].substr(0,digitsPerInd) == monoms[i].substr(digitsPerInd,digitsPerInd)) {
				squaredMonomInds[std::stoi(monoms[i].substr(0,digitsPerInd))] = i;
			}
		}

		// Calculate the linear constraints based on the min/maxes 
		for (int i=0; i<varMinMax.size(); i++) {

			// Given two points, find ax+by+c=0
			std::vector<double> point1 = {varMinMax[i].first, varMinMax[i].first*varMinMax[i].first};
			std::vector<double> point2 = {varMinMax[i].second, varMinMax[i].second*varMinMax[i].second};
			double a = point2[1]-point1[1];
			double b = point1[0]-point2[0];
			double c = -a*point1[0]-b*point1[1];

			// Add this as a linear pos con
			Polynomial<polyType> newCon1(varsTotal);
			newCon1.addTerm(a, {firstMonomInds[i]});
			newCon1.addTerm(b, {squaredMonomInds[i]});
			newCon1.addTerm(c, {});
			conPositiveLinear.push_back(newCon1);

		}

		// Convert the linear equality constraints to MOSEK form
		std::vector<int> ARows;
		std::vector<int> ACols;
		std::vector<polyType> AVals;
		for (int i=0; i<conZeroLinear.size(); i++) {
			for (auto const &pair: conZeroLinear[i].coeffs) {
				ARows.push_back(i);
				if (pair.first == "") {
					ACols.push_back(oneIndex);
				} else {
					ACols.push_back(std::stoi(pair.first));
				}
				AVals.push_back(pair.second);
			}
		}
		auto AM = mosek::fusion::Matrix::sparse(conZeroLinear.size(), varsTotal, monty::new_array_ptr<int>(ARows), monty::new_array_ptr<int>(ACols), monty::new_array_ptr<polyType>(AVals));

		// Convert the linear positivity constraints to MOSEK form
		std::vector<int> BRows;
		std::vector<int> BCols;
		std::vector<polyType> BVals;
		for (int i=0; i<conPositiveLinear.size(); i++) {
			for (auto const &pair: conPositiveLinear[i].coeffs) {
				BRows.push_back(i);
				if (pair.first == "") {
					BCols.push_back(oneIndex);
				} else {
					BCols.push_back(std::stoi(pair.first));
				}
				BVals.push_back(pair.second);
			}
		}
		auto BM = mosek::fusion::Matrix::sparse(conPositiveLinear.size(), varsTotal, monty::new_array_ptr<int>(BRows), monty::new_array_ptr<int>(BCols), monty::new_array_ptr<polyType>(BVals));

		// The box constraints, given our box
		std::vector<double> mins(monoms.size(), -1);
		std::vector<double> maxs(monoms.size(), 1);
		for (int i=0; i<varMinMax.size(); i++) {
			mins[firstMonomInds[i]] = varMinMax[i].first;
			mins[squaredMonomInds[i]] = 0;
			maxs[firstMonomInds[i]] = varMinMax[i].second;
			maxs[squaredMonomInds[i]] = std::max(varMinMax[i].first*varMinMax[i].first, varMinMax[i].second*varMinMax[i].second);
		}

		// Create a model
		mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});

		// DEBUG
		//M->setLogHandler([=](const std::string & msg){std::cout << msg << std::flush;});

		// Create the variable
		mosek::fusion::Variable::t xM = M->variable(varsTotal, mosek::fusion::Domain::inRange(monty::new_array_ptr<double>(mins), monty::new_array_ptr<double>(maxs)));

		// The first element of the vector should be one
		M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));

		// Linear equality constraints
		M->constraint(mosek::fusion::Expr::mul(AM, xM), mosek::fusion::Domain::equalsTo(0.0));

		// Linear positivity constraints
		M->constraint(mosek::fusion::Expr::mul(BM, xM), mosek::fusion::Domain::greaterThan(0));

		// SDP constraints
		for (int i=0; i<shouldBePSD.size(); i++) {
			int matDim = std::sqrt(shouldBePSD[i]->size());
			M->constraint(xM->pick(shouldBePSD[i])->reshape(matDim, matDim), mosek::fusion::Domain::inPSDCone(matDim));
		}

		// Quadratic cones
		for (int i=0; i<qCones.size(); i++) {
			M->constraint(mosek::fusion::Expr::vstack(std::sqrt(std::get<0>(qCones[i])), xM->index(firstMonomInds[std::get<1>(qCones[i])]), xM->index(firstMonomInds[std::get<2>(qCones[i])])), mosek::fusion::Domain::inQCone(3));
		}

		// Solve the problem
		M->solve();
		auto statProb = M->getProblemStatus();
		auto statSol = M->getPrimalSolutionStatus();

		// If valid, return this fact
		if (statProb == mosek::fusion::ProblemStatus::PrimalInfeasible || statSol == mosek::fusion::SolutionStatus::Unknown) {
			return std::pair<bool,std::vector<polyType>>(false, {});

		// Otherwise, extract the result
		} else {

			// Get the solution values
			auto sol = *(xM->level());
			polyType outer = M->primalObjValue();

			// Output the relevent moments
			std::vector<polyType> solVec(xM->getSize());
			for (int i=0; i<solVec.size(); i++) {
				solVec[i] = sol[i];
			}

			return std::pair<bool,std::vector<polyType>>(true, solVec);

		}

	}

	// Solve a SDP program given an objective and zero/positive constraints
	std::pair<polyType,std::vector<polyType>> solveSDP(Polynomial<polyType>& objLinear, std::vector<Polynomial<polyType>>& conZeroLinear, std::vector<Polynomial<polyType>> conPositiveLinear, std::vector<std::string>& monoms, std::vector<std::vector<Polynomial<polyType>>>& monomProducts) {

		// Create the PSD matrices from this list
		std::vector<std::shared_ptr<monty::ndarray<int,1>>> shouldBePSD;
		for (int j=0; j<monomProducts.size(); j++) {

			// Get the list of all monomial locations for the PSD matrix
			std::vector<int> monLocs;
			for (int i=0; i<monomProducts[j].size(); i++) {
				for (int k=0; k<monomProducts[j].size(); k++) {

					// Calculate the product
					std::string monString = (monomProducts[j][i]*monomProducts[j][k]).getMonomials()[0];

					// Find this in the monomial list
					auto loc = std::find(monoms.begin(), monoms.end(), monString);
					if (loc != monoms.end()) {
						monLocs.push_back(loc - monoms.begin());
					} else {
						monLocs.push_back(monoms.size());
						monoms.push_back(monString);
					}

				}
			}

			// This (when reformatted) should be positive-semidefinite
			shouldBePSD.push_back(monty::new_array_ptr<int>(monLocs));

		}

		// Set some vars
		int oneIndex = 0;
		int varsTotal = monoms.size();

		// Convert the objective to MOSEK form
		std::vector<polyType> c(varsTotal);
		for (auto const &pair: objLinear.coeffs) {
			int ind = oneIndex;
			if (pair.first != "") {
				ind = std::stoi(pair.first);
			}
			c[ind] = pair.second;
		}
		auto cM = monty::new_array_ptr<polyType>(c);

		// Convert the linear equality constraints to MOSEK form
		std::vector<int> ARows;
		std::vector<int> ACols;
		std::vector<polyType> AVals;
		for (int i=0; i<conZeroLinear.size(); i++) {
			for (auto const &pair: conZeroLinear[i].coeffs) {
				ARows.push_back(i);
				if (pair.first == "") {
					ACols.push_back(oneIndex);
				} else {
					ACols.push_back(std::stoi(pair.first));
				}
				AVals.push_back(pair.second);
			}
		}
		auto AM = mosek::fusion::Matrix::sparse(conZeroLinear.size(), varsTotal, monty::new_array_ptr<int>(ARows), monty::new_array_ptr<int>(ACols), monty::new_array_ptr<polyType>(AVals));

		// Convert the linear positivity constraints to MOSEK form
		std::vector<int> BRows;
		std::vector<int> BCols;
		std::vector<polyType> BVals;
		for (int i=0; i<conPositiveLinear.size(); i++) {
			for (auto const &pair: conPositiveLinear[i].coeffs) {
				BRows.push_back(i);
				if (pair.first == "") {
					BCols.push_back(oneIndex);
				} else {
					BCols.push_back(std::stoi(pair.first));
				}
				BVals.push_back(pair.second);
			}
		}
		auto BM = mosek::fusion::Matrix::sparse(conPositiveLinear.size(), varsTotal, monty::new_array_ptr<int>(BRows), monty::new_array_ptr<int>(BCols), monty::new_array_ptr<polyType>(BVals));

		// Create a model
		mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});

		// DEBUG
		//M->setLogHandler([=](const std::string & msg){std::cout << msg << std::flush;});

		// Create the variable
		mosek::fusion::Variable::t xM = M->variable(varsTotal, mosek::fusion::Domain::inRange(-1, 1));

		// The first element of the vector should be one
		M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));

		// Linear equality constraints
		M->constraint(mosek::fusion::Expr::mul(AM, xM), mosek::fusion::Domain::equalsTo(0.0));

		// Linear positivity constraints
		M->constraint(mosek::fusion::Expr::mul(BM, xM), mosek::fusion::Domain::greaterThan(0));

		// SDP constraints
		for (int i=0; i<shouldBePSD.size(); i++) {
			int matDim = std::sqrt(shouldBePSD[i]->size());
			M->constraint(xM->pick(shouldBePSD[i])->reshape(matDim, matDim), mosek::fusion::Domain::inPSDCone(matDim));
		}

		// Objective is to minimize the sum of the original linear terms
		M->objective(mosek::fusion::ObjectiveSense::Minimize, mosek::fusion::Expr::dot(cM, xM));

		// Solve the problem
		M->solve();

		// Get the solution values
		auto sol = *(xM->level());
		polyType outer = M->primalObjValue();

		// Output the relevent moments
		std::vector<polyType> solVec(xM->getSize());
		for (int i=0; i<solVec.size(); i++) {
			solVec[i] = sol[i];
		}

		return std::pair<polyType,std::vector<polyType>>(outer, solVec);

	}

	// Return combinations 
	std::vector<std::vector<int>> getAllMonomials(int startingVar, int maxVar, int dimension) {

		// Stop when asking for single order monomials
		std::vector<std::vector<int>> toReturn;
		if (dimension == 1) {
			for (int i=startingVar; i<maxVar; i++) {
				toReturn.push_back({i});
			}
			return toReturn;
		}

		// For each var, consider this var as the first and then recurse
		for (int i=startingVar; i<maxVar; i++) {
			std::vector<std::vector<int>> y = getAllMonomials(i, maxVar, dimension-1);
			for (int j=0; j<y.size(); j++) {
				y[j].insert(y[j].begin(), i);
			}
			toReturn.insert(toReturn.end(), y.begin(), y.end());
		}

		return toReturn;

	}

	// Add every monom of a certain order (e.g. for order=2, {12}, {13}, {14}...)
	void addMonomsOfOrder(std::vector<std::string>& monoms, int order) {

		// Get all combinations of this of size=order
		std::vector<std::vector<int>> toAdd = getAllMonomials(0, maxVariables, order);

		// Convert each of these to monomials and add
		for (const auto& a : toAdd) {
			std::string newMonom = "";
			for (const auto& j : a) {
				std::string newSection = std::to_string(j);
				newSection.insert(0, digitsPerInd-newSection.size(), ' ');
				newMonom += newSection;
			}
			if (std::find(monoms.begin(), monoms.end(), newMonom) == monoms.end()) {
				monoms.push_back(newMonom);
			}
		}

	}

	// Find any linear equalities and simplify everything else
	PolynomialProblem<polyType> removeLinear() {

		// Copy this problem
		PolynomialProblem<polyType> newProb = *this;

		// For each equality
		for (int i=0; i<newProb.conZero.size(); i++) {

			// If it's first order (i.e. linear)
			if (newProb.conZero[i].getDegree() == 1) {

				// Make sure there are no zeros, these can cause problems
				newProb.conZero[i] = newProb.conZero[i].prune();

				// Pick a random variable
				std::string varToRemove = "";
				polyType scalingCoeff = 0;
				for (auto const &pair: newProb.conZero[i].coeffs) {
					if (pair.first != "") {
						varToRemove = pair.first;
						scalingCoeff = -pair.second;
						break;
					}
				}

				// The polynomial this var is equal to 
				Polynomial<polyType> equalPoly(maxVariables);
				for (auto const &pair: newProb.conZero[i].coeffs) {
					if (pair.first != varToRemove) {
						equalPoly.addTerm(pair.second/scalingCoeff, pair.first);
					}
				}

				// Replace this
				for (int i=0; i<newProb.conZero.size(); i++) {
					newProb.conZero[i] = newProb.conZero[i].replaceWithPoly(varToRemove, equalPoly);
				}

			}

		}

		return newProb;

	}

	// Collapse to the minimum number of variables
	PolynomialProblem<polyType> collapse() {

		// The map from old indices to new (minified) indices
		std::unordered_map<int,int> indMap;
		int nextInd = 0;
		for (int i=0; i<maxVariables; i++) {

			// Check the objective
			if (obj.contains(i)) {
				indMap[i] = nextInd;
				nextInd++;
			}

			// Check the equality cons if still not found
			if (indMap.find(i) == indMap.end()) {
				for (int j=0; j<conZero.size(); j++) {
					if (conZero[j].contains(i)) {
						indMap[i] = nextInd;
						nextInd++;
						break;
					}
				}
			}

			// Check the inequality cons if still not found
			if (indMap.find(i) == indMap.end()) {
				for (int j=0; j<conPositive.size(); j++) {
					if (conPositive[j].contains(i)) {
						indMap[i] = nextInd;
						nextInd++;
						break;
					}
				}
			}

		}

		// Replace the objective
		int newNumVars = nextInd;
		Polynomial<polyType> newObj = obj.replaceWithVariable(indMap).changeMaxVariables(newNumVars);
		std::vector<Polynomial<polyType>> newConZero;

		// Replace the equality cons
		for (int i=0; i<conZero.size(); i++) {
			Polynomial<polyType> newPoly = conZero[i].replaceWithVariable(indMap).changeMaxVariables(newNumVars);
			if (newPoly.size() > 0) { 
				newConZero.push_back(newPoly);
			}
		}

		// Replace the inequality cons
		std::vector<Polynomial<polyType>> newConPositive;
		for (int i=0; i<conPositive.size(); i++) {
			Polynomial<polyType> newPoly = conPositive[i].replaceWithVariable(indMap).changeMaxVariables(newNumVars);
			if (newPoly.size() > 0) { 
				newConPositive.push_back(newPoly);
			}
		}

		// Create the new poly problem and return
		return PolynomialProblem<polyType>(newObj, newConZero, newConPositive);

	}

	// Attempt find a feasible point of this problem
	std::vector<polyType> findFeasiblePoint(int zeroInd=-1, double alpha=0.1, double tolerance=1e-10, int maxIters=10000000, int threads=4, bool verbose=false) {

		// Combine these to create a single polynomial
		Polynomial<polyType> poly(maxVariables);
		for (int i=0; i<conZero.size(); i++) {
			poly += conZero[i]*conZero[i];
		}

		// If no index specified, add a var and use that
		if (zeroInd == -1) {
			poly = poly.changeMaxVariables(maxVariables+1);
			zeroInd = poly.maxVariables-1;
		}

		// Find a root of this polynomial
		std::vector<polyType> x = poly.findRoot(zeroInd, alpha, tolerance, maxIters, threads, verbose);
		return x;
		
	}

	// Attempt to find a series of constraints that show this is infeasible
	void proveInfeasible(int maxIters=100000000, bool verbose=false) {

		// Get the monomial list and sort it
		std::vector<std::string> monoms = getMonomials();
		std::sort(monoms.begin(), monoms.end(), [](const std::string& first, const std::string& second){return first.size() < second.size();});

		// List of semdefinite matrices
		std::vector<std::vector<Polynomial<polyType>>> monomProducts;

		// Random seed
		std::srand(time(0));

		// Make sure we have all first and second order moments
		addMonomsOfOrder(monoms, 1);

		// First monom should always be 1
		auto loc = std::find(monoms.begin(), monoms.end(), "");
		if (loc != monoms.end()) {
			monoms.erase(loc);
		}
		monoms.insert(monoms.begin(), "");
		int numOGMonoms = monoms.size();
		
		// Also get the monomials as polynomials and prepare for fast eval
		std::vector<Polynomial<polyType>> monomsAsPolys(monoms.size());
		for (int i=0; i<monoms.size(); i++) {
			monomsAsPolys[i] = Polynomial<polyType>(maxVariables, 1, monoms[i]);
			monomsAsPolys[i].prepareEvalMixed();
		}

		// Create the mapping from monomials to indices (to linearize)
		std::unordered_map<std::string,std::string> mapping;
		int digitsPerIndAfterLinear = std::ceil(std::log10(monoms.size()+1));
		for (int i=1; i<monoms.size(); i++) {
			std::string newInd = std::to_string(i);
			newInd.insert(0, digitsPerIndAfterLinear-newInd.size(), ' ');
			mapping[monoms[i]] = newInd;
		}

		// Linearize the problem
		Polynomial<polyType> objLinear = obj.replaceWithVariable(mapping);
		std::vector<Polynomial<polyType>> conZeroLinear(conZero.size());
		for (int i=0; i<conZero.size(); i++) {
			conZeroLinear[i] = conZero[i].replaceWithVariable(mapping);
		}
		std::vector<Polynomial<polyType>> conPositiveLinear(conPositive.size());
		for (int i=0; i<conPositive.size(); i++) {
			conPositiveLinear[i] = conPositive[i].replaceWithVariable(mapping);
		}
		
		// Add first order SDP cons (1, i, j)
		for (int i=0; i<monoms.size(); i++) {
			if (monoms[i].size() == 2*digitsPerInd) {
				int ind1 = std::stoi(monoms[i].substr(0,digitsPerInd));
				int ind2 = std::stoi(monoms[i].substr(digitsPerInd,digitsPerInd));
				monomProducts.push_back({Polynomial<polyType>(maxVariables, 1), Polynomial<polyType>(maxVariables, 1, {ind1}), Polynomial<polyType>(maxVariables, 1, {ind2})});
			}
		}

		// Add second order cone for x^2+y^2=1/d
		int d = 0;
		std::vector<std::tuple<double,int,int>> qCones;
		for (int i=0; i<conZero.size(); i++) {

			// First make sure we have a constant and 3 terms total
			if (conZero[i][""] < 0 && conZero[i].size() == 3) {

				// Get vars and monoms from this
				auto varsInThisPoly = conZero[i].getVariables();
				auto monomsInThisPoly = conZero[i].getMonomials();

				// Check to make sure all vars are squares
				bool allSquares = true;
				for (int j=0; j<monomsInThisPoly.size(); j++) {
					if (monomsInThisPoly[j].size() >= 2*digitsPerInd) {
						if (monomsInThisPoly[j].substr(0,digitsPerInd) != monomsInThisPoly[j].substr(digitsPerInd,digitsPerInd)) { 
							allSquares = false;
						}
					}
				}

				// If all valid, add to the quadratic cones list
				if (allSquares) {
					qCones.push_back({-conZero[i][""], varsInThisPoly[0], varsInThisPoly[1]});
					d = std::round(-1.0 / conZero[i][""]);
					conZero.erase(conZero.begin()+i);
					i--;
				}

			}
		}

		// Start with the most general area
		std::vector<std::vector<std::pair<double,double>>> toProcess;
		std::vector<std::pair<double,double>> varMinMax(obj.maxVariables);
		double overRtD = 1.0 / std::sqrt(d);
		for (int i=0; i<obj.maxVariables; i++) {
			varMinMax[i].first = -overRtD;
			varMinMax[i].second = overRtD;
		}
		toProcess.push_back(varMinMax);

		// Generate a series of random points
		std::vector<std::vector<double>> points;
		int numPointsToGen = 100;
		if (d >= 3) {
			numPointsToGen = 100000;
		}
		for (int i=0; i<numPointsToGen; i++) {
			std::vector<double> newPoint(maxVariables);
			for (int j=0; j<maxVariables; j++) {
				newPoint[j] = 2.0*overRtD*(double(rand())/(RAND_MAX))-overRtD;
			}
			points.push_back(newPoint);
		}

		// Get the inds of the first order monomials and their squares
		std::vector<int> firstMonomInds(maxVariables, -1);
		std::vector<int> squaredMonomInds(maxVariables, -1);
		for (int i=0; i<monoms.size(); i++) {
			if (monoms[i].size() == digitsPerInd) {
				firstMonomInds[std::stoi(monoms[i])] = i;
			} else if (monoms[i].size() == 2*digitsPerInd && monoms[i].substr(0,digitsPerInd) == monoms[i].substr(digitsPerInd,digitsPerInd)) {
				squaredMonomInds[std::stoi(monoms[i].substr(0,digitsPerInd))] = i;
			}
		}

		// Keep splitting until all fail
		int iter = 0;
		auto toProcessBackup = toProcess;
		toProcess = toProcessBackup;
		iter = 0;
		while (toProcess.size() > 0) {

			// Test this subregion
			auto cutRes = solveSDPWithCuts(toProcess[0], conZeroLinear, conPositiveLinear, monoms, monomProducts, qCones);

			// If feasbile, split the region
			if (cutRes.first) {

				// Check the resulting vector for a good place to split TODO
				std::vector<double> errors(maxVariables);
				for (int i=0; i<maxVariables; i++) {
					errors[i] = std::abs(cutRes.second[squaredMonomInds[i]] - cutRes.second[firstMonomInds[i]]*cutRes.second[firstMonomInds[i]]);
				}

				// Find the biggest error TODO could try boltz probabilties
				double biggestError = -10000;
				int bestInd = -1;
				for (int i=0; i<maxVariables; i++) {
					if (errors[i] > biggestError) {
						biggestError = errors[i];
						bestInd = i;
					}
				}

				// Find the biggest section
				//double biggestDiff = -10000;
				//for (int i=0; i<maxVariables; i++) {
					//double diff = toProcess[0][i].second-toProcess[0][i].first;
					//if (diff > biggestDiff) {
						//biggestDiff = diff;
					//}
				//}

				// If we've converged
				if (biggestError < 1e-5) {
					std::cout << "converged in " << iter << " iters to region " << toProcess[0] << std::endl;
					std::cout << "with solution " << cutRes.second << std::endl;
					break;
				}

				// Split it
				double midPoint = (toProcess[0][bestInd].first + toProcess[0][bestInd].second) / 2.0;
				auto copyLeft = toProcess[0];
				auto copyRight = toProcess[0];
				copyLeft[bestInd].second = midPoint;
				copyRight[bestInd].first = midPoint;

				// Add the new paths to the queue
				toProcess.insert(toProcess.begin()+1, copyLeft);
				toProcess.insert(toProcess.begin()+1, copyRight);

			} else {

				// See how many points we can remove
				for (int i=0; i<points.size(); i++) {
					bool inRegion = true;
					for (int k=0; k<toProcess[0].size(); k++) {
						if (points[i][k] < toProcess[0][k].first || points[i][k] > toProcess[0][k].second) {
							inRegion = false;
							break;
						}
					}
					if (inRegion) {
						points.erase(points.begin()+i);
						i--;
					}
				}

			}

			// Per-iteration output
			std::cout << iter << " " << cutRes.first << " " << toProcess.size() << " " << 1.0 - (double(points.size()) / numPointsToGen) << "            \r" << std::flush;
			//std::cout << iter << " " << cutRes.first << " " << toProcess.size() << " " << 1.0 - (double(points.size()) / numPointsToGen) << std::endl;

			// Remove the one we just processed
			toProcess.erase(toProcess.begin());

			// Keep track of the iteration number
			iter++;
			if (iter > 100000000) {
				break;
			}

		}
		std::cout << std::endl;

		// Benchmarks
		// d2n4 58 iterations 0.4s
		// d3n5 estimated ~1400000 iterations in ~10 hours

	}

};

#endif
	
