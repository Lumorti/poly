#ifndef _poly
#define _poly

// Standard includes
#include <limits>
#include <iostream>
#include <chrono>
#include <vector>
#include <complex>
#include <math.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <random>
#include <algorithm>
#include <iomanip>
#include <fstream>

// Use Eigen for matrix/vector ops
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

// OpenMP for parallelisation
#include <omp.h>

// Optim
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

// MOSEK
#include "fusion.h"

// SCS
#include "scs.h"

// Allow complex literals like 1i
using namespace std::complex_literals;

// -----------------------------------------------------------------------------
// --------------------------- Helper Functions --------------------------------
// -----------------------------------------------------------------------------

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
		if (i < vSize-1) {
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

	// Use fixed precision
	int precision = output.precision();
	output << std::showpos << std::fixed;

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
	
	// Reset formatting
	output << std::noshowpos << std::defaultfloat;
	
	return output;

}

// -----------------------------------------------------------------------------
// --------------------------- Polynomial Class --------------------------------
// -----------------------------------------------------------------------------

// Class allowing manipulation of polynomials
template <class polyType>
class Polynomial {

private:

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

	// Same as above but converts to Polynomial
	std::vector<Polynomial> getAllMonomialsAsPoly(int numVars, int dimension) {

		// Get the monomials
		std::vector<std::vector<int>> toReturn;
		for (int d=1; d<=dimension; d++) {
			std::vector<std::vector<int>> nextDim = getAllMonomials(0, numVars, d);
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

	// --------------------------------------------
	//              Constructors
	// --------------------------------------------
	
	// Default constructor (not valid since no variable count given)
	Polynomial() {
		isNaN = true;
		maxVariables = 0;
		digitsPerInd = 0;
	}

	// Simplest constructor
	Polynomial(int maxVariables_) {

		// Copy metadata
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));

	}

	// Constructor with a single term
	Polynomial(int maxVariables_, polyType coeff, std::vector<int> inds) {

		// Copy metadata
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));

		// Add this one term
		addTerm(coeff, inds);

	}

	// Constructor with a single term
	Polynomial(int maxVariables_, polyType coeff, std::initializer_list<int> inds) {

		// Copy metadata
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));

		// Add this one term
		addTerm(coeff, inds);

	}

	// Constructor with a single constant term
	Polynomial(int maxVariables_, polyType coeff) {

		// Copy metadata
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));

		// Add this one term
		addTerm(coeff);

	}

	// Constructor with a single term from string index
	Polynomial(int maxVariables_, polyType coeff, std::string inds) {

		// Copy metadata
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));

		// Add this one term
		addTerm(coeff, inds);

	}

	// Constructor from the string form
	Polynomial(int maxVariables_, std::string asString) {

		// Copy metadata
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

	// --------------------------------------------
	//              Data retrieval
	// --------------------------------------------
	
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

	// --------------------------------------------
	//              Adding terms
	// --------------------------------------------

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

	// --------------------------------------------
	//              Simplification
	// --------------------------------------------

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

	// --------------------------------------------
	//              Replacing
	// --------------------------------------------

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
	Polynomial replaceWithPoly(std::unordered_map<std::string, Polynomial<polyType>> indMap) {

		// Get the number of new max variables from the indMap
		int newMaxVars = 0;
		for (auto const &pair: indMap) {
			newMaxVars = std::max(newMaxVars, pair.second.maxVariables);
		}

		// For each element in this polynomial
		Polynomial newPoly(newMaxVars);
		for (auto const &pair: coeffs) {

			// Loop through the polynomial
			Polynomial toMultiply2(newMaxVars, 1, {});
			std::string remainingIndex = "";
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {

				// See if this index is in the map
				if (indMap.find(pair.first.substr(i, digitsPerInd)) != indMap.end()) {
					toMultiply2 *= indMap[pair.first.substr(i, digitsPerInd)];

				// If not, add it to the remaining index
				} else {
					remainingIndex += pair.first.substr(i, digitsPerInd);
				}

			}

			// Get the remaining index as a polynomial
			Polynomial toMultiply1(newMaxVars, pair.second, remainingIndex);

			// Add to the new poly
			newPoly += toMultiply1*toMultiply2;

		}

		return newPoly;

	}

	// Substitute a variable for a polynomial
	Polynomial replaceWithPoly(std::unordered_map<int, Polynomial<polyType>> indMapInt) {

		// Get the number of new max variables from the indMap
		int newMaxVars = 0;
		for (auto const &pair: indMapInt) {
			newMaxVars = std::max(newMaxVars, pair.second.maxVariables);
		}

		// Convert the int map to a string map
		std::unordered_map<std::string, Polynomial<polyType>> indMap;
		for (auto const &pair: indMapInt) {
			std::string fromString = std::to_string(pair.first);
			fromString.insert(0, digitsPerInd-fromString.size(), ' ');
			indMap[fromString] = pair.second;
		}

		// For each element in this polynomial
		Polynomial newPoly(newMaxVars);
		for (auto const &pair: coeffs) {

			// Loop through the polynomial
			Polynomial toMultiply2(newMaxVars, 1, {});
			std::string remainingIndex = "";
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {

				// See if this index is in the map
				if (indMap.find(pair.first.substr(i, digitsPerInd)) != indMap.end()) {
					toMultiply2 *= indMap[pair.first.substr(i, digitsPerInd)];

				// If not, add it to the remaining index
				} else {
					remainingIndex += pair.first.substr(i, digitsPerInd);
				}

			}

			// Get the remaining index as a polynomial
			Polynomial toMultiply1(newMaxVars, pair.second, remainingIndex);

			// Add to the new poly
			newPoly += toMultiply1*toMultiply2;

		}

		return newPoly;

	}

	// Substitute a int variable for a value  
	Polynomial replaceWithValue(std::unordered_map<int, polyType> indMapInt) {

		// Convert the int map to a string map
		std::unordered_map<std::string, polyType> indMap;
		for (auto const &pair: indMapInt) {
			std::string fromString = std::to_string(pair.first);
			fromString.insert(0, digitsPerInd-fromString.size(), ' ');
			indMap[fromString] = pair.second;
		}

		// For each element in this polynomial
		Polynomial newPoly(maxVariables);
		for (auto const &pair: coeffs) {

			// Loop through the term
			std::string newKey = "";
			polyType newVal = pair.second;
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {

				// See if this index is in the map
				if (indMap.find(pair.first.substr(i, digitsPerInd)) != indMap.end()) {
					newVal *= indMap[pair.first.substr(i, digitsPerInd)];

				// If not, add it to the remaining index
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
	Polynomial replaceWithValue(std::unordered_map<std::string, polyType> indMap) {

		// For each element in this polynomial
		Polynomial newPoly(maxVariables);
		for (auto const &pair: coeffs) {

			// Loop through the term
			std::string newKey = "";
			polyType newVal = pair.second;
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {

				// See if this index is in the map
				if (indMap.find(pair.first.substr(i, digitsPerInd)) != indMap.end()) {
					newVal *= indMap[pair.first.substr(i, digitsPerInd)];

				// If not, add it to the remaining index
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

	// --------------------------------------------
	//              Integration
	// --------------------------------------------

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

	// --------------------------------------------
	//               Addition
	// --------------------------------------------

	// Overload the addition operator with another polynomial
	Polynomial operator+(const Polynomial& other) const {

		// Start with one equation
		Polynomial result = other;

		// Reserve some space for the new map
		result.coeffs.reserve(std::max(coeffs.size(), other.coeffs.size()));

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
	// Overload for addition with constant
	template <class otherType>
	Polynomial operator+(const otherType& other) const {
		
		// If it's new add it, otherwise combine with the existing
		Polynomial newPoly = *this;
		if (newPoly.coeffs.find("") != newPoly.coeffs.end()) {
			newPoly.coeffs[""] += polyType(other);
			if (std::abs(newPoly.coeffs[""]) < newPoly.zeroTol) {
				newPoly.coeffs.erase("");
			}
		} else {
			newPoly.coeffs[""] = polyType(other);
		}

		return newPoly;

	}

	// Overload for in-place addition with constant
	template <class otherType>
	Polynomial& operator+=(const otherType& other) {

		// If it's new add it, otherwise combine with the existing
		if (coeffs.find("") != coeffs.end()) {
			coeffs[""] += polyType(other);
			if (std::abs(coeffs[""]) < zeroTol) {
				coeffs.erase("");
			}
		} else {
			coeffs[""] = polyType(other);
		}

		return *this;

	}

	// Overload for in-place addition with another poly
	Polynomial& operator+=(const Polynomial& other) {

		// For each term of the other
		for (auto const &pair: other.coeffs) {

			// If it's new add it, otherwise combine with the existing
			if (coeffs.find(pair.first) != coeffs.end()) {
				coeffs[pair.first] += pair.second;
				if (std::abs(coeffs[pair.first]) < zeroTol) {
					coeffs.erase(pair.first);
				}
			} else {
				coeffs[pair.first] = pair.second;
			}

		}

		return *this;

	}

	// --------------------------------------------
	//               Subtraction
	// --------------------------------------------
	
	// Overload the self-subtraction operator
	Polynomial operator-() const {
		Polynomial negPoly(maxVariables);
		for (auto const &pair: coeffs) {
			negPoly.coeffs[pair.first] = -pair.second;
		}
		return negPoly;
	}

	// Overload for subtraction with a constant
	template <class otherType>
	Polynomial operator-(const otherType& other) const {
		
		// If it's new add it, otherwise combine with the existing
		Polynomial newPoly = *this;
		if (newPoly.coeffs.find("") != newPoly.coeffs.end()) {
			newPoly.coeffs[""] -= polyType(other);
			if (std::abs(newPoly.coeffs[""]) < newPoly.zeroTol) {
				newPoly.coeffs.erase("");
			}
		} else {
			newPoly.coeffs[""] = polyType(other);
		}

		return newPoly;

	}

	// Overload the addition operator with another polynomial
	Polynomial operator-(const Polynomial& other) const {

		// Start with one equation
		Polynomial result = *this;

		// Reserve some space for the new map
		result.coeffs.reserve(std::max(coeffs.size(), other.coeffs.size()));

		// For each term of the other
		for (auto const &pair: other.coeffs) {

			// If it's new add it, otherwise combine with the existing
			if (result.coeffs.find(pair.first) != result.coeffs.end()) {
				result.coeffs[pair.first] -= pair.second;
				if (std::abs(result.coeffs[pair.first]) < zeroTol) {
					result.coeffs.erase(pair.first);
				}
			} else {
				result.coeffs[pair.first] = -pair.second;
			}

		}

		return result;

	}

	// Overload for in-place subtraction with a constant
	template <class otherType>
	Polynomial& operator-=(const otherType& other) {

		// If it's new add it, otherwise combine with the existing
		if (coeffs.find("") != coeffs.end()) {
			coeffs[""] -= polyType(other);
			if (std::abs(coeffs[""]) < zeroTol) {
				coeffs.erase("");
			}
		} else {
			coeffs[""] = -polyType(other);
		}

		return *this;

	}

	// Overload for in-place subtraction with another poly
	Polynomial& operator-=(const Polynomial& other) {

		// For each term of the other
		for (auto const &pair: other.coeffs) {

			// If it's new add it, otherwise combine with the existing
			if (coeffs.find(pair.first) != coeffs.end()) {
				coeffs[pair.first] -= pair.second;
				if (std::abs(coeffs[pair.first]) < zeroTol) {
					coeffs.erase(pair.first);
				}
			} else {
				coeffs[pair.first] = -pair.second;
			}

		}

		return *this;

	}

	// --------------------------------------------
	//               Multiplication
	// --------------------------------------------

	// Lazy overload for in-place multiplication
	template <class otherType>
	Polynomial& operator*=(const otherType& other) {
		*this = (*this) * other;
		return *this;
	}

	// Overload the multiplication operator with a constant
	template <class otherType>
	Polynomial operator*(const otherType& other) const {

		// For each term of this equation
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

	// --------------------------------------------
	//               Division
	// --------------------------------------------

	// Lazy overload for in-place division
	template <class otherType>
	Polynomial& operator/=(const otherType& other) {
		*this = (*this) / other;
		return *this;
	}

	// Overload the division operator with a constant
	template <class otherType>
	Polynomial operator/(const otherType& other) const {

		// For each term of this equation
		Polynomial result(maxVariables);
		polyType otherConverted = polyType(other);
		for (auto const &pair: coeffs) {

			// Multiply by the constant
			result.coeffs[pair.first] = otherConverted/pair.second;

		}

		return result;

	}
	
	// Overload the division operator with another poly (but only size 1)
	Polynomial operator/(const Polynomial& other) const {

		// If dividing by a more complex poly, return a NaN
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

	// --------------------------------------------
	//             Evaluation
	// --------------------------------------------

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

		// Flag that ready for fast eval now
		fastEvalReady = true;

	}

	// Prepares this equation for ONLY fast eval (saves memory)
	void prepareEvalFastOnly() {
		
		// Prepare for fast eval
		prepareEvalFast();

		// Free some space
		coeffs.clear();

	}

	// Evaluate a polynomial with x values
	template <typename type>
	polyType eval(type x) {

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

	// Evaluate a polynomial with x values, assuming first prepareEvalFast()
	template <typename type>
	polyType evalFast(type x) {

		// For each term in the polynomial
		polyType soFar = 0;
		for (int i=0; i<vals.size(); i++) {

			// Multiply all the values for this term
			polyType sub = 1;
			for (int j=0; j<inds[i].size(); j++) {
				sub *= x[inds[i][j]];
			}

			// Add to the total
			soFar += sub*vals[i];

		}

		return soFar;

	}

	// --------------------------------------------
	//             Output
	// --------------------------------------------

	// Get the output width of the polynomial
	int getOutputWidth() const {

		// Create a copy of the output stream
		std::ostringstream streamCopy;
		streamCopy.copyfmt(std::cout);
		auto startLoc = streamCopy.tellp();

		// Write the poly to the new stream
		streamCopy << *this;

		// The output width is the change in length of this
		auto length = streamCopy.tellp() - startLoc;
		return length;

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
	std::string asAMPL() const {

		// For each term
		std::string toReturn = "";
		int numSoFar = 0;
		for (auto const &pair: coeffs) {
            
            // Make sure it's not zero
            if (std::abs(pair.second) < zeroTol) {
                continue;
            }

			// If it contains at least one variable
			if (pair.first != "") {

				// First the coeff
				if (pair.second != 1) {
					toReturn += std::to_string(pair.second);
					toReturn += "*";
				}

				// Then the indices
				for (int i=0; i<pair.first.size(); i+=digitsPerInd) {
					toReturn += "x";
					toReturn += pair.first.substr(i, digitsPerInd);
					if (i < pair.first.size()-digitsPerInd) {
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
				toReturn += "+";
			}

		}

        // Filter out any spaces
        std::string withoutSpaces = "";
        for (int i=0; i<toReturn.size(); i++) {
            if (toReturn[i] != ' ') {
                withoutSpaces += toReturn[i];
            }
        }
        toReturn = withoutSpaces;

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
		std::stringstream toReturn;
		int numSoFar = 0;
		for (auto const &pair: coeffs) {

			// The coeff and the indices
			toReturn << pair.second << "*{" << pair.first << "}";

			// Output an addition on everything but the last
			numSoFar += 1;
			if (numSoFar < coeffs.size()) {
				toReturn << " + ";
			}

		}

		return toReturn.str();

	}

	// Output, transforming into a trig poly of a single variable
	// (Requires that each pair of variables has a fixed magnitude)
	std::string asMonovariableTrig(std::vector<std::tuple<int,int,double>> norms) {

		// The mapping to transform the vars to trig functions
		std::unordered_map<std::string, std::string> varToTrig;
		std::unordered_map<std::string, double> varMags;
		for (int i=0; i<norms.size(); i++) {

			// One variable will become sin
			std::string padded1 = std::to_string(std::get<0>(norms[i]));
			padded1.insert(0, digitsPerInd-padded1.size(), ' ');
			if (i == 0) {
				varToTrig[padded1] = "sin(x)";
			} else if (i == 1) {
				varToTrig[padded1] = "sin(x*w)";
			} else {
				varToTrig[padded1] = "sin(x*w**"+std::to_string(i)+")";
			}
			varMags[padded1] = std::sqrt(std::get<2>(norms[i]));

			// The other will become cos
			std::string padded2 = std::to_string(std::get<1>(norms[i]));
			padded2.insert(0, digitsPerInd-padded2.size(), ' ');
			if (i == 0) {
				varToTrig[padded2] = "cos(x)";
			} else if (i == 1) {
				varToTrig[padded2] = "cos(x*w)";
			} else {
				varToTrig[padded2] = "cos(x*w**"+std::to_string(i)+")";
			}
			varMags[padded2] = std::sqrt(std::get<2>(norms[i]));

		}

		// For each term
		std::stringstream toReturn;
		toReturn << std::setprecision(16);
		int numSoFar = 0;
		for (auto const &pair: coeffs) {

			std::string vars = "";
			double coeff = pair.second;

			// Loop through each variable
			for (int i=0; i<pair.first.size(); i+=digitsPerInd) {

				// Get the variable
				std::string var = pair.first.substr(i, digitsPerInd);

				// If it's not the constant
				if (var != "") {

					// Get the trig function
					std::string trig = varToTrig[var];

					// Get the magnitude
					double mag = varMags[var];

					// Update everything
					vars += trig;
					if (i < pair.first.size()-digitsPerInd) {
						vars += "*";
					}
					coeff *= mag;

				}

			}

			// Only output if the coeff is non-zero
			numSoFar += 1;
			if (std::abs(coeff) > 1e-13) {

				// The coeff and the indices
				if (vars.size() == 0) {
					toReturn << coeff;
				} else {
					toReturn << coeff << "*" << vars;
				}

				// Output an addition on everything but the last
				if (numSoFar < coeffs.size()) {
					toReturn << " + ";
				}

			}

		}

		return toReturn.str();

	}

	// When doing std::cout << Polynomial
	friend std::ostream &operator<<(std::ostream &output, const Polynomial &other) {
		output << other.asString();
		return output;
	}

	// --------------------------------------------
	//         Assigning and Comparison
	// --------------------------------------------

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

			// Search through the indices
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

	// --------------------------------------------
	//                  Other 
	// --------------------------------------------
	
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

	// Check if the polynomial can be represented as a sum-of-squares
	bool isSOS(int verbosity=1) {

		// Get the degree and half it
		int d = getDegree();
		int dHalf = d / 2;

		// Generate the vector of monoms
		std::vector<int> vars = getVariables();
		std::vector<Polynomial> x = getAllMonomialsAsPoly(vars.size(), dHalf);
		x.insert(x.begin(), Polynomial(maxVariables, 0));
		if (verbosity >= 2) {
			std::cout << x << std::endl;
		}

		// Which matrix elements must sum to each coefficient
		std::unordered_map<std::string,std::vector<std::vector<int>>> cons;
		std::vector<std::vector<int>> cons2(x.size(), std::vector<int>(x.size(), 0));
		for (int i=0; i<x.size(); i++) {
			for (int j=0; j<x.size(); j++) {
				cons[(x[i]*x[j]).getMonomials()[0]].push_back({i,j});
				cons2[i][j] += coeffs[(x[i]*x[j]).getMonomials()[0]];
			}
		}

		// Print the constraints
		if (verbosity >= 1) {
			std::cout << "sdp size: " << x.size() << " by " << x.size() << std::endl;
			std::cout << "num constraints: " << cons.size() << std::endl;
		}
		if (verbosity >= 2) {
			for (auto const &pair: cons) {
				std::cout << "the sum of ";
				for (auto const &pair2: pair.second) {
					std::cout << "(" << pair2[0] << "," << pair2[1] << ") + ";
				}
				std::cout << " must equal " << coeffs[pair.first] << std::endl;
			}
		}

		// Create a model
		mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});

		// Create the variable
		mosek::fusion::Variable::t xM = M->variable(mosek::fusion::Domain::inPSDCone(x.size()));
		mosek::fusion::Variable::t lambda = M->variable(mosek::fusion::Domain::greaterThan(0.0));

		// For each el + el + el = coeff
		for (auto const &pair: cons) {
			if (pair.first == "") {
				M->constraint(mosek::fusion::Expr::add(lambda, mosek::fusion::Expr::sum(xM->pick(monty::new_array_ptr(pair.second)))), mosek::fusion::Domain::equalsTo(coeffs[pair.first]));
			} else {
				M->constraint(mosek::fusion::Expr::sum(xM->pick(monty::new_array_ptr(pair.second))), mosek::fusion::Domain::equalsTo(coeffs[pair.first]));
			}
		}

		// Objective is to see how big a term we can add
		M->objective(mosek::fusion::ObjectiveSense::Maximize, lambda);

		// Solve the problem
		M->solve();

		// If it's infeasible, it's not SOS
		if (M->getProblemStatus() == mosek::fusion::ProblemStatus::PrimalInfeasible || M->getProblemStatus() == mosek::fusion::ProblemStatus::DualInfeasible) {
			return false;
		} else {

			// Ouput lambda
			auto lambdaLevel = *(lambda->level());
			if (verbosity >= 1) {
				std::cout << "lower bound of poly: " << lambdaLevel[0] << std::endl;
			}

			// Output the matrix
			if (verbosity >= 2) {
				auto sol = *(xM->level());
				std::vector<std::vector<double>> sol2(x.size(), std::vector<double>(x.size(), 0));
				for (int i=0; i<x.size(); i++) {
					for (int j=0; j<x.size(); j++) {
						sol2[i][j] = sol[i*x.size()+j];
					}
				}
				std::cout << sol2 << std::endl;
			}

			// Original poly was SOS
			return true;

		}

	}

	// --------------------------------------------
	//              Optimization
	// --------------------------------------------

	// Try to find a root, with one variable being optimized towards zero
	std::vector<polyType> findRoot(int zeroInd=0, double alpha=0.9, double tolerance=1e-10, int maxIters=-1, int threads=4, int verbosity=1, double maxMag=1, double stabilityTerm=1e-13, std::vector<polyType> startX={}) {
		return integrate(zeroInd).findStationaryPoint(alpha, tolerance, maxIters, threads, verbosity, maxMag, stabilityTerm, startX);
	}

	// Cost/gradient function for optim
	static double gradFunction(const Eigen::VectorXd& x, Eigen::VectorXd* gradOut, void* optData) {
		
		// Recast this generic pointer into the correct format
		struct dataType {
			int maxVariables;
			Polynomial obj;
			std::vector<Polynomial> gradient;
		};
		dataType* optDataRecast = reinterpret_cast<dataType*>(optData);

		// Calculate the gradient
		if (gradOut) {
			for (int i=0; i<optDataRecast->maxVariables; i++) {
				(*gradOut)(i) = optDataRecast->gradient[i].evalFast(x);
			}
		}

		// Return the evaluated objective function
		return (optDataRecast->obj).evalFast(x);

	}

	// Use the optim library to minimize the polynomial
	std::vector<polyType> minimize(double alpha=0.9, double tolerance=1e-10, int maxIters=-1, int threads=4, int verbosity=1, double maxMag=1, double stabilityTerm=1e-13) {

		// Prepare everything for parallel computation
		omp_set_num_threads(threads);
		Eigen::setNbThreads(threads);

		// Start a timer
		auto begin = std::chrono::steady_clock::now();

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
		prepareEvalFast();
		#pragma omp parallel for
		for (int i=0; i<maxVariables; i++) {
			gradient[i].prepareEvalFast();
			for (int j=i; j<maxVariables; j++) {
				hessian[i][j].prepareEvalFast();
			}
		}

		// Check a bunch of random xs to see which is the best
		Eigen::VectorXd x = maxMag*Eigen::VectorXd::Random(maxVariables);

		// The format of a our data object
		struct dataType {
			int maxVariables;
			Polynomial obj;
			std::vector<Polynomial> gradient;
		};

		// The data object to be passed to the gradient function
		dataType optData;
		optData.maxVariables = maxVariables;
		optData.gradient = gradient;
		optData.hessian = hessian;
		optData.obj = *this;
 
		// Settings for optim
		optim::algo_settings_t settings;
		settings.print_level = 0;
		settings.iter_max = 100000;
		settings.grad_err_tol = 1e-08;
		settings.rel_sol_change_tol = 1e-14;
		settings.rel_objfn_change_tol = 1e-08;

		// Run the optimisation from the optim library
		bool success = optim::lbfgs(x, gradFunction, &optData, settings);

		// Convert to a std::vector and return
		std::vector<polyType> toReturn2(maxVariables);
		for (int i=0; i<maxVariables; i++) {
			toReturn2[i] = x(i);
		}
		return toReturn2;

	}

	// Use the Newton method to find a stationary point
	std::vector<polyType> findStationaryPoint(double alpha=0.9, double tolerance=1e-10, int maxIters=-1, int threads=4, int verbosity=1, double maxMag=1, double stabilityTerm=1e-13, std::vector<polyType> startX={}) {

		// Prepare everything for parallel computation
		omp_set_num_threads(threads);
		Eigen::setNbThreads(threads);

		// Start a timer
		auto begin = std::chrono::steady_clock::now();

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
		prepareEvalFast();
		#pragma omp parallel for
		for (int i=0; i<maxVariables; i++) {
			gradient[i].prepareEvalFast();
			for (int j=i; j<maxVariables; j++) {
				hessian[i][j].prepareEvalFast();
			}
		}

		// Start from a given point
		Eigen::VectorXd x = maxMag*Eigen::VectorXd::Random(maxVariables);
		if (startX.size() > 0) {
			for (int i=0; i<startX.size(); i++) {
				x(i) = startX[i];
			}
		}

		// Perform gradient descent using this info
		Eigen::MatrixXd inv(maxVariables, maxVariables);
		Eigen::VectorXd p(maxVariables);
		Eigen::MatrixXd H(maxVariables, maxVariables);
		Eigen::VectorXd g(maxVariables);
		Eigen::VectorXd bestX = x;
		double maxX = 0;
		int iter = 0;
		int precisionIters = 0;
		double norm = 1;
		double prevNorm = norm;
		double minVal = 10000000;
		double minValPrePrecision = minVal;
		std::vector<Eigen::VectorXd> pastX;
		while (true) {

			// Increse the iter count, also the precision count if we're in precision mode
			if (precisionIters == 0) {
				iter++;
			} else {
				precisionIters++;
			}

			// Convert x to a C++ array for faster eval performance
			double *xFast = x.data();

			// Calculate the gradient
			#pragma omp parallel for
			for (int i=0; i<maxVariables; i++) {
				g(i) = gradient[i].evalFast(xFast);
			}
			prevNorm = norm;
			norm = std::abs(g.norm());

			// Calculate the Hessian
			#pragma omp parallel for
			for (int i=0; i<maxVariables; i++) {
				for (int j=i; j<maxVariables; j++) {
					H(i,j) = hessian[i][j].evalFast(xFast);
					H(j,i) = H(i,j);
				}
			}

			// Add some diagonal for a bit of numerical stability
			if (stabilityTerm >= 0) {
				for (int i=0; i<maxVariables; i++) {
					H(i,i) += stabilityTerm;
				}
			}

			// Determine the direction
			if (precisionIters == 0) {
				p = -H.partialPivLu().solve(g);
			} else {
				p = -H.fullPivHouseholderQr().solve(g);
			}

			// Perform the update
			x += alpha*p;

			// Keep track of the best we've found
			if (norm < minVal) {
				minVal = norm;
				bestX = x;
			}

			// Per-iteration output
			if (precisionIters == 0) {
				if (verbosity >= 2) {
					std::cout << iter << " " << norm << " " << minVal << "\n" << std::flush;
				} else if (verbosity >= 1) {
					std::cout << iter << " " << norm << " " << minVal << "          \r" << std::flush;
				}
			} else {
				if (verbosity >= 2) {
					std::cout << (precisionIters-1) << " (precision) " << norm << " " << minVal << " (from " << minValPrePrecision << ")\n" << std::flush;
				} else if (verbosity >= 1) {
					std::cout << (precisionIters-1) << " (precision) " << norm << " " << minVal << " (from " << minValPrePrecision << ")          \r" << std::flush;
				}
			}

			// Before finishing, see if we can improve the solution a bit
			if (precisionIters == 0 && (norm < tolerance || (iter > maxIters && maxIters > 0))) {
				x = bestX;
				alpha = 0.95;
				stabilityTerm = -1;
				precisionIters = 1;
				minValPrePrecision = minVal;
			}

			// Stop when we've improved the solution a bit
			if (precisionIters > 100 || maxIters == 1) {
				break;
			}

			// Jump if we're stalling a bit
			if (precisionIters == 0 && (p.norm() <= 1e-10 || norm > 1e20 || isnan(norm) || iter % 10000 == 0)) {
				x = maxMag*Eigen::VectorXd::Random(maxVariables);
			}

		}
		x = bestX;

		// If asking for everything
		if (verbosity >= 2) {
			std::cout << std::endl;
			std::cout << "x: " << x.transpose() << std::endl;
			std::cout << std::endl;
			std::cout << "grad: " << g.transpose() << std::endl;
			std::cout << std::endl;
		}

		// Stop the timer and report
		std::cout << std::defaultfloat;
		auto end = std::chrono::steady_clock::now();
		double timeTaken = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000.0;
		if (verbosity == 1) {
			std::cout << std::endl;
		}
		if (verbosity >= 1) {
			std::cout << "finished in " << iter << " iterations (" << timeTaken << " seconds total)" << std::endl;
		}

		// Convert the eigen vec into a normal vec
		std::vector<polyType> toReturn(maxVariables);
		for (int i=0; i<maxVariables; i++) {
			toReturn[i] = x(i);
		}
		return toReturn;

	}

};

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
	return (-poly)+other;
}

// Specific overload for double
Polynomial<double> real(const Polynomial<std::complex<double>>& poly) {
	Polynomial<double> newPoly(poly.maxVariables);
	for (auto const &pair: poly.coeffs) {
		if (std::abs(std::real(pair.second)) > poly.zeroTol) {
			newPoly.coeffs[pair.first] = std::real(pair.second);
		}
	}
	return newPoly;
}

// Specific overload for double
Polynomial<double> imag(const Polynomial<std::complex<double>>& poly) {
    Polynomial<double> newPoly(poly.maxVariables);
    for (auto const &pair: poly.coeffs) {
        if (std::abs(std::imag(pair.second)) > poly.zeroTol) {
            newPoly.coeffs[pair.first] = std::imag(pair.second);
        }
    }
    return newPoly;
}

// Overload std::real
template <typename polyType2, typename polyType>
Polynomial<polyType2> std::real(const Polynomial<polyType>& poly) {
	Polynomial<polyType2> newPoly(poly.maxVariables);
	for (auto const &pair: poly.coeffs) {
		if (std::abs(std::real(pair.second)) > poly.zeroTol) {
			newPoly.coeffs[pair.first] = std::real(pair.second);
		}
	}
	return newPoly;
}

// Overload std::imag
template <typename polyType2, typename polyType>
Polynomial<polyType2> std::imag(const Polynomial<polyType>& poly) {
	Polynomial<polyType2> newPoly(poly.maxVariables);
	for (auto const &pair: poly.coeffs) {
		if (std::abs(std::imag(pair.second)) > poly.zeroTol) {
			newPoly.coeffs[pair.first] = std::imag(pair.second);
		}
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

// ----------------------------------------------------------------------------------
// --------------------------- PolynomialProblem Class ------------------------------
// ----------------------------------------------------------------------------------

// For minimizing a polynomial of vars subject to constraints
template <class polyType>
class PolynomialProblem {

private:

	// Calculate the equation of a line given 2 points
	// In the form a + bx + cy = 0
	std::vector<double> getLineFromPoints(std::vector<double> point1, std::vector<double> point2) {
		double b = point2[1]-point1[1];
		double c = point1[0]-point2[0];
		double a = -b*point1[0]-c*point1[1];
		return {a, b, c};
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

	// Make a number of a seconds look nicer for output
	std::string representTime(double numSeconds) {
		if (numSeconds <= 0) {
			return "?s";
		} else if (numSeconds < 1) {
			return std::to_string(int(numSeconds * 1.0e3)) + "ms";
		} else if (numSeconds < 60) {
			return std::to_string(int(numSeconds)) + "s";
		} else if (numSeconds < 3600) {
			return std::to_string(int(numSeconds / 60.0)) + "m";
		} else if (numSeconds < 86400) {
			return std::to_string(int(numSeconds / 3600.0)) + "h";
		} else if (numSeconds < 604800) {
			return std::to_string(int(numSeconds / 86400.0)) + "d";
		} else if (numSeconds < 2.628e6) {
			return std::to_string(int(numSeconds / 604800.0)) + "w";
		} else if (numSeconds < 3.156e+7) {
			return std::to_string(int(numSeconds / 2.628e6)) + "mo";
		} else {
			double numYears = numSeconds / 3.156e+7;
			int mag = int(std::log10(numYears));
			if (mag < 0) {
				return "?s";
			} else if (numYears < 1e6) {
				return std::to_string(long(numYears)) + "y";
			} else {
				return "10^" + std::to_string(mag) + "y";
			}
		}
	}

	// Same as above but converts to Polynomial
	std::vector<Polynomial<double>> getAllMonomialsAsPoly(int numVars, int dimension) {

		// Get the monomials
		std::vector<std::vector<int>> toReturn;
		for (int d=1; d<=dimension; d++) {
			std::vector<std::vector<int>> nextDim = getAllMonomials(0, numVars, d);
			toReturn.insert(toReturn.end(), nextDim.begin(), nextDim.end());
		}

		// Convert to Polynomial
		std::vector<Polynomial<double>> toReturnPoly(toReturn.size(), Polynomial<double>(numVars));
		for (int i=0; i<toReturn.size(); i++) {
			toReturnPoly[i].addTerm(1, toReturn[i]);
		}

		return toReturnPoly;

	}

public:

	// The things definiting the problem
	int maxVariables = 1;
	int digitsPerInd = 1;
	Polynomial<polyType> obj;
	std::vector<Polynomial<polyType>> conZero;
	std::vector<Polynomial<polyType>> conPositive;
	std::vector<std::vector<Polynomial<polyType>>> conPSD;
	std::vector<bool> varIsBinary;
	std::vector<std::pair<polyType,polyType>> varBounds;
	double zeroTol = 1e-15;
	double largeTol = 1e10;

	// Default constructor
	PolynomialProblem() {}

	// Main constructor 
	PolynomialProblem(int maxVariables_) {
		maxVariables = maxVariables_;
		digitsPerInd = std::ceil(std::log10(maxVariables+1));
		varIsBinary = std::vector<bool>(maxVariables, false);
		varBounds = std::vector<std::pair<polyType,polyType>>(maxVariables, std::make_pair(-largeTol, largeTol));
		obj = Polynomial<polyType>(maxVariables);
	}

	// Constructor from other problem
	PolynomialProblem(const PolynomialProblem& other) {
		maxVariables = other.maxVariables;
		digitsPerInd = other.digitsPerInd;
		obj = other.obj;
		conZero = other.conZero;
		conPositive = other.conPositive;
		conPSD = other.conPSD;
		zeroTol = other.zeroTol;
		varIsBinary = other.varIsBinary;
		varBounds = other.varBounds;
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
		for (int i=0; i<conPSD.size(); i++) {
			for (int j=0; j<conPSD[i].size(); j++) {
				tempList = conPSD[i][j].getMonomials();
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

	// Given a list of var indices and values, replace everything
	template <typename otherType>
	PolynomialProblem replaceWithValue(std::unordered_map<int,otherType> map) {

		// Convert the map to the correct type
		std::unordered_map<int,polyType> convertedMap;
		for (auto it=map.begin(); it!=map.end(); it++) {
			convertedMap[it->first] = polyType(it->second);
		}

		// Create a copy of the problem
		PolynomialProblem newProb(*this);

		// Substitute in the values for the objective
		newProb.obj = newProb.obj.replaceWithValue(convertedMap);

		// Substitute in the values for the zero constraints
		for (int i=0; i<conZero.size(); i++) {
			newProb.conZero[i] = newProb.conZero[i].replaceWithValue(convertedMap);
		}

		// Substitute in the values for the positivity constraints
		for (int i=0; i<newProb.conPositive.size(); i++) {
			newProb.conPositive[i] = newProb.conPositive[i].replaceWithValue(convertedMap);
		}

		// Substitute in the values for the PSD constraint
		for (int i=0; i<conPSD.size(); i++) {
			for (int j=0; j<conPSD[i].size(); j++) {
				newProb.conPSD[i][j] = newProb.conPSD[i][j].replaceWithValue(convertedMap);
			}
		}

		return newProb;

	}

	// Apply a map from variables to variables
	PolynomialProblem<polyType> replaceWithVariable(std::unordered_map<int,int> indMap) {

		// Copy the problem
		int newNumVars = indMap.size();
		PolynomialProblem<polyType> newProb(newNumVars);

		// Copy the bounds that are still relevant
		newProb.varIsBinary = std::vector<bool>(newNumVars, false);
		newProb.varBounds = std::vector<std::pair<polyType,polyType>>(newNumVars, std::make_pair(-largeTol, largeTol));
		for (auto it=indMap.begin(); it!=indMap.end(); it++) {
			newProb.varIsBinary[it->second] = varIsBinary[it->first];
			newProb.varBounds[it->second] = varBounds[it->first];
		}

		// Replace the objective
		newProb.obj = obj.replaceWithVariable(indMap).changeMaxVariables(newNumVars);

		// Replace the equality cons
		for (int i=0; i<conZero.size(); i++) {
			Polynomial<polyType> newPoly = conZero[i].replaceWithVariable(indMap).changeMaxVariables(newNumVars);
			if (newPoly.size() > 0) { 
				newProb.conZero.push_back(newPoly);
			}
		}

		// Replace the inequality cons
		for (int i=0; i<conPositive.size(); i++) {
			Polynomial<polyType> newPoly = conPositive[i].replaceWithVariable(indMap).changeMaxVariables(newNumVars);
			if (newPoly.size() > 0) { 
				newProb.conPositive.push_back(newPoly);
			}
		}

		// Replace the PSD con
		for (int i=0; i<conPSD.size(); i++) {
			newProb.conPSD.push_back({});
			for (int j=0; j<conPSD[i].size(); j++) {
				Polynomial<polyType> newPoly = conPSD[i][j].replaceWithVariable(indMap).changeMaxVariables(newNumVars);
				newProb.conPSD[i].push_back(newPoly);
			}
		}

		// Create the new poly problem and return
		return newProb;

	}

	// When doing std::cout << PolynomialProblem
	friend std::ostream &operator<<(std::ostream &output, const PolynomialProblem &other) {

		// If it's optimisation
		if (other.obj.size() > 0) {
			output << "Task - Minimize:" << std::endl;
			output << other.obj << std::endl << std::endl;

		// If it's constraint satisfaction
		} else {
			output << "Task - Find a satisfying point" << std::endl << std::endl;
		}

		// If any variable is bounded, then output the bounds
		bool hasBounds = false;
		for (int i=0; i<other.maxVariables; i++) {
			if (other.varIsBinary[i]) {
				output << "Variable " << i << " is binary in {" << other.varBounds[i].first << ", " << other.varBounds[i].second << "}" << std::endl;
				hasBounds = true;
			} else if (other.varBounds[i].first > -other.largeTol || other.varBounds[i].second < other.largeTol) {
				output << "Variable " << i << " is in range [" << other.varBounds[i].first << ", " << other.varBounds[i].second << "]" << std::endl;
				hasBounds = true;
			}
		}
		if (hasBounds) {
			output << std::endl;
		}

		// Output each constraint
		if (other.conZero.size() + other.conPositive.size() + other.conPSD.size() > 0) {
			output << "Subject to: " << std::endl;
			for (int i=0; i<other.conZero.size(); i++) {
				output << std::endl << other.conZero[i] << " = 0 " << std::endl;
			}
			for (int i=0; i<other.conPositive.size(); i++) {
				output << std::endl << other.conPositive[i] << " > 0 " << std::endl;
			}
			if (other.conPSD.size() > 0) {
				output << std::endl << other.conPSD << " >= 0" << std::endl;
			}
		}

		// Output the number of variables and constraints
		output << "Total variables: " << other.maxVariables << std::endl;
		output << "Total constraints: " << other.conZero.size() + other.conPositive.size() + other.conPSD.size();

		return output;

	}

	// Return the max power of any monomial of any equation
	int getDegree() const {
		int maxDegree = obj.getDegree();
		for (auto const &eqn : conZero) {
			maxDegree = std::max(maxDegree, eqn.getDegree());
		}
		for (auto const &eqn : conPositive) {
			maxDegree = std::max(maxDegree, eqn.getDegree());
		}
		for (auto const &eqn1 : conPSD) {
			for (auto const &eqn2 : eqn1) {
				maxDegree = std::max(maxDegree, eqn2.getDegree());
			}
		}
		return maxDegree;
	}

	// Find any linear equalities and simplify everything else
	PolynomialProblem<polyType> removeLinear() {

		// Copy this problem
		PolynomialProblem<polyType> newProb = *this;

		// For each equality
		for (int i=0; i<newProb.conZero.size(); i++) {

			// Make sure there are no zeros, these can cause problems
			newProb.conZero[i] = newProb.conZero[i].prune();

			// If it's first order (i.e. linear)
			if (newProb.conZero[i].getDegree() == 1) {

				// Pick a random variable
				std::string varToRemove = "";
				int bestInd = -1;
				polyType scalingCoeff = 0;
				for (auto const &pair: newProb.conZero[i].coeffs) {
					if (pair.first != "") {
						if (std::stoi(pair.first) > bestInd) {
							varToRemove = pair.first;
							scalingCoeff = -pair.second;
							bestInd = std::stoi(pair.first);
						}
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
				std::unordered_map<std::string,Polynomial<polyType>> map;
				map[varToRemove] = equalPoly;
				newProb.obj = newProb.obj.replaceWithPoly(map);
				for (int j=0; j<newProb.conZero.size(); j++) {
					newProb.conZero[j] = newProb.conZero[j].replaceWithPoly(map);
				}
				for (int j=0; j<newProb.conPositive.size(); j++) {
					newProb.conPositive[j] = newProb.conPositive[j].replaceWithPoly(map);
				}
                for (int j=0; j<newProb.conPSD.size(); j++) {
                    for (int k=0; k<newProb.conPSD[j].size(); k++) {
                        newProb.conPSD[j][k] = newProb.conPSD[j][k].replaceWithPoly(map);
                    }
                }

			}

		}

		return newProb;

	}

	// Get a map going to the minimum amount of variables
	std::unordered_map<int,int> getMinimalMap() {

		// The map from old indices to new (minified) indices
		std::unordered_map<int,int> indMap;
		int nextInd = 0;
		for (int i=0; i<maxVariables; i++) {
			bool foundThis = false;

			// Check the objective
			if (obj.contains(i)) {
				indMap[i] = nextInd;
				foundThis = true;
				nextInd++;
			}

			// Check the equality cons if still not found
			if (!foundThis) {
				for (int j=0; j<conZero.size(); j++) {
					if (conZero[j].contains(i)) {
						indMap[i] = nextInd;
						foundThis = true;
						nextInd++;
						break;
					}
				}
			}

			// Check the inequality cons if still not found
			if (!foundThis) {
				for (int j=0; j<conPositive.size(); j++) {
					if (conPositive[j].contains(i)) {
						indMap[i] = nextInd;
						foundThis = true;
						nextInd++;
						break;
					}
				}
			}

			// Check the PSD cons if still not found
			if (!foundThis) {
				for (int j=0; j<conPSD.size(); j++) {
					for (int k=0; k<conPSD[j].size(); k++) {
						if (conPSD[j][k].contains(i)) {
							indMap[i] = nextInd;
							foundThis = true;
							nextInd++;
							break;
						}
					}
                    if (foundThis) {
                        break;
                    }
				}
			}

		}

		return indMap;
		
	}

    // Check if the problem is already linear (i.e. standard SDP)
    bool isLinear() {
        return getDegree() <= 1;
    }

	// Attempt find a feasible point of this problem
	std::vector<polyType> findFeasibleEqualityPoint(int zeroInd=-1, double alpha=0.9, double tolerance=1e-10, int maxIters=-1, int threads=4, int verbosity=1, double maxMag=1, double stabilityTerm=1e-13, std::vector<polyType> startX={}, double addConstant=0.0) {

		// Combine these to create a single polynomial
		Polynomial<polyType> poly(maxVariables);
		for (int i=0; i<conZero.size(); i++) {
			poly += conZero[i]*conZero[i];
		}
		poly += addConstant;

		// If no index specified, add a var and use that
		if (zeroInd == -1) {
			poly = poly.changeMaxVariables(maxVariables+1);
			zeroInd = poly.maxVariables-1;
		}

		// If verbose, print the polynomial
		if (verbosity >= 2) {
			std::cout << "---------------------" << std::endl;
			std::cout << "As single polynomial:" << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << std::endl;
			std::cout << poly << " = 0" << std::endl;
			std::cout << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << "Optimization:" << std::endl;
			std::cout << "---------------------" << std::endl;
			std::cout << std::endl;
		}

		// Find a root of this polynomial
		return poly.findRoot(zeroInd, alpha, tolerance, maxIters, threads, verbosity, maxMag, stabilityTerm, startX);
		
	}

	// Find the exact solution through brute force 
	std::pair<polyType,std::vector<int>> bruteForce() {

		// If there's no objective, just look for a single valid point
		bool noObj = (obj.size() == 0);

		// Create a vector listing all the variables
		std::vector<int> inds(maxVariables);
		for (int i=0; i<maxVariables; i++) {
			inds[i] = i;
		}

		// Prep fast eval
		obj.prepareEvalFast();
		for (int i=0; i<conZero.size(); i++) {
			conZero[i].prepareEvalFast();
		}
		for (int i=0; i<conPositive.size(); i++) {
			conPositive[i].prepareEvalFast();
		}

		// For each possible set of variables
		long numCombs = std::pow(2, maxVariables);
		polyType bestVal = 10000000;
		std::vector<int> bestSol(maxVariables, 0);
		long outputEvery = std::max(long(1), numCombs / 10000);
		for (long int k=0; k<numCombs; k++) {

			// Convert to -1/1
			std::vector<int> sol(maxVariables, -1);
			for (int i=0; i<maxVariables; i++) {
				if ((k >> i) & 1) {
					sol[maxVariables-i-1] = 1;
				}
			}

			// If it satisfies all constraints
			double conError = 0;
			for (int i=0; i<conZero.size(); i++) {
				conError = std::max(conError, std::abs(conZero[i].evalFast(sol)));
			}
			for (int i=0; i<conPositive.size(); i++) {
				conError = std::max(conError, -std::min(conPositive[i].evalFast(sol), 0.0));
			}

			// If it's valid, check if the objective value is better
			if (conError < zeroTol) {

				// Get the objective value
				double objVal = obj.eval(sol);

				// If there's no objective, finding a feasible point is enough
				if (noObj) {
					std::cout << std::endl;
					return {objVal, sol};
				}

				// Otherwise keep track of the best 
				if (objVal < bestVal) {
					bestVal = objVal;
					bestSol = sol;
				}

			// Progress output every so-many iterations
			} else if (k % outputEvery == 0) {
				std::cout << k << "/" << numCombs << "\r" << std::flush;

			}

		}

		std::cout << std::endl;
		return {bestVal, bestSol};
	
	}

	// Attempt to find a series of constraints that show this is infeasible
	void proveInfeasibleSCS(int maxIters=-1, std::string level="1f", double bound=1, std::string logFileName="", int verbosity=1, int numVarsToSplit=0) {

		// Get the monomial list and sort it
		std::vector<std::string> monoms = getMonomials();
		std::sort(monoms.begin(), monoms.end(), [](const std::string& first, const std::string& second){return first.size() < second.size();});

		// List of semdefinite matrices
		std::vector<std::vector<Polynomial<polyType>>> monomProducts;

		// Get the order of the set of equations
		int degree = getDegree();

		// Make sure we have all first order moments
		addMonomsOfOrder(monoms, 1);

		// First monom should always be 1
		auto loc = std::find(monoms.begin(), monoms.end(), "");
		if (loc != monoms.end()) {
			monoms.erase(loc);
		}
		monoms.insert(monoms.begin(), "");

		// Add all square monoms
		if (degree >= 2) {
			for (int i=0; i<maxVariables; i++) {
				std::string newInd = std::to_string(i);
				newInd.insert(0, digitsPerInd-newInd.size(), ' ');
				if (std::find(monoms.begin(), monoms.end(), newInd+newInd) == monoms.end()) {
					monoms.push_back(newInd+newInd);
				}
			}
		}

		// Add all quartic monoms
		if (degree >= 4) {
			for (int i=0; i<maxVariables; i++) {
				std::string newInd = std::to_string(i);
				newInd.insert(0, digitsPerInd-newInd.size(), ' ');
				if (std::find(monoms.begin(), monoms.end(), newInd+newInd+newInd+newInd) == monoms.end()) {
					monoms.push_back(newInd+newInd+newInd+newInd);
				}
			}
		}
		
		// Create the mapping from monomials to indices (to linearize)
		int numOGMonoms = monoms.size();
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

		// Process the level string e.g. "1+2f,3p"
		std::vector<int> levelsToInclude = {};
		std::string currentThing = "";
		for (int i=0; i<level.size(); i++) {

			// For multiple levels e.g. 1+2
			if (level[i] == '+') {
				levelsToInclude.push_back(std::stoi(currentThing));
				currentThing = "";

			// p for a partial level
			} else if (level[i] == 'p') {

				// Add whatever numbers are left
				levelsToInclude.push_back(std::stoi(currentThing));

				// Get the list monomials that can appear on the top row
				std::vector<std::string> possibleMonoms;
				for (int j=0; j<levelsToInclude.size(); j++) {
					addMonomsOfOrder(possibleMonoms, levelsToInclude[j]);
				}

				// Get all combinations of 3x3 moment matrices
				for (int j=0; j<possibleMonoms.size(); j++) {
					for (int k=j+1; k<possibleMonoms.size(); k++) {
						monomProducts.push_back({Polynomial<polyType>(maxVariables, 1), Polynomial<polyType>(maxVariables, 1, possibleMonoms[j]), Polynomial<polyType>(maxVariables, 1, possibleMonoms[k])});
					}
				}

				// Reset stuff
				currentThing = "";
				levelsToInclude = {};

			// f for a full level
			} else if (level[i] == 'f') {

				// Add whatever numbers are left
				levelsToInclude.push_back(std::stoi(currentThing));

				// Get the list monomials that can appear on the top row
				std::vector<std::string> possibleMonoms;
				for (int j=0; j<levelsToInclude.size(); j++) {
					addMonomsOfOrder(possibleMonoms, levelsToInclude[j]);
				}

				// Add these all to the top row of a moment matrix
				std::vector<Polynomial<double>> toAdd;
				toAdd.push_back(Polynomial<polyType>(maxVariables, 1));
				for (int j=0; j<possibleMonoms.size(); j++) {
					toAdd.push_back(Polynomial<polyType>(maxVariables, 1, possibleMonoms[j]));
				}
				monomProducts.push_back(toAdd);

				// Reset stuff
				currentThing = "";
				levelsToInclude = {};

			// Otherwise add this (probably digit) to the string to process
			} else if (level[i] != ',' && level[i] != ' ') {
				currentThing += level[i];
			}

		}

		// Output for debugging
		if (verbosity >= 2) {
			std::cout << monomProducts << std::endl;
			std::cout << "num in top row of SD mat: " << monomProducts[0].size() << std::endl;
		}

		// Start with the most general area
		double maxArea = 1;
		std::vector<std::vector<std::pair<double,double>>> toProcess;
		std::vector<std::pair<double,double>> varMinMax(obj.maxVariables);
		for (int i=0; i<obj.maxVariables; i++) {
			varMinMax[i].first = -bound;
			varMinMax[i].second = bound;
			maxArea *= 2*bound;
		}
		toProcess.push_back(varMinMax);
				
		// If told to start by splitting at various vars 
		if (numVarsToSplit > 0) {

			// Create the lists of vars to split, starting from the end
			toProcess = {};
			std::vector<int> varsToSplit = {};
			for (int i=maxVariables-1; i>maxVariables-numVarsToSplit-1; i--) {
				varsToSplit.push_back(i);
			}

			// It's exponential in the number of vars, do every combination
			for (int i=0; i<std::pow(2, varsToSplit.size()); i++) {
				auto regionCopy = varMinMax;
				for (int j=0; j<varsToSplit.size(); j++) {
					if (i >> j & 1) {
						regionCopy[varsToSplit[j]].second = 0;
					} else {
						regionCopy[varsToSplit[j]].first = 0;
					}
				}
				toProcess.push_back(regionCopy);
			}

		}

		// Create the PSD matrices from this list
		std::vector<std::vector<int>> shouldBePSD;
		std::vector<std::vector<double>> shouldBePSDCoeffs;
		std::unordered_map<std::string,int> monomsInverted;
		for (int j=0; j<monoms.size(); j++) {
			monomsInverted[monoms[j]] = j;
		}
		for (int j=0; j<monomProducts.size(); j++) {

			// Get the list of all monomial locations for the PSD matrix
			std::vector<int> monLocs;
			std::vector<double> monCoeffs = {};
			for (int i=0; i<monomProducts[j].size(); i++) {
				for (int k=i; k<monomProducts[j].size(); k++) {

					// Calculate the product
					std::string monString = (monomProducts[j][i]*monomProducts[j][k]).getMonomials()[0];

					// Find this in the monomial list
					auto loc = monomsInverted.find(monString);
					if (loc != monomsInverted.end()) {
						monLocs.push_back(monomsInverted[monString]);
					} else {
						monomsInverted[monString] = monoms.size();
						monLocs.push_back(monoms.size());
						monoms.push_back(monString);
					}

					// The coeff for mosek's svec
					if (i != k) {
						monCoeffs.push_back(std::sqrt(2.0));
					} else {
						monCoeffs.push_back(1.0);
					}

				}
			}

			// This (when reformatted) should be positive-semidefinite
			shouldBePSD.push_back(monLocs);
			shouldBePSDCoeffs.push_back(monCoeffs);

		}

		// Get the inds of the first order monomials and their squares
		std::vector<int> firstMonomInds(varMinMax.size(), -1);
		std::vector<std::vector<int>> quadraticMonomInds(varMinMax.size(), std::vector<int>(varMinMax.size(), -1));
		for (int i=0; i<monoms.size(); i++) {
			if (monoms[i].size() == digitsPerInd) {
				firstMonomInds[std::stoi(monoms[i])] = i;
			} else if (monoms[i].size() == 2*digitsPerInd) {
				int ind1 = std::stoi(monoms[i].substr(0,digitsPerInd));
				int ind2 = std::stoi(monoms[i].substr(digitsPerInd,digitsPerInd));
				quadraticMonomInds[ind1][ind2] = i;
				quadraticMonomInds[ind2][ind1] = i;
			}
		}

		// Set some vars
		int oneIndex = 0;
		int varsTotal = monoms.size() + 1;
		int lambdaInd = varsTotal - 1;

		// The box constraints, given our box
		std::vector<double> mins(monoms.size(), -1);
		std::vector<double> maxs(monoms.size(), 1);
		for (int i=0; i<monoms.size(); i++) {

			// Check if the monom is a power of a single var
			bool allSame = true;
			int monomDegree = monoms[i].size() / digitsPerInd;
			for (int j=1; j<monomDegree; j++) {
				if (monoms[i].substr(j*digitsPerInd, digitsPerInd) != monoms[i].substr(0, digitsPerInd)) {
					allSame = false;
					break;
				}
			}

			// If it's an even power, it's positive
			if (allSame && monomDegree % 2 == 0) {
				mins[i] = 0;
			} else {
				mins[i] = -std::pow(bound, monomDegree);
			}
			maxs[i] = std::pow(bound, monomDegree);

		}

		// Convert the linear equality constraints to SCS form
		std::vector<int> ARowsSCS;
		std::vector<int> AColsSCS;
		std::vector<polyType> AValsSCS;
		int nextI = 0;
		ARowsSCS.push_back(nextI);
		AColsSCS.push_back(oneIndex);
		AValsSCS.push_back(1);
		nextI++;
		for (int i=0; i<conZeroLinear.size(); i++) {
			for (auto const &pair: conZeroLinear[i].coeffs) {
				ARowsSCS.push_back(nextI);
				if (pair.first == "") {
					AColsSCS.push_back(oneIndex);
				} else {
					AColsSCS.push_back(std::stoi(pair.first));
				}
				AValsSCS.push_back(pair.second);
			}
			nextI++;
		}

		// The original positivity constraints
		for (int i=0; i<conPositiveLinear.size(); i++) {
			for (auto const &pair: conPositiveLinear[i].coeffs) {
				ARowsSCS.push_back(nextI);
				if (pair.first == "") {
					AColsSCS.push_back(oneIndex);
				} else {
					AColsSCS.push_back(std::stoi(pair.first));
				}
				AValsSCS.push_back(pair.second);
			}
			nextI++;
		}

		// The cutting constraints 
		std::vector<std::vector<int>> paramLocsCutting;
		std::vector<std::vector<int>> paramLocsMin;
		std::vector<std::vector<int>> paramLocsMax;
		for (int i=0; i<toProcess[0].size(); i++) {

			// ax + bx^2 + c >= 0
			ARowsSCS.push_back(nextI);
			AColsSCS.push_back(oneIndex);
			AValsSCS.push_back(42);
			ARowsSCS.push_back(nextI);
			AColsSCS.push_back(firstMonomInds[i]);
			AValsSCS.push_back(42);
			ARowsSCS.push_back(nextI);
			AColsSCS.push_back(quadraticMonomInds[i][i]);
			AValsSCS.push_back(42);
			paramLocsCutting.push_back({int(AValsSCS.size()-3), int(AValsSCS.size()-2), int(AValsSCS.size()-1)});
			nextI++;

		}
		int numParamCons = paramLocsCutting.size();

		// Box constraints
		for (int i=0; i<monoms.size(); i++) {
			ARowsSCS.push_back(nextI);
			AColsSCS.push_back(i);
			AValsSCS.push_back(1);
			nextI++;
		}

		// The SDP constraints
		for (int i=0; i<shouldBePSD.size(); i++) {
			for (int j=0; j<shouldBePSD[i].size(); j++) {
				ARowsSCS.push_back(nextI);
				AColsSCS.push_back(shouldBePSD[i][j]);
				AValsSCS.push_back(shouldBePSDCoeffs[i][j]);
				if (std::abs(shouldBePSDCoeffs[i][j]-1.0) < 1e-5) {
					ARowsSCS.push_back(nextI);
					AColsSCS.push_back(lambdaInd);
					AValsSCS.push_back(1);
				}
				nextI++;
			}
		}

		// Sort the A matrix by columns
		std::vector<int> AOrdering(AValsSCS.size());
		for (int i=0; i<AValsSCS.size(); i++) {
			AOrdering[i] = i;
		}
		std::sort(AOrdering.begin(), AOrdering.end(),
			[&](const int& a, const int& b) {
				return (AColsSCS[a]*ARowsSCS.size() + ARowsSCS[a] < AColsSCS[b]*ARowsSCS.size() + ARowsSCS[b]);
			}
		);

		// Get the inverse map and apply it to the cutting locations
		std::vector<int> AOrderingInv(AOrdering.size());
		for (int i=0; i<AOrdering.size(); i++) {
			AOrderingInv[AOrdering[i]] = i;
		}
		for (int i=0; i<paramLocsCutting.size(); i++) {
			paramLocsCutting[i][0] = AOrderingInv[paramLocsCutting[i][0]];
			paramLocsCutting[i][1] = AOrderingInv[paramLocsCutting[i][1]];
			paramLocsCutting[i][2] = AOrderingInv[paramLocsCutting[i][2]];
		}
		for (int i=0; i<paramLocsMax.size(); i++) {
			paramLocsMax[i][0] = AOrderingInv[paramLocsMax[i][0]];
			paramLocsMax[i][1] = AOrderingInv[paramLocsMax[i][1]];
		}
		for (int i=0; i<paramLocsMin.size(); i++) {
			paramLocsMin[i][0] = AOrderingInv[paramLocsMin[i][0]];
			paramLocsMin[i][1] = AOrderingInv[paramLocsMin[i][1]];
		}

		// The A matrix for SCS 
		double* A_x = new double[AValsSCS.size()];
		for (int i=0; i<AValsSCS.size(); i++) {
			A_x[i] = -AValsSCS[AOrdering[i]];
		}
		int* A_i = new int[ARowsSCS.size()];
		int* A_p = new int[varsTotal+1];
		for (int i=0; i<varsTotal+1; i++) {
			A_p[i] = -1;
		}
		A_p[varsTotal] = AValsSCS.size();
		for (int i=0; i<AOrdering.size(); i++) {
			A_i[i] = ARowsSCS[AOrdering[i]];
			if (A_p[AColsSCS[AOrdering[i]]] == -1) {
				A_p[AColsSCS[AOrdering[i]]] = i;
			}
		}

		// Params needed for SCS
		int* SDPSizes = new int[monomProducts.size()];
		for (int i=0; i<monomProducts.size(); i++) {
			SDPSizes[i] = monomProducts[i].size();
		}
		double* boxConsMin = new double[monoms.size()-1];
		double* boxConsMax = new double[monoms.size()-1];
		for (int i=0; i<monoms.size()-1; i++) {
			boxConsMin[i] = mins[i];
			boxConsMax[i] = maxs[i];
		}
		int numVarsSCS = varsTotal;
		int numConsSCS = nextI;

		// The b vector for SCS
		double* b = new double[numConsSCS];
		for (int i=0; i<numConsSCS; i++) {
			b[i] = 0;
		}
		b[oneIndex] = -1;

		// The c vector for SCS
		double* c = new double[numVarsSCS];
		for (int i=0; i<numVarsSCS; i++) {
			c[i] = 0;
		}
		c[lambdaInd] = 1;

		// Set up the SCS system
		ScsCone *coneSCS = (ScsCone *)calloc(1, sizeof(ScsCone));
		ScsData *dataSCS = (ScsData *)calloc(1, sizeof(ScsData));
		ScsSettings *stgs = (ScsSettings *)calloc(1, sizeof(ScsSettings));
		ScsSolution *sol = (ScsSolution *)calloc(1, sizeof(ScsSolution));
		ScsInfo *info = (ScsInfo *)calloc(1, sizeof(ScsInfo));
		dataSCS->m = numConsSCS;
		dataSCS->n = numVarsSCS;
		dataSCS->b = b;
		dataSCS->c = c;
		dataSCS->P = SCS_NULL;
		auto A = ScsMatrix({A_x, A_i, A_p, dataSCS->m, dataSCS->n});
		dataSCS->A = &A;
		coneSCS->z = 1 + conZeroLinear.size();
		coneSCS->l = numParamCons + conPositiveLinear.size();
		coneSCS->bl = boxConsMin;
		coneSCS->bu = boxConsMax;
		coneSCS->bsize = monoms.size();
		coneSCS->ssize = 1;
		coneSCS->s = SDPSizes;

		// Solver parameters
		scs_set_default_settings(stgs);
		if (verbosity >= 2) {
			stgs->verbose = true;
		} else {
			stgs->verbose = false;
		}
		//stgs->max_iters = 1000;
		//stgs->normalize = false;
		stgs->acceleration_lookback = 0;
		//stgs->adaptive_scale = false;
		//stgs->scale = 0.01;
		//stgs->rho_x = 1e-8;
		stgs->eps_abs = 1e-4;
		stgs->eps_rel = stgs->eps_abs;
		//stgs->eps_infeas = 0.5;

		// Open a file if told to record the points
		std::ofstream logFile;
		if (logFileName.size() > 0) {
			logFile.open(logFileName);
			for (int i=0; i<toProcess[0].size(); i++) {
				logFile << i << ", ";
			}
			logFile << "feasibility";
			logFile << std::endl;
		}

		// Keep splitting until all fail
		int iter = 1;
		iter = 0;
		double totalArea = 0;
		double areaCovered = 1;
		auto begin = std::chrono::steady_clock::now();
		auto end = std::chrono::steady_clock::now();
		double secondsPerIter = 0;
		int numIllPosed = 0;
		while (toProcess.size() > 0) {

			// The area taken up by this section
			areaCovered = 1;
			for (int k=0; k<toProcess[0].size(); k++) {
				areaCovered *= toProcess[0][k].second - toProcess[0][k].first;
			}

			// Update the linear constraints on the quadratics
			int nextInd = 0;
			for (int i=0; i<toProcess[0].size(); i++) {

				// Given two points, find ax+by+c=0
				std::vector<double> point1 = {toProcess[0][i].first, toProcess[0][i].first*toProcess[0][i].first};
				std::vector<double> point2 = {toProcess[0][i].second, toProcess[0][i].second*toProcess[0][i].second};
				std::vector<double> coeffs = getLineFromPoints(point1, point2);

				// Add this as a linear pos con
				dataSCS->A->x[paramLocsCutting[i][0]] = -coeffs[0];
				dataSCS->A->x[paramLocsCutting[i][1]] = -coeffs[1];
				dataSCS->A->x[paramLocsCutting[i][2]] = -coeffs[2];

			}

			// Solve the SDP and then free memory
			ScsWork *scs_work = scs_init(dataSCS, coneSCS, stgs);
			int exitFlag = scs_solve(scs_work, sol, info, 0);
			double objPrimal = info->pobj;
			double objDual = info->dobj;
			double gap = info->gap;
			std::vector<double> solVec(dataSCS->n, 0);
			for (int i=0; i<dataSCS->n; i++) {
				solVec[i] = sol->x[i];
			}
			scs_finish(scs_work);

			// If infeasible, good
			if (verbosity >= 2) {
				std::cout << "    primal: " << objPrimal << ", dual: " << objDual << std::endl;
			}
			if (exitFlag <= 0 || std::min(objDual, objPrimal) > 0.0) {

				// Stop if we ctrl-c'd
				if (exitFlag == -5) {
					return;
				}

				// Write to a file if told to
				if (logFileName.size() > 0) {
					for (int i=0; i<toProcess[0].size(); i++) {
						logFile << solVec[firstMonomInds[i]] << ", ";
					}
					logFile << "no";
					logFile << std::endl;
				}

				// Update the total area count
				totalArea += areaCovered;

			// Otherwise, extract the result and figure out where to split
			} else {

				// Check the resulting vector for a good place to split 
				std::vector<double> errors(maxVariables);
				for (int i=0; i<monoms.size(); i++) {
					if (monoms[i].size() == 2*digitsPerInd) {
						int ind1 = std::stoi(monoms[i].substr(0, digitsPerInd));
						int ind2 = std::stoi(monoms[i].substr(digitsPerInd, digitsPerInd));
						errors[ind1] += std::pow(solVec[i] - solVec[firstMonomInds[ind1]]*solVec[firstMonomInds[ind2]], 2);
						errors[ind2] += std::pow(solVec[i] - solVec[firstMonomInds[ind1]]*solVec[firstMonomInds[ind2]], 2);
					}
				}

				// Find the biggest error
				double biggestError = -10000;
				int bestInd = -1;
				for (int i=0; i<maxVariables; i++) {
					if (errors[i] > biggestError) {
						biggestError = errors[i];
						bestInd = i;
					}
				}

				// If we've converged
				if (biggestError < 1e-6) {

					// If logging, write to file and continue
					if (logFileName.size() > 0) {
						for (int i=0; i<toProcess[0].size(); i++) {
							logFile << solVec[firstMonomInds[i]] << ", ";
						}
						logFile << "yes";
						logFile << std::endl;

					// Otheriwse write to std::out and stop
					} else {
                        if (verbosity >= 1) {
                            std::cout << std::endl;
                            std::cout << "converged in " << iter << " iters to region " << toProcess[0] << std::endl;
                        }
						break;
					}

				// If there's still space to split
				} else {

					// Split it
					double minPoint = toProcess[0][bestInd].first;
					double maxPoint = toProcess[0][bestInd].second;
					double midPoint = (minPoint + maxPoint) / 2.0;
					double mostFeasiblePoint = solVec[firstMonomInds[bestInd]];
					double distanceBetween = mostFeasiblePoint - midPoint;
					double splitPoint = mostFeasiblePoint;
					auto copyLeft = toProcess[0];
					auto copyRight = toProcess[0];
					copyLeft[bestInd].second = splitPoint;
					copyRight[bestInd].first = splitPoint;

					// Add the new paths to the queue
					toProcess.insert(toProcess.begin()+1, copyLeft);
					toProcess.insert(toProcess.begin()+1, copyRight);

				}

			}

			// Time estimation
			end = std::chrono::steady_clock::now();
			secondsPerIter = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / (iter * 1.0e6);
			double areaPerIter = totalArea / iter;
			double itersRemaining = (maxArea - totalArea) / areaPerIter;
			if (totalArea < 1e-8) {
				itersRemaining = -1;
			}
			double secondsRemaining = itersRemaining * secondsPerIter;

			// Per-iteration output
			std::cout << std::defaultfloat;
			if (verbosity >= 2) {
				std::cout << iter << "i  " << 100.0 * totalArea / maxArea << "%  " << 100.0 * areaPerIter / maxArea << "%/i  " << representTime(secondsPerIter) << "/i  " << numIllPosed << "  " << representTime(secondsRemaining) << "  " << 100.0 * areaCovered / maxArea << "%\n" << std::flush;
			} else if (verbosity >= 1) {
				std::cout << iter << "i  " << 100.0 * totalArea / maxArea << "%  " << 100.0 * areaPerIter / maxArea << "%/i  " << representTime(secondsPerIter) << "/i  " << numIllPosed << "  " << representTime(secondsRemaining) << "  " << 100.0 * areaCovered / maxArea << "%                  \r" << std::flush;
			}

			// Remove the one we just processed
			toProcess.erase(toProcess.begin());

			// Keep track of the iteration number
			iter++;
			if (maxIters >= 0 && iter > maxIters) {
				break;
			}	

		}
        if (verbosity >= 1) {
            std::cout << std::endl;
            std::cout << representTime(iter * secondsPerIter) << std::endl;
        }

		// Free the memory used by SCS
		free(coneSCS);
		free(dataSCS);
		free(stgs);
		free(info);
		free(sol->x);
		free(sol->y);
		free(sol->s);
		free(sol);

	}
	
	// Attempt to find a series of constraints that show this is infeasible
	void proveInfeasibleBinary(int maxIters=-1, std::string level="1f", double bound=1, std::string logFileName="", int verbosity=1, int numVarsToSplit=0) {

		// Get the monomial list and sort it
		std::vector<std::string> monoms = getMonomials();
		std::sort(monoms.begin(), monoms.end(), [](const std::string& first, const std::string& second){return first.size() < second.size();});

		// Get the order of the set of equations
		int degree = getDegree();

		// Make sure we have all first order moments
		addMonomsOfOrder(monoms, 1);

		// First monom should always be 1
		auto loc = std::find(monoms.begin(), monoms.end(), "");
		if (loc != monoms.end()) {
			monoms.erase(loc);
		}
		monoms.insert(monoms.begin(), "");

		// Create the mapping from monomials to indices (to linearize)
		int numOGMonoms = monoms.size();
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
		std::vector<std::vector<Polynomial<polyType>>> conPSDLinear;
		for (int i=0; i<conPSD.size(); i++) {
			conPSDLinear.push_back({});
			for (int j=0; j<conPSD[i].size(); j++) {
				conPSDLinear[i].push_back(conPSD[i][j].replaceWithVariable(mapping));
			}
		}

		// Start with the most general area
		std::vector<std::vector<int>> toProcess;
		std::vector<int> startingState(obj.maxVariables, 0);
		toProcess.push_back(startingState);
		double maxArea = std::pow(2, startingState.size());

		// Since it's binary we have linear constraints rather than SDP
		std::vector<std::string> toMaybeAdd;
		if (level == "1f") {
			addMonomsOfOrder(toMaybeAdd, 2);
		} else if (level == "2f") {
			addMonomsOfOrder(toMaybeAdd, 3);
		} else if (level == "3f") {
			addMonomsOfOrder(toMaybeAdd, 4);
		}

		// Don't use monoms with duplicates
		std::vector<std::string> monomsToMakeMats;
		for (int i=0; i<toMaybeAdd.size(); i++) {
			bool valid = true;
			for (int j=0; j<toMaybeAdd[i].size()-digitsPerInd; j+=digitsPerInd) {
				if (toMaybeAdd[i].substr(j, digitsPerInd) == toMaybeAdd[i].substr(j+digitsPerInd, digitsPerInd)) {
					valid = false;
					break;
				}
			}
			if (valid) {
				monomsToMakeMats.push_back(toMaybeAdd[i]);
			}
		}

		// Create a hash table of the monoms for better performance
		std::unordered_map<std::string,int> monomsInverted;
		for (int j=0; j<monoms.size(); j++) {
			monomsInverted[monoms[j]] = j;
		}

		// Create the linear cons from this list
		for (int j=0; j<monomsToMakeMats.size(); j++) {

			// For each combination of the components of this monom
			std::vector<std::string> baseComponents;
			for (int i=0; i<monomsToMakeMats[j].size(); i+=digitsPerInd) {
				baseComponents.push_back(monomsToMakeMats[j].substr(i, digitsPerInd));
			}

			// Now add the various combinations of these base components
			std::vector<std::string> monomialsToAffect = baseComponents;
			monomialsToAffect.insert(monomialsToAffect.begin(), "");
			if (monomsToMakeMats[j].size() >= 2*digitsPerInd) {
				for (int i1=0; i1<baseComponents.size(); i1++) {
					for (int i2=i1+1; i2<baseComponents.size(); i2++) {
						monomialsToAffect.push_back(baseComponents[i1] + baseComponents[i2]);
					}
				}
			}
			if (monomsToMakeMats[j].size() >= 3*digitsPerInd) {
				for (int i1=0; i1<baseComponents.size(); i1++) {
					for (int i2=i1+1; i2<baseComponents.size(); i2++) {
						for (int i3=i2+1; i3<baseComponents.size(); i3++) {
							monomialsToAffect.push_back(baseComponents[i1] + baseComponents[i2] + baseComponents[i3]);
						}
					}
				}
			}
			if (monomsToMakeMats[j].size() >= 4*digitsPerInd) {
				for (int i1=0; i1<baseComponents.size(); i1++) {
					for (int i2=i1+1; i2<baseComponents.size(); i2++) {
						for (int i3=i2+1; i3<baseComponents.size(); i3++) {
							for (int i4=i3+1; i4<baseComponents.size(); i4++) {
								monomialsToAffect.push_back(baseComponents[i1] + baseComponents[i2] + baseComponents[i3] + baseComponents[i4]);
							}
						}
					}
				}
			}

			// Make sure these are all present in the monoms list
			for (int k=0; k<monomialsToAffect.size(); k++) {
				if (monomsInverted.find(monomialsToAffect[k]) == monomsInverted.end()) {
					monomsInverted[monomialsToAffect[k]] = monoms.size();
					monoms.push_back(monomialsToAffect[k]);
				}
			}

			// Iterate over all of the base components
			int numCombs = std::pow(2, baseComponents.size());
			for (int k=0; k<numCombs; k++) {

				// The base components should iterate over -1/+1
				std::vector<int> coeffs(monomialsToAffect.size(), -1);
				coeffs[0] = 1;
				for (int i=0; i<baseComponents.size(); i++) {
					if ((k >> i) & 1) {
						coeffs[1+baseComponents.size()-i-1] = 1;
					}
				}

				// Then everything else should be formed based on these
				for (int i=1+baseComponents.size(); i<monomialsToAffect.size(); i++) {
					coeffs[i] = 1;
					for (int l=0; l<monomialsToAffect[i].size(); l+=digitsPerInd) {
						std::string indString = monomialsToAffect[i].substr(l, digitsPerInd);
						int loc = std::find(monomialsToAffect.begin(), monomialsToAffect.end(), indString) - monomialsToAffect.begin();
						coeffs[i] *= coeffs[loc];
					}
				}

				// Add this constraint, pre-linearized
				Polynomial<polyType> newCon(monoms.size());
				for (int i=0; i<monomialsToAffect.size(); i++) {
					newCon.addTerm(coeffs[i], {monomsInverted[monomialsToAffect[i]]});
				}
				conPositiveLinear.push_back(newCon);

			}

		}

		// Get the inds of the first order monomials and their squares
		std::vector<int> firstMonomInds(startingState.size(), -1);
		std::vector<std::vector<int>> quadraticMonomInds(startingState.size(), std::vector<int>(startingState.size(), -1));
		for (int i=0; i<monoms.size(); i++) {
			if (monoms[i].size() == digitsPerInd) {
				firstMonomInds[std::stoi(monoms[i])] = i;
			} else if (monoms[i].size() == 2*digitsPerInd) {
				int ind1 = std::stoi(monoms[i].substr(0, digitsPerInd));
				int ind2 = std::stoi(monoms[i].substr(digitsPerInd, digitsPerInd));
				quadraticMonomInds[ind1][ind2] = i;
				quadraticMonomInds[ind2][ind1] = i;
			}
		}

		// Set some vars
		int oneIndex = 0;
		int varsTotal = monoms.size();

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

		// Determine the largest amount of terms per element of the PSD matrix
		int numMatsToSum = 1;
		for (int i=0; i<conPSDLinear.size(); i++) {
			for (int j=i; j<conPSDLinear[i].size(); j++) {
				numMatsToSum = std::max(numMatsToSum, int(conPSDLinear[i][j].coeffs.size()));
			}
		}

		// Convert the linear PSD constraints to MOSEK form
		int upperTriangSize = (conPSDLinear.size()*(conPSDLinear.size()+1)) / 2;
		std::vector<std::vector<int>> psdLocs(numMatsToSum, std::vector<int>(upperTriangSize, oneIndex));
		std::vector<std::vector<double>> psdCoeffs(numMatsToSum, std::vector<double>(upperTriangSize, 0.0));
		int whichEl = 0;
		for (int i=0; i<conPSDLinear.size(); i++) {
			for (int j=i; j<conPSDLinear[i].size(); j++) {
				int whichMat = 0;
				if (conPSDLinear[i][j].coeffs.size() > 0) {
					for (auto const &pair: conPSDLinear[i][j].coeffs) {
						if (pair.first == "") {
							psdLocs[whichMat][whichEl] = oneIndex;
						} else {
							psdLocs[whichMat][whichEl] = std::stoi(pair.first);
						}
						if (i == j) {
							psdCoeffs[whichMat][whichEl] = pair.second;
						} else {
							psdCoeffs[whichMat][whichEl] = std::sqrt(2.0)*pair.second;
						}
						whichMat++;
					}
				} else {
					psdLocs[whichMat][whichEl] = oneIndex;
					psdCoeffs[whichMat][whichEl] = 0;
				}
				whichEl++;
			}
		}
		std::vector<std::shared_ptr<monty::ndarray<int,1>>> psdLocsM(numMatsToSum);
		std::vector<std::shared_ptr<monty::ndarray<double,1>>> psdCoeffsM(numMatsToSum);
		for (int i=0; i<numMatsToSum; i++) {
			psdLocsM[i] = monty::new_array_ptr<int>(psdLocs[i]);
			psdCoeffsM[i] = monty::new_array_ptr<double>(psdCoeffs[i]);
		}

		// Create a model
		mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});

		// Create the variable
		mosek::fusion::Variable::t xM = M->variable(varsTotal, mosek::fusion::Domain::inRange(-1, 1));

		// If asking for everything
		if (verbosity >= 3) {
			M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
		}

		// The first element of the vector should be one
		M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));

		// Linear equality constraints
		M->constraint(mosek::fusion::Expr::mul(AM, xM), mosek::fusion::Domain::equalsTo(0.0));

		// Linear positivity constraints
		M->constraint(mosek::fusion::Expr::mul(BM, xM), mosek::fusion::Domain::greaterThan(0.0));

		// Seems to help
		M->objective(mosek::fusion::ObjectiveSense::Maximize, mosek::fusion::Expr::sum(xM));

		// Linear PSD constraints
		if (conPSDLinear.size() > 0) {
			mosek::fusion::Expression::t matSum = mosek::fusion::Expr::mulElm(psdCoeffsM[0], xM->pick(psdLocsM[0]));
			for (int i=1; i<numMatsToSum; i++) {
				matSum = mosek::fusion::Expr::add(matSum, mosek::fusion::Expr::mulElm(psdCoeffsM[i], xM->pick(psdLocsM[i])));
			}
			M->constraint(matSum, mosek::fusion::Domain::inSVecPSDCone());
		}

		// Open a file if told to record the points
		std::ofstream logFile;
		if (logFileName.size() > 0) {
			logFile.open(logFileName);
			for (int i=0; i<toProcess[0].size(); i++) {
				logFile << i << ", ";
			}
			logFile << "feasibility";
			logFile << std::endl;
		}

		// Keep splitting until all fail
		int iter = 1;
		double totalArea = 0;
		double areaCovered = 1;
		auto begin = std::chrono::steady_clock::now();
		auto end = std::chrono::steady_clock::now();
		double secondsPerIter = 0;
		int numIllPosed = 0;
		std::vector<monty::rc_ptr<mosek::fusion::Constraint>> fixingCons;
		while (toProcess.size() > 0) {

			// Debug output
			if (verbosity >= 2) {
				std::cout << "    checking region: " << toProcess[0] << std::endl;
			}

			// The area taken up by this section
			areaCovered = 1;
			for (int k=0; k<toProcess[0].size(); k++) {
				if (toProcess[0][k] == 0) {
					areaCovered *= 2;
				}
			}

			// The new parameterized constraint vector and objective
			//std::vector<std::vector<double>> newD(numParamCons, std::vector<double>(varsTotal, 0));

			// Update the linear constraints fixing certain variables
			int nextInd = 0;
			for (int i=0; i<fixingCons.size(); i++) {
				fixingCons[i]->remove();
			}
			fixingCons = {};
			for (int i=0; i<toProcess[0].size(); i++) {
				if (toProcess[0][i] != 0) {
					fixingCons.push_back(M->constraint(xM->index(firstMonomInds[i]), mosek::fusion::Domain::equalsTo(toProcess[0][i])));
				}
			}
			for (int i=0; i<toProcess[0].size(); i++) {
				for (int j=i; j<toProcess[0].size(); j++) {

					// Only constrain terms we actually reference
					if (quadraticMonomInds[i][j] < 0) {
						continue;
					}

					// The square terms are one
					if (i == j) {
						fixingCons.push_back(M->constraint(xM->index(quadraticMonomInds[i][j]), mosek::fusion::Domain::equalsTo(1.0)));

					// If fixing one variable, then the quadratic is equal to the non-fixed variable
					} else if (toProcess[0][i] == 1 && toProcess[0][j] == 0) {
						fixingCons.push_back(M->constraint(mosek::fusion::Expr::sub(xM->index(quadraticMonomInds[i][j]), xM->index(firstMonomInds[j])), mosek::fusion::Domain::equalsTo(0.0)));
					} else if (toProcess[0][i] == -1 && toProcess[0][j] == 0) {
						fixingCons.push_back(M->constraint(mosek::fusion::Expr::add(xM->index(quadraticMonomInds[i][j]), xM->index(firstMonomInds[j])), mosek::fusion::Domain::equalsTo(0.0)));

					// If fixing one variable, then the quadratic is equal to the non-fixed variable
					} else if (toProcess[0][i] == 0 && toProcess[0][j] == 1) {
						fixingCons.push_back(M->constraint(mosek::fusion::Expr::sub(xM->index(quadraticMonomInds[i][j]), xM->index(firstMonomInds[i])), mosek::fusion::Domain::equalsTo(0.0)));
					} else if (toProcess[0][i] == 0 && toProcess[0][j] == -1) {
						fixingCons.push_back(M->constraint(mosek::fusion::Expr::add(xM->index(quadraticMonomInds[i][j]), xM->index(firstMonomInds[i])), mosek::fusion::Domain::equalsTo(0.0)));

					// If fixing both variables, then the quadratic is equal to the product
					} else if (toProcess[0][i] != 0 && toProcess[0][j] != 0) {
						fixingCons.push_back(M->constraint(xM->index(quadraticMonomInds[i][j]), mosek::fusion::Domain::equalsTo(toProcess[0][i]*toProcess[0][j])));
					}

				}
			}

			// Solve the problem
			M->solve();
			auto statProb = M->getProblemStatus();
			auto statSol = M->getPrimalSolutionStatus();

			// If infeasible, good
			if (statProb == mosek::fusion::ProblemStatus::PrimalInfeasible) {

				// Write to a file if told to
				if (logFileName.size() > 0 && statSol != mosek::fusion::SolutionStatus::Undefined && statSol != mosek::fusion::SolutionStatus::Unknown && statProb != mosek::fusion::ProblemStatus::PrimalInfeasible) {
					auto sol = *(xM->level());
					std::vector<polyType> solVec(xM->getSize());
					for (int i=0; i<solVec.size(); i++) {
						solVec[i] = sol[i];
					}
					for (int i=0; i<toProcess[0].size(); i++) {
						logFile << solVec[firstMonomInds[i]] << ", ";
					}
					logFile << "no";
					logFile << std::endl;
				}

				// Update the total area count
				totalArea += areaCovered;

			// Otherwise, extract the result and figure out where to split
			} else {

				// Get the solution values
				auto sol = *(xM->level());
				polyType objPrimal = M->primalObjValue();
				polyType objDual = M->dualObjValue();

				// Output the relevant moments
				std::vector<polyType> solVec(xM->getSize());
				for (int i=0; i<solVec.size(); i++) {
					solVec[i] = sol[i];
				}

				// Output each monom along with it's solution value
				if (verbosity >= 2) {
					for (int i=0; i<monoms.size(); i++) {
						std::cout << monoms[i] << ": " << solVec[i] << std::endl;
					}
				}

				// Check the resulting vector for a good place to split
				std::vector<double> errors(maxVariables);
				for (int i=0; i<monoms.size(); i++) {
					if (monoms[i].size() == 2*digitsPerInd) {
						int ind1 = std::stoi(monoms[i].substr(0, digitsPerInd));
						int ind2 = std::stoi(monoms[i].substr(digitsPerInd, digitsPerInd));
						errors[ind1] += std::abs(solVec[i] - solVec[firstMonomInds[ind1]]*solVec[firstMonomInds[ind2]]);
						errors[ind2] += std::abs(solVec[i] - solVec[firstMonomInds[ind1]]*solVec[firstMonomInds[ind2]]);
					}
				}

				// Find the biggest error
				double biggestError = -10000;
				double smallestError = 10000;
				int bestInd = -1;
				for (int i=0; i<maxVariables; i++) {
					if (errors[i] > biggestError && toProcess[0][i] == 0) {
						biggestError = errors[i];
						bestInd = i;
					}
					if (errors[i] < smallestError && toProcess[0][i] == 0) {
						smallestError = errors[i];
					}
				}

				// Output the largest error
				if (verbosity >= 2) {
					std::cout << "    biggest error: " << biggestError << " at index " << bestInd << std::endl;
				}

				// If we've converged
				if (biggestError < 1e-8) {

					// If logging, write to file and continue
					if (logFileName.size() > 0) {
						for (int i=0; i<toProcess[0].size(); i++) {
							logFile << solVec[firstMonomInds[i]] << ", ";
						}
						logFile << "yes";
						logFile << std::endl;

					// Otheriwse write to std::out and stop
					} else {
						std::cout << std::endl;
						std::vector<double> vals(maxVariables);
						for (int i=0; i<maxVariables; i++) {
							vals[i] = solVec[firstMonomInds[i]];
						}
						std::cout << "converged in " << iter << " iters to: " << vals << std::endl;
						break;
					}

				// If there's still space to split
				} else {

					// Split it 
					if (verbosity >= 2) {
						std::cout << "    splitting var " << bestInd << std::endl;
					}
					auto copyLeft = toProcess[0];
					auto copyRight = toProcess[0];
					copyLeft[bestInd] = 1;
					copyRight[bestInd] = -1;

					// Add the new paths to the queue
					toProcess.insert(toProcess.begin()+1, copyLeft);
					toProcess.insert(toProcess.begin()+1, copyRight);

				}

			}

			// Time estimation
			end = std::chrono::steady_clock::now();
			secondsPerIter = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / (iter * 1.0e6);
			double areaPerIter = totalArea / iter;
			double itersRemaining = (maxArea - totalArea) / areaPerIter;
			double secondsRemaining = itersRemaining * secondsPerIter;

			// Per-iteration output
			std::cout << std::defaultfloat;
			if (verbosity >= 2) {
				std::cout << iter << "i  " << 100.0 * totalArea / maxArea << "%  " << 100.0 * areaPerIter / maxArea << "%/i  " << representTime(secondsPerIter) << "/i  " << numIllPosed << "  " << representTime(secondsRemaining) << "  " << 100.0 * areaCovered / maxArea << "%\n" << std::flush;
			} else if (verbosity >= 1) {
				std::cout << iter << "i  " << 100.0 * totalArea / maxArea << "%  " << 100.0 * areaPerIter / maxArea << "%/i  " << representTime(secondsPerIter) << "/i  " << numIllPosed << "  " << representTime(secondsRemaining) << "  " << 100.0 * areaCovered / maxArea << "%                  \r" << std::flush;
			}

			// Remove the one we just processed
			toProcess.erase(toProcess.begin());

			// Keep track of the iteration number
			iter++;
			if (maxIters >= 0 && iter > maxIters) {
				break;
			}	

		}
		if (verbosity >= 1) { 
			std::cout << std::endl;
			std::cout << representTime(iter * secondsPerIter) << std::endl;
		}

	}

	// Given a point, check if it meets all of the constraints
	bool isFeasible(std::vector<polyType> x, double customTol=1e-7) {

		// Equality constraints
		for (int i=0; i<conZero.size(); i++) {
            polyType res = conZero[i].eval(x);
			if (std::abs(res) > customTol) {
				std::cout << "equality constraint " << i << " violated: " << res << std::endl;
				return false;
			}
		}

		// Inequality constraints
		for (int i=0; i<conPositive.size(); i++) {
            polyType res = conPositive[i].eval(x);
			if (res < -customTol) {
				std::cout << "positive constraint " << i << " violated: " << res << std::endl;
				return false;
			}
		}

		// Variable bounds
		for (int i=0; i<x.size(); i++) {
			if (x[i] < varBounds[i].first - customTol || x[i] > varBounds[i].second + customTol) {
				std::cout << "variable bound " << i << " violated: " << x[i] << std::endl;
				return false;
			}
		}

		// Binary constraints
		for (int i=0; i<x.size(); i++) {
			if (varIsBinary[i] && std::abs(x[i] - varBounds[i].first) > customTol && std::abs(x[i] - varBounds[i].second) > customTol) {
				std::cout << "binary constraint " << i << " violated: " << x[i] << std::endl;
				return false;
			}
		}

		// Positive semi-definiteness constraints
		if (conPSD.size() > 0) {
			Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(conPSD.size(), conPSD[0].size());
			for (int i=0; i<conPSD.size(); i++) {
				for (int j=0; j<conPSD[0].size(); j++) {
					mat(i,j) = conPSD[i][j].eval(x);
				}
			}
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mat);
            double smallestEigenvalue = es.eigenvalues().minCoeff();
			if (smallestEigenvalue < -customTol) {
                std::cout << "PSD constraint violated: " << smallestEigenvalue << std::endl;
				return false;
			}
		}

		// Otherwise I guess it's fine
		return true;

	}

	// Given a point, find the nearest feasible point
	std::vector<polyType> getNearestFeasiblePoint(std::vector<polyType> x, std::vector<std::pair<polyType, polyType>> region) {
	
		// Start with a copy of the point
		std::vector<polyType> xFeasible = x;

		// For each variable
		for (int i=0; i<x.size(); i++) {

			// If it's a binary variable, we just need to round
			if (varIsBinary[i]) {
				double midPoint = (varBounds[i].first + varBounds[i].second) / 2.0;
				if (x[i] < midPoint) {
					xFeasible[i] = varBounds[i].first;
				} else {
					xFeasible[i] = varBounds[i].second;
				}

			// Otherwise for now just use the point
			} else {
				xFeasible[i] = x[i];

			}

		}

		// If it's feasible, return it
		if (isFeasible(xFeasible)) {
			return xFeasible;
		} else {
			return std::vector<polyType>();
		}

	}

    // Get a bunch of points in the feasible region
    std::vector<std::vector<polyType>> manyFeasible(int num, int verbosity=1) {

        // Save the objective
        Polynomial<polyType> prevObj = obj;

        // For the number that we need
        std::vector<std::vector<polyType>> points;
        for (int i=0; i<num; i++) {
            if (verbosity >= 1) {
                std::cout << i << " / " << num << "\r" << std::flush;
            }

            // Random objective
            obj = Polynomial<polyType>(maxVariables);
            for (int j=0; j<maxVariables; j++) {
                obj.addTerm(double(rand()) / double(RAND_MAX), {j});
            }

            // Solve and get the solution
            auto sol = optimize(1, 0, 1);
            points.push_back(sol.second);

        }
        if (verbosity >= 1) {
            std::cout << std::endl;
        }

        // Restore the previous objective
        obj = prevObj;

        return points;

    }

    // Get a bunch of points near the edge
    std::vector<std::vector<polyType>> manyFeasibleMixed(int num1, int num2, int verbosity=1) {

        // Generate a bunch of random feasible points
        std::vector<std::vector<polyType>> feasiblePoints = manyFeasible(num1, verbosity);

        // For the number that we need
        std::vector<std::vector<polyType>> points;
        for (int i=0; i<num2; i++) {

            // Start from a random feasible point
            int feasInd = rand() % num1;
            std::vector<polyType> x = feasiblePoints[feasInd];

            // Pick a random point to mix with
            int randInd = rand() % num1;

            // Get a random combination of the feasible points
            double t = double(rand()) / double(RAND_MAX);
            for (int k=0; k<maxVariables; k++) {
                x[k] = (1-t) * x[k] + t * feasiblePoints[randInd][k];
            }

            // Add this to the list
            points.push_back(x);

        }

        return points;

    }

	// Minimize using branch and bound plus SDP
	std::pair<polyType, std::vector<polyType>> optimize(int level=1, int verbosity=1, int maxIters=-1) {

        // If it's linear, one iteration is enough
        bool isLin = isLinear();
        if (isLin) {
            maxIters = 1;
        }

		// Get the monomial list and sort it
		std::vector<std::string> monoms = getMonomials();
		std::sort(monoms.begin(), monoms.end(), [](const std::string& first, const std::string& second){return first.size() < second.size();});

		// List of semdefinite matrices
		std::vector<std::vector<Polynomial<polyType>>> monomProducts;

		// Get the order of the set of equations
		int degree = getDegree();

		// Make sure we have all first order moments
		addMonomsOfOrder(monoms, 1);

		// First monom should always be 1
		auto loc = std::find(monoms.begin(), monoms.end(), "");
		if (loc != monoms.end()) {
			monoms.erase(loc);
		}
		monoms.insert(monoms.begin(), "");

		// Add all square monoms
		if (degree >= 2) {
			for (int i=0; i<maxVariables; i++) {
				std::string newInd = std::to_string(i);
				newInd.insert(0, digitsPerInd-newInd.size(), ' ');
				if (std::find(monoms.begin(), monoms.end(), newInd+newInd) == monoms.end()) {
					monoms.push_back(newInd+newInd);
				}
			}
		}

		// Add all quartic monoms
		if (degree >= 4) {
			for (int i=0; i<maxVariables; i++) {
				std::string newInd = std::to_string(i);
				newInd.insert(0, digitsPerInd-newInd.size(), ' ');
				if (std::find(monoms.begin(), monoms.end(), newInd+newInd+newInd+newInd) == monoms.end()) {
					monoms.push_back(newInd+newInd+newInd+newInd);
				}
			}
		}
		
		// Create the mapping from monomials to indices (to linearize)
		int numOGMonoms = monoms.size();
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
		std::vector<std::vector<Polynomial<polyType>>> conPSDLinear;
		for (int i=0; i<conPSD.size(); i++) {
			conPSDLinear.push_back({});
			for (int j=0; j<conPSD[i].size(); j++) {
				conPSDLinear[i].push_back(conPSD[i][j].replaceWithVariable(mapping));
			}
		}

        // If we're doing a moment based approximation
        if (!isLin) {

            // Get the list monomials that can appear on the top row
            std::vector<std::string> possibleMonoms;
            for (int j=0; j<level; j++) {
                addMonomsOfOrder(possibleMonoms, j+1);
            }

            // Make sure that your terms in your objective appear somewhere TODO
            for (int j=0; j<obj.size(); j++) {
                std::string monom = obj.getMonomials()[j];
                monom = monom.substr(0, monom.size()-digitsPerInd);
                if (std::find(possibleMonoms.begin(), possibleMonoms.end(), monom) == possibleMonoms.end()) {
                    if (verbosity >= 2) {
                        std::cout << "Adding " << monom << " to the moment matrix" << std::endl;
                    }
                    possibleMonoms.push_back(monom);
                }
            }

            // Only allow certain terms in the moment matrix
            std::vector<std::string> allowedMonoms = {" 0", " 6", "15"};
            if (allowedMonoms.size() > 0) {
                for (int j=0; j<possibleMonoms.size(); j++) {
                    bool containsAllowed = false;
                    for (int k=0; k<allowedMonoms.size(); k++) {
                        if (possibleMonoms[j].find(allowedMonoms[k]) != std::string::npos) {
                            containsAllowed = true;
                            break;
                        }
                    }
                    if (!containsAllowed) {
                        if (verbosity >= 2) {
                            std::cout << "Removing " << possibleMonoms[j] << " from the moment matrix" << std::endl;
                        }
                        possibleMonoms.erase(possibleMonoms.begin()+j);
                        j--;
                    }
                }
            }

            // Add these all to the top row of a moment matrix
            std::vector<Polynomial<double>> toAdd;
            if (possibleMonoms.size() > 0) {
                toAdd.push_back(Polynomial<polyType>(maxVariables, 1));
                for (int j=0; j<possibleMonoms.size(); j++) {
                    toAdd.push_back(Polynomial<polyType>(maxVariables, 1, possibleMonoms[j]));
                }
                monomProducts.push_back(toAdd);
            }

        }

		// Determine the outer bound for all variables
		std::pair<double,double> generalBounds = {0,0};
		for (int i=0; i<obj.maxVariables; i++) {
			if (varBounds[i].second > generalBounds.second) {
				generalBounds.second = varBounds[i].second;
			}
			if (varBounds[i].first < generalBounds.first) {
				generalBounds.first = varBounds[i].first;
			}
		}

		// Start with the most general area
		double maxArea = 1;
		std::vector<std::vector<std::pair<double,double>>> toProcess;
		for (int i=0; i<obj.maxVariables; i++) {
			if (varIsBinary[i]) {
				maxArea *= 2;
			} else {
				maxArea *= varBounds[i].second-varBounds[i].first;
			}
		}
		toProcess.push_back(varBounds);

		// Verbose output
		if (verbosity >= 2 && monomProducts.size() > 0) {
			std::cout << "Moment matrices: " << std::endl;
			for (int i=0; i<monomProducts.size(); i++) {
				std::cout << monomProducts[i] << std::endl;
			}
			std::cout << std::endl;
		}

		// Create the PSD matrices from this list
		std::vector<std::shared_ptr<monty::ndarray<int,1>>> shouldBePSD;
		std::vector<std::shared_ptr<monty::ndarray<double,1>>> shouldBePSDCoeffs;
		std::vector<std::shared_ptr<monty::ndarray<double,1>>> identityAsSVec;
		std::vector<double> idenCoeffs;
		std::vector<double> monCoeffs;
		std::vector<int> monLocs;
		std::unordered_map<std::string,int> monomsInverted;
		for (int j=0; j<monoms.size(); j++) {
			monomsInverted[monoms[j]] = j;
		}
		for (int j=0; j<monomProducts.size(); j++) {

			// Get the list of all monomial locations for the PSD matrix
			monLocs = {};
			monCoeffs = {};
			idenCoeffs = {};
			for (int i=0; i<monomProducts[j].size(); i++) {
				for (int k=i; k<monomProducts[j].size(); k++) {

					// Calculate the product
					std::string monString = (monomProducts[j][i]*monomProducts[j][k]).getMonomials()[0];

					// Find this in the monomial list
					auto loc = monomsInverted.find(monString);
					if (loc != monomsInverted.end()) {
						monLocs.push_back(monomsInverted[monString]);
					} else {
						monomsInverted[monString] = monoms.size();
						monLocs.push_back(monoms.size());
						monoms.push_back(monString);
					}

					// The coeff for mosek's svec
					if (i != k) {
						monCoeffs.push_back(std::sqrt(2.0));
						idenCoeffs.push_back(0.0);
					} else {
						monCoeffs.push_back(1.0);
						idenCoeffs.push_back(1.0);
					}

				}
			}

			// This (when reformatted) should be positive-semidefinite
			shouldBePSD.push_back(monty::new_array_ptr<int>(monLocs));
			shouldBePSDCoeffs.push_back(monty::new_array_ptr<double>(monCoeffs));
			identityAsSVec.push_back(monty::new_array_ptr<double>(idenCoeffs));

			// Clear some memory
			monomProducts.erase(monomProducts.begin());
			j--;

		}

		// Get the inds of the first order monomials and their squares
		std::vector<int> firstMonomInds(maxVariables, -1);
		std::vector<std::vector<int>> quadraticMonomInds(maxVariables, std::vector<int>(maxVariables, -1));
		for (int i=0; i<monoms.size(); i++) {
			if (monoms[i].size() == digitsPerInd) {
				firstMonomInds[std::stoi(monoms[i])] = i;
			} else if (monoms[i].size() == 2*digitsPerInd) {
				int ind1 = std::stoi(monoms[i].substr(0,digitsPerInd));
				int ind2 = std::stoi(monoms[i].substr(digitsPerInd,digitsPerInd));
				quadraticMonomInds[ind1][ind2] = i;
				quadraticMonomInds[ind2][ind1] = i;
			}
		}

		// Set some vars
		int oneIndex = 0;
		int varsTotal = monoms.size();

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

		// Convert the linear objective function to MOSEK form
		std::vector<double> objCoeffs(varsTotal, 0.0);
		for (auto const &pair: objLinear.coeffs) {
			if (pair.first == "") {
				objCoeffs[oneIndex] = pair.second;
			} else {
				objCoeffs[std::stoi(pair.first)] = pair.second;
			}
		}
		auto cM = mosek::fusion::Matrix::dense(1, varsTotal, monty::new_array_ptr<double>(objCoeffs));

		// Determine the largest amount of terms per element of the PSD matrix
		int numMatsToSum = 1;
		for (int i=0; i<conPSDLinear.size(); i++) {
			for (int j=i; j<conPSDLinear[i].size(); j++) {
				numMatsToSum = std::max(numMatsToSum, int(conPSDLinear[i][j].coeffs.size()));
			}
		}

		// Convert the linear PSD constraints to MOSEK form
		int upperTriangSize = (conPSDLinear.size()*(conPSDLinear.size()+1)) / 2;
		std::vector<std::vector<int>> psdLocs(numMatsToSum, std::vector<int>(upperTriangSize, oneIndex));
		std::vector<std::vector<double>> psdCoeffs(numMatsToSum, std::vector<double>(upperTriangSize, 0.0));
		int whichEl = 0;
		for (int i=0; i<conPSDLinear.size(); i++) {
			for (int j=i; j<conPSDLinear[i].size(); j++) {
				int whichMat = 0;
				if (conPSDLinear[i][j].coeffs.size() > 0) {
					for (auto const &pair: conPSDLinear[i][j].coeffs) {
						if (pair.first == "") {
							psdLocs[whichMat][whichEl] = oneIndex;
						} else {
							psdLocs[whichMat][whichEl] = std::stoi(pair.first);
						}
						if (i == j) {
							psdCoeffs[whichMat][whichEl] = pair.second;
						} else {
							psdCoeffs[whichMat][whichEl] = std::sqrt(2.0)*pair.second;
						}
						whichMat++;
					}
				} else {
					psdLocs[whichMat][whichEl] = oneIndex;
					psdCoeffs[whichMat][whichEl] = 0;
				}
				whichEl++;
			}
		}
		std::vector<std::shared_ptr<monty::ndarray<int,1>>> psdLocsM(numMatsToSum);
		std::vector<std::shared_ptr<monty::ndarray<double,1>>> psdCoeffsM(numMatsToSum);
		for (int i=0; i<numMatsToSum; i++) {
			psdLocsM[i] = monty::new_array_ptr<int>(psdLocs[i]);
			psdCoeffsM[i] = monty::new_array_ptr<double>(psdCoeffs[i]);
		}

		// Create a model
		mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});

		// Create the variable
		mosek::fusion::Variable::t xM = M->variable(varsTotal, mosek::fusion::Domain::inRange(generalBounds.first, generalBounds.second));

		// If the variables should be constrained
		for (int i=0; i<maxVariables; i++) {

			// If there are bounds
			if (varBounds[i].first > -largeTol && varBounds[i].second < largeTol) {
				M->constraint(xM->index(firstMonomInds[i]), mosek::fusion::Domain::inRange(varBounds[i].first, varBounds[i].second));
			} else if (varBounds[i].first > -largeTol) {
				M->constraint(xM->index(firstMonomInds[i]), mosek::fusion::Domain::greaterThan(varBounds[i].first));
			} else if (varBounds[i].second < largeTol) {
				M->constraint(xM->index(firstMonomInds[i]), mosek::fusion::Domain::lessThan(varBounds[i].second));
			}

			// If it's binary, either -1/1 or 0/1
			if (varIsBinary[i]) {
				if (varBounds[i].first == -1 && varBounds[i].second == 1) {
					M->constraint(xM->index(quadraticMonomInds[i][i]), mosek::fusion::Domain::equalsTo(1.0));
				} else if (varBounds[i].first == 0 && varBounds[i].second == 1) {
					M->constraint(mosek::fusion::Expr::sub(xM->index(quadraticMonomInds[i][i]), xM->index(firstMonomInds[i])), mosek::fusion::Domain::equalsTo(0.0));
				}
			}

		}

		// Parameterized linear positivity constraints
		int numParamPosCons = 0;
		std::vector<long> sparsity;
		for (int i=0; i<toProcess[0].size(); i++) {
            if (quadraticMonomInds[i][i] != -1) {
                sparsity.push_back(i*varsTotal + firstMonomInds[i]);
                sparsity.push_back(i*varsTotal + quadraticMonomInds[i][i]);
                sparsity.push_back(i*varsTotal + oneIndex);
                numParamPosCons++;
            }
		}
		std::sort(sparsity.begin(), sparsity.end());
        mosek::fusion::Parameter::t DM = M->parameter(monty::new_array_ptr<int>({numParamPosCons, varsTotal}), monty::new_array_ptr<long>(sparsity));
        M->constraint(mosek::fusion::Expr::mul(DM, xM), mosek::fusion::Domain::greaterThan(0));

		// Parameterized linear equality constraints
		int numParamEqCons = 0;
		std::vector<long> sparsity2;
		for (int i=0; i<toProcess[0].size(); i++) {
			sparsity2.push_back(i*varsTotal + firstMonomInds[i]);
			sparsity2.push_back(i*varsTotal + oneIndex);
			numParamEqCons++;
		}
		std::sort(sparsity2.begin(), sparsity2.end());
        mosek::fusion::Parameter::t EM = M->parameter(monty::new_array_ptr<int>({numParamEqCons, varsTotal}), monty::new_array_ptr<long>(sparsity2));
        M->constraint(mosek::fusion::Expr::mul(EM, xM), mosek::fusion::Domain::equalsTo(0.0));

		// The first element of the vector should be one
		M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));

		// Linear equality constraints
		M->constraint(mosek::fusion::Expr::mul(AM, xM), mosek::fusion::Domain::equalsTo(0.0));

		// Linear positivity constraints
		M->constraint(mosek::fusion::Expr::mul(BM, xM), mosek::fusion::Domain::greaterThan(0));

		// Linear objective function
		M->objective(mosek::fusion::ObjectiveSense::Minimize, mosek::fusion::Expr::dot(cM, xM));

		// Moment matrix constraints
		for (int i=0; i<shouldBePSD.size(); i++) {
			M->constraint(mosek::fusion::Expr::mulElm(shouldBePSDCoeffs[i], xM->pick(shouldBePSD[i])), mosek::fusion::Domain::inSVecPSDCone());
		}

		// Linear PSD constraints
		if (conPSDLinear.size() > 0) {
			mosek::fusion::Expression::t matSum = mosek::fusion::Expr::mulElm(psdCoeffsM[0], xM->pick(psdLocsM[0]));
			for (int i=1; i<numMatsToSum; i++) {
				matSum = mosek::fusion::Expr::add(matSum, mosek::fusion::Expr::mulElm(psdCoeffsM[i], xM->pick(psdLocsM[i])));
			}
			M->constraint(matSum, mosek::fusion::Domain::inSVecPSDCone());
		}

		// Keep splitting
		int iter = 1;
		double totalArea = 0;
		double totalAreaInfeasible = 0;
		double areaCovered = 1;
		auto begin = std::chrono::steady_clock::now();
		auto end = std::chrono::steady_clock::now();
		double secondsPerIter = 0;
		int numIllPosed = 0;
		double upperBound = largeTol;
		double lowerBound = -largeTol;
		std::vector<double> lowerBounds = {lowerBound};
		std::vector<double> bestFeasiblePoint;
		while (toProcess.size() > 0) {

			// Debug output
			if (verbosity >= 2) {
				std::cout << "    checking region: " << toProcess[0] << std::endl;
			}

			// The area taken up by this section
			areaCovered = 1;
			for (int k=0; k<toProcess[0].size(); k++) {
				if (varIsBinary[k]) {
					if (toProcess[0][k].first == toProcess[0][k].second) {
						areaCovered *= 1.0;
					} else {
						areaCovered *= 2.0;
					}
				} else {
					areaCovered *= toProcess[0][k].second - toProcess[0][k].first;
				}
			}

			// Update the linear constraints on the quadratics
            if (!isLin) {
                std::vector<std::vector<double>> newD(numParamPosCons, std::vector<double>(varsTotal, 0));
                std::vector<std::vector<double>> newE(numParamEqCons, std::vector<double>(varsTotal, 0));
                for (int i=0; i<toProcess[0].size(); i++) {

                    // If it's a binary variable
                    if (varIsBinary[i]) {

                        // If it's been constrained to a value
                        if (toProcess[0][i].first == toProcess[0][i].second) {

                            // Add this as a linear equality con
                            newE[i][oneIndex] = -toProcess[0][i].first;
                            newE[i][firstMonomInds[i]] = 1;

                        }

                    } else {

                        // Given two points, find ax+by+c=0
                        std::vector<double> point1 = {toProcess[0][i].first, toProcess[0][i].first*toProcess[0][i].first};
                        std::vector<double> point2 = {toProcess[0][i].second, toProcess[0][i].second*toProcess[0][i].second};
                        std::vector<double> coeffs = getLineFromPoints(point1, point2);

                        // Add this as a linear positivity con
                        newD[i][oneIndex] = coeffs[0];
                        newD[i][firstMonomInds[i]] = coeffs[1];
                        newD[i][quadraticMonomInds[i][i]] = coeffs[2];

                    }

                }
                DM->setValue(monty::new_array_ptr<double>(newD));
                EM->setValue(monty::new_array_ptr<double>(newE));
            }

			// Solve the problem
			M->solve();
			auto statProb = M->getProblemStatus();
			auto statSol = M->getPrimalSolutionStatus();

			// If the problem was feasible, extract the result and figure out where to split
			if (statProb == mosek::fusion::ProblemStatus::PrimalAndDualFeasible && statSol == mosek::fusion::SolutionStatus::Optimal) {

				// Get the solution values
				auto sol = *(xM->level());
				polyType objPrimal = M->primalObjValue();
				polyType objDual = M->dualObjValue();

				// Output the relevant moments
				std::vector<polyType> solVec(xM->getSize());
				for (int i=0; i<solVec.size(); i++) {
					solVec[i] = sol[i];
				}

				// Get the nearest feasible point
				std::vector<double> x(maxVariables, 0);
				for (int i=0; i<maxVariables; i++) {
					x[i] = solVec[firstMonomInds[i]];
				}
				std::vector<double> nearestFeasiblePoint(x.size());
                if (isLin) {
                    nearestFeasiblePoint = x;
                } else {
                    nearestFeasiblePoint = getNearestFeasiblePoint(x, toProcess[0]);
                }

				// Output if verbose
				if (verbosity >= 2) {
					std::cout << "    obj primal: " << objPrimal << std::endl;
					std::cout << "    x: " << x << std::endl;
					std::cout << "    nearest feasible point: " << nearestFeasiblePoint << std::endl;
				}

                // If the point is feasible
                if (nearestFeasiblePoint.size() > 0) {

                    // Use this to give an upper bound
                    double objEval = obj.eval(nearestFeasiblePoint);
                    if (objEval < upperBound) {
                        upperBound = objEval;
                        bestFeasiblePoint = nearestFeasiblePoint;
                    }

                    // Output if verbose
                    if (verbosity >= 2) {
                        std::cout << "    eval at this point: " << objEval << std::endl;
                    }

                }

				// Check the resulting vector for a good place to split 
				std::vector<double> errors(maxVariables);
                //if (nearestFeasiblePoint.size() > 0) {
                    //for (int i=0; i<maxVariables; i++) {
                        //errors[i] = std::pow(nearestFeasiblePoint[i] - x[i], 2);
                    //}
                //}
                for (int i=0; i<maxVariables; i++) {
                    errors[i] = std::pow(solVec[quadraticMonomInds[i][i]] - solVec[firstMonomInds[i]]*solVec[firstMonomInds[i]], 2);
                }

				// Find the biggest error
				double biggestError = -10000;
				int bestInd = -1;
				for (int i=0; i<maxVariables; i++) {
					if (errors[i] > biggestError && toProcess[0][i].first != toProcess[0][i].second) {
						biggestError = errors[i];
						bestInd = i;
					}
				}

				// Output if verbose
				if (verbosity >= 2) {
					std::cout << "    errors: " << errors << std::endl;
				}

				// If the area could still contain the min
				if (objPrimal < upperBound) {

					// Split it 
					double minPoint = toProcess[0][bestInd].first;
					double maxPoint = toProcess[0][bestInd].second;
					double midPoint = (minPoint + maxPoint) / 2.0;
					double mostFeasiblePoint = solVec[firstMonomInds[bestInd]];
					double splitPoint = mostFeasiblePoint;
					if (verbosity >= 2) {
						std::cout << "    splitting var " << bestInd << " at " << splitPoint << std::endl;
					}
					auto copyLeft = toProcess[0];
					auto copyRight = toProcess[0];

					// If it's binary, just split it to the two points
					if (varIsBinary[bestInd]) {
						copyLeft[bestInd].second = toProcess[0][bestInd].first;
						copyRight[bestInd].first = toProcess[0][bestInd].second;

					// Otherwise at the point
					} else {
						copyLeft[bestInd].second = splitPoint;
						copyRight[bestInd].first = splitPoint;
					}

					// Add the new paths to the queue, lowest lowerBound first
					int queueLoc = toProcess.size();
					for (int i=1; i<toProcess.size(); i++) {
						if (objPrimal < lowerBounds[i]) {
							queueLoc = i;
							break;
						}
					}
					toProcess.insert(toProcess.begin()+queueLoc, copyLeft);
					toProcess.insert(toProcess.begin()+queueLoc, copyRight);
					lowerBounds.insert(lowerBounds.begin()+queueLoc, objPrimal);
					lowerBounds.insert(lowerBounds.begin()+queueLoc, objPrimal);

				// Otherwise we can say that this region is done
				} else {
					totalArea += areaCovered;

				}

			// If it's infeasible
			} else {
				totalArea += areaCovered;
				totalAreaInfeasible += areaCovered;
			}

			// Remove the one we just processed
			toProcess.erase(toProcess.begin());
			lowerBounds.erase(lowerBounds.begin());

			// The lower bound is the lowest lower bound
			lowerBound = lowerBounds[0];

            // If it's linear
            if (isLin) {
                lowerBound = upperBound;
            }

			// Time estimation
			end = std::chrono::steady_clock::now();
			secondsPerIter = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / (iter * 1.0e6);
			double areaPerIter = totalArea / iter;
			double itersRemaining = (maxArea - totalArea) / areaPerIter;
			double secondsRemaining = itersRemaining * secondsPerIter;

			// Per-iteration output
			std::cout << std::defaultfloat;
			if (verbosity >= 1) {
				std::cout << iter << "i"
							<< " * " << representTime(secondsPerIter) << "/i" 
							<< " = " << representTime(iter*secondsPerIter)
							<< "   cov:" << 100.0 * totalArea / maxArea << "%" 
							<< "   est:" << representTime(secondsRemaining) << "  " 
							<< lowerBound << " <= " << upperBound;
				if (verbosity >= 2) {
					std::cout << std::endl;
				} else {
					std::cout << "\r" << std::flush;
				}
			}
            
			// If we've converged
			if (std::abs(lowerBound - upperBound) < 1e-5) {
				if (verbosity >= 1) {
					std::cout << std::endl;
				}
				return {upperBound, bestFeasiblePoint};
			}

			// Keep track of the iteration number
			iter++;
			if (maxIters >= 0 && iter > maxIters) {
				break;
			}	

		}

		// If we got here it's infeasible
		if (verbosity >= 1) {
			std::cout << std::endl;
		}
		if (std::abs(totalAreaInfeasible - maxArea) < 1e-5) {
			return {0.0, {}};
		} else {
			return {upperBound, bestFeasiblePoint};
		}

	}

	// Attempt to find a series of constraints that show this is infeasible
	void proveInfeasible(int maxIters=-1, std::string level="1f", double bound=1, std::string logFileName="", int verbosity=1, int numVarsToSplit=0) {

		// Get the monomial list and sort it
		std::vector<std::string> monoms = getMonomials();
		std::sort(monoms.begin(), monoms.end(), [](const std::string& first, const std::string& second){return first.size() < second.size();});

		// List of semdefinite matrices
		std::vector<std::vector<Polynomial<polyType>>> monomProducts;

		// Get the order of the set of equations
		int degree = getDegree();

		// Make sure we have all first order moments
		addMonomsOfOrder(monoms, 1);

		// First monom should always be 1
		auto loc = std::find(monoms.begin(), monoms.end(), "");
		if (loc != monoms.end()) {
			monoms.erase(loc);
		}
		monoms.insert(monoms.begin(), "");

		// Add all square monoms
		if (degree >= 2) {
			for (int i=0; i<maxVariables; i++) {
				std::string newInd = std::to_string(i);
				newInd.insert(0, digitsPerInd-newInd.size(), ' ');
				if (std::find(monoms.begin(), monoms.end(), newInd+newInd) == monoms.end()) {
					monoms.push_back(newInd+newInd);
				}
			}
		}

		// Add all quartic monoms
		if (degree >= 4) {
			for (int i=0; i<maxVariables; i++) {
				std::string newInd = std::to_string(i);
				newInd.insert(0, digitsPerInd-newInd.size(), ' ');
				if (std::find(monoms.begin(), monoms.end(), newInd+newInd+newInd+newInd) == monoms.end()) {
					monoms.push_back(newInd+newInd+newInd+newInd);
				}
			}
		}
		
		// Create the mapping from monomials to indices (to linearize)
		int numOGMonoms = monoms.size();
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
		std::vector<std::vector<Polynomial<polyType>>> conPSDLinear;
		for (int i=0; i<conPSD.size(); i++) {
			conPSDLinear.push_back({});
			for (int j=0; j<conPSD[i].size(); j++) {
				conPSDLinear[i].push_back(conPSD[i][j].replaceWithVariable(mapping));
			}
		}

		// Process the level string e.g. "1+2f,3p"
		std::vector<int> levelsToInclude = {};
		std::string currentThing = "";
		for (int i=0; i<level.size(); i++) {

			// For multiple levels e.g. 1+2
			if (level[i] == '+') {
				levelsToInclude.push_back(std::stoi(currentThing));
				currentThing = "";

			// If told to do a certain level for each individual constraint
			} else if (level[i] == 's') {

				// Add whatever numbers are left
				levelsToInclude.push_back(std::stoi(currentThing));

                // Get the maximum level
                int l = 1;
                for (int j=0; j<levelsToInclude.size(); j++) {
                    l = std::max(l, levelsToInclude[j]);
                }

                // For each linear constraint
                for (int j=0; j<conZero.size(); j++) {
                    std::vector<int> possibleVars = conZero[j].getVariables();
                    std::vector<Polynomial<polyType>> asPolys = {Polynomial<polyType>(maxVariables, 1)};
                    if (l >= 1) {
                        for (int k=0; k<possibleVars.size(); k++) {
                            asPolys.push_back(Polynomial<polyType>(maxVariables, 1, {possibleVars[k]}));
                        }
                    }
                    if (l >= 2) {
                        for (int k=0; k<possibleVars.size(); k++) {
                            for (int m=k; m<possibleVars.size(); m++) {
                                asPolys.push_back(Polynomial<polyType>(maxVariables, 1, {possibleVars[k], possibleVars[m]}));
                            }
                        }
                    }
                    if (l >= 3) {
                        for (int k=0; k<possibleVars.size(); k++) {
                            for (int m=k; m<possibleVars.size(); m++) {
                                for (int n=m; n<possibleVars.size(); n++) {
                                    asPolys.push_back(Polynomial<polyType>(maxVariables, 1, {possibleVars[k], possibleVars[m], possibleVars[n]}));
                                }
                            }
                        }
                    }
                    monomProducts.push_back(asPolys);
                }

                // Reset things
                currentThing = "";
                levelsToInclude = {};

			// p for a partial level
			} else if (level[i] == 'p') {

				// Add whatever numbers are left
				levelsToInclude.push_back(std::stoi(currentThing));

				// Get the list monomials that can appear on the top row
				std::vector<std::string> possibleMonoms;
				for (int j=0; j<levelsToInclude.size(); j++) {
					addMonomsOfOrder(possibleMonoms, levelsToInclude[j]);
				}

				// Get all combinations of 3x3 moment matrices
				for (int j=0; j<possibleMonoms.size(); j++) {
					for (int k=j+1; k<possibleMonoms.size(); k++) {
						monomProducts.push_back({Polynomial<polyType>(maxVariables, 1), Polynomial<polyType>(maxVariables, 1, possibleMonoms[j]), Polynomial<polyType>(maxVariables, 1, possibleMonoms[k])});
					}
				}

				// Reset stuff
				currentThing = "";
				levelsToInclude = {};

			// f for a full level
			} else if (level[i] == 'f') {

				// Add whatever numbers are left
				levelsToInclude.push_back(std::stoi(currentThing));

				// Get the list monomials that can appear on the top row
				std::vector<std::string> possibleMonoms;
				for (int j=0; j<levelsToInclude.size(); j++) {
					addMonomsOfOrder(possibleMonoms, levelsToInclude[j]);
				}

				// Add these all to the top row of a moment matrix
				std::vector<Polynomial<double>> toAdd;
				toAdd.push_back(Polynomial<polyType>(maxVariables, 1));
				for (int j=0; j<possibleMonoms.size(); j++) {
					toAdd.push_back(Polynomial<polyType>(maxVariables, 1, possibleMonoms[j]));
				}
				monomProducts.push_back(toAdd);

				// Reset stuff
				currentThing = "";
				levelsToInclude = {};

			// Otherwise add this (probably digit) to the string to process
			} else if (level[i] != ',' && level[i] != ' ') {
				currentThing += level[i];
			}

		}

		// Start with the most general area
		double maxArea = 1;
		std::vector<std::vector<std::pair<double,double>>> toProcess;
		std::vector<std::pair<double,double>> varMinMax(obj.maxVariables);
		for (int i=0; i<obj.maxVariables; i++) {
			varMinMax[i].first = -bound;
			varMinMax[i].second = bound;
			maxArea *= 2*bound;
		}
		toProcess.push_back(varMinMax);

		// If told to start by splitting at various vars
		if (numVarsToSplit > 0) {

			// Create the lists of vars to split, starting from the end
			toProcess = {};
			std::vector<int> varsToSplit = {};
			for (int i=maxVariables-1; i>maxVariables-numVarsToSplit-1; i--) {
				varsToSplit.push_back(i);
			}

			// It's exponential in the number of vars, do every combination
			for (int i=0; i<std::pow(2, varsToSplit.size()); i++) {
				auto regionCopy = varMinMax;
				for (int j=0; j<varsToSplit.size(); j++) {
					if (i >> j & 1) {
						regionCopy[varsToSplit[j]].second = 0;
					} else {
						regionCopy[varsToSplit[j]].first = 0;
					}
				}
				toProcess.push_back(regionCopy);
			}

		}

		// Verbose output
		if (verbosity >= 2) {
			std::cout << "extra PSD matrices: " << std::endl;
			for (int i=0; i<monomProducts.size(); i++) {
				std::cout << monomProducts[i] << std::endl;
			}
		}

		// Create the PSD matrices from this list
		std::vector<std::shared_ptr<monty::ndarray<int,1>>> shouldBePSD;
		std::vector<std::shared_ptr<monty::ndarray<double,1>>> shouldBePSDCoeffs;
		std::vector<std::shared_ptr<monty::ndarray<double,1>>> identityAsSVec;
		std::vector<double> idenCoeffs;
		std::vector<double> monCoeffs;
		std::vector<int> monLocs;
		std::unordered_map<std::string,int> monomsInverted;
		for (int j=0; j<monoms.size(); j++) {
			monomsInverted[monoms[j]] = j;
		}
		for (int j=0; j<monomProducts.size(); j++) {

			// Get the list of all monomial locations for the PSD matrix
			monLocs = {};
			monCoeffs = {};
			idenCoeffs = {};
			for (int i=0; i<monomProducts[j].size(); i++) {
				for (int k=i; k<monomProducts[j].size(); k++) {

					// Calculate the product
					std::string monString = (monomProducts[j][i]*monomProducts[j][k]).getMonomials()[0];

					// Find this in the monomial list
					auto loc = monomsInverted.find(monString);
					if (loc != monomsInverted.end()) {
						monLocs.push_back(monomsInverted[monString]);
					} else {
						monomsInverted[monString] = monoms.size();
						monLocs.push_back(monoms.size());
						monoms.push_back(monString);
					}

					// The coeff for mosek's svec
					if (i != k) {
						monCoeffs.push_back(std::sqrt(2.0));
						idenCoeffs.push_back(0.0);
					} else {
						monCoeffs.push_back(1.0);
						idenCoeffs.push_back(1.0);
					}

				}
			}

			// This (when reformatted) should be positive-semidefinite
			shouldBePSD.push_back(monty::new_array_ptr<int>(monLocs));
			shouldBePSDCoeffs.push_back(monty::new_array_ptr<double>(monCoeffs));
			identityAsSVec.push_back(monty::new_array_ptr<double>(idenCoeffs));

			// Clear some memory
			monomProducts.erase(monomProducts.begin());
			j--;

		}

		// Get the inds of the first order monomials and their squares
		std::vector<int> firstMonomInds(varMinMax.size(), -1);
		std::vector<std::vector<int>> quadraticMonomInds(varMinMax.size(), std::vector<int>(varMinMax.size(), 1));
		for (int i=0; i<monoms.size(); i++) {
			if (monoms[i].size() == digitsPerInd) {
				firstMonomInds[std::stoi(monoms[i])] = i;
			} else if (monoms[i].size() == 2*digitsPerInd) {
				int ind1 = std::stoi(monoms[i].substr(0,digitsPerInd));
				int ind2 = std::stoi(monoms[i].substr(digitsPerInd,digitsPerInd));
				quadraticMonomInds[ind1][ind2] = i;
				quadraticMonomInds[ind2][ind1] = i;
			}
		}

		// Set some vars
		int oneIndex = 0;
		int varsTotal = monoms.size();

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

		// Determine the largest amount of terms per element of the PSD matrix
		int numMatsToSum = 1;
		for (int i=0; i<conPSDLinear.size(); i++) {
			for (int j=i; j<conPSDLinear[i].size(); j++) {
				numMatsToSum = std::max(numMatsToSum, int(conPSDLinear[i][j].coeffs.size()));
			}
		}

		// Convert the linear PSD constraints to MOSEK form
		int upperTriangSize = (conPSDLinear.size()*(conPSDLinear.size()+1)) / 2;
		std::vector<std::vector<int>> psdLocs(numMatsToSum, std::vector<int>(upperTriangSize, oneIndex));
		std::vector<std::vector<double>> psdCoeffs(numMatsToSum, std::vector<double>(upperTriangSize, 0.0));
		int whichEl = 0;
		for (int i=0; i<conPSDLinear.size(); i++) {
			for (int j=i; j<conPSDLinear[i].size(); j++) {
				int whichMat = 0;
				if (conPSDLinear[i][j].coeffs.size() > 0) {
					for (auto const &pair: conPSDLinear[i][j].coeffs) {
						if (pair.first == "") {
							psdLocs[whichMat][whichEl] = oneIndex;
						} else {
							psdLocs[whichMat][whichEl] = std::stoi(pair.first);
						}
						if (i == j) {
							psdCoeffs[whichMat][whichEl] = pair.second;
						} else {
							psdCoeffs[whichMat][whichEl] = std::sqrt(2.0)*pair.second;
						}
						whichMat++;
					}
				} else {
					psdLocs[whichMat][whichEl] = oneIndex;
					psdCoeffs[whichMat][whichEl] = 0;
				}
				whichEl++;
			}
		}
		std::vector<std::shared_ptr<monty::ndarray<int,1>>> psdLocsM(numMatsToSum);
		std::vector<std::shared_ptr<monty::ndarray<double,1>>> psdCoeffsM(numMatsToSum);
		for (int i=0; i<numMatsToSum; i++) {
			psdLocsM[i] = monty::new_array_ptr<int>(psdLocs[i]);
			psdCoeffsM[i] = monty::new_array_ptr<double>(psdCoeffs[i]);
		}

		// Create a model
		mosek::fusion::Model::t M = new mosek::fusion::Model(); auto _M = monty::finally([&]() {M->dispose();});

		// Create the variable
		mosek::fusion::Variable::t xM = M->variable(varsTotal);

		// Use an extra variable to minimize violation of SD
		mosek::fusion::Variable::t lambda = M->variable();

		// Parameterized linear positivity constraints
		int numParamCons = maxVariables;
		std::vector<long> sparsity;
		for (int i=0; i<toProcess[0].size(); i++) {
			sparsity.push_back(i*varsTotal + firstMonomInds[i]);
			sparsity.push_back(i*varsTotal + quadraticMonomInds[i][i]);
			sparsity.push_back(i*varsTotal + oneIndex);
		}
		std::sort(sparsity.begin(), sparsity.end());
		mosek::fusion::Parameter::t DM = M->parameter(monty::new_array_ptr<int>({numParamCons, varsTotal}), monty::new_array_ptr<long>(sparsity));
		M->constraint(mosek::fusion::Expr::mul(DM, xM), mosek::fusion::Domain::greaterThan(0));

		// The first element of the vector should be one
		M->constraint(xM->index(oneIndex), mosek::fusion::Domain::equalsTo(1.0));

		// Linear equality constraints
		M->constraint(mosek::fusion::Expr::mul(AM, xM), mosek::fusion::Domain::equalsTo(0.0));

		// Linear positivity constraints
		M->constraint(mosek::fusion::Expr::mul(BM, xM), mosek::fusion::Domain::greaterThan(0));

		// Try to violate the SDP constraints the least
		M->objective(mosek::fusion::ObjectiveSense::Minimize, lambda);

		// Moment matrix constraints
		for (int i=0; i<shouldBePSD.size(); i++) {
			M->constraint(mosek::fusion::Expr::add(mosek::fusion::Expr::mulElm(shouldBePSDCoeffs[i], xM->pick(shouldBePSD[i])), mosek::fusion::Expr::mul(lambda, identityAsSVec[i])), mosek::fusion::Domain::inSVecPSDCone());
		}

		// Linear PSD constraints
		if (conPSDLinear.size() > 0) {
			mosek::fusion::Expression::t matSum = mosek::fusion::Expr::mulElm(psdCoeffsM[0], xM->pick(psdLocsM[0]));
			for (int i=1; i<numMatsToSum; i++) {
				matSum = mosek::fusion::Expr::add(matSum, mosek::fusion::Expr::mulElm(psdCoeffsM[i], xM->pick(psdLocsM[i])));
			}
			M->constraint(matSum, mosek::fusion::Domain::inSVecPSDCone());
		}

		// Open a file if told to record the points
		std::ofstream logFile;
		if (logFileName.size() > 0) {
			logFile.open(logFileName);
			for (int i=0; i<toProcess[0].size(); i++) {
				logFile << i << ", ";
			}
			logFile << "feasibility";
			logFile << std::endl;
		}

		// Keep splitting until all fail
		int iter = 1;
		double totalArea = 0;
		double areaCovered = 1;
		auto begin = std::chrono::steady_clock::now();
		auto end = std::chrono::steady_clock::now();
		double secondsPerIter = 0;
		int numIllPosed = 0;
		while (toProcess.size() > 0) {

			// Debug output
			if (verbosity >= 2) {
				std::cout << "    checking region: " << toProcess[0] << std::endl;
			}

			// The area taken up by this section
			areaCovered = 1;
			for (int k=0; k<toProcess[0].size(); k++) {
				areaCovered *= toProcess[0][k].second - toProcess[0][k].first;
			}

			// The new parameterized constraint vector and objective
			std::vector<std::vector<double>> newD(numParamCons, std::vector<double>(varsTotal, 0));

			// Update the linear constraints on the quadratics
			int nextInd = 0;
			for (int i=0; i<toProcess[0].size(); i++) {

				// Given two points, find ax+by+c=0
				std::vector<double> point1 = {toProcess[0][i].first, toProcess[0][i].first*toProcess[0][i].first};
				std::vector<double> point2 = {toProcess[0][i].second, toProcess[0][i].second*toProcess[0][i].second};
				std::vector<double> coeffs = getLineFromPoints(point1, point2);

				// Add this as a linear pos con
				newD[nextInd][oneIndex] = coeffs[0];
				newD[nextInd][firstMonomInds[i]] = coeffs[1];
				newD[nextInd][quadraticMonomInds[i][i]] = coeffs[2];
				nextInd++;

			}

			// Solve the problem
			DM->setValue(monty::new_array_ptr<double>(newD));
			M->solve();
			auto statProb = M->getProblemStatus();
			auto statSol = M->getPrimalSolutionStatus();

			// If infeasible, good
			if (statProb == mosek::fusion::ProblemStatus::PrimalInfeasible || statSol == mosek::fusion::SolutionStatus::Undefined || statSol == mosek::fusion::SolutionStatus::Unknown || M->dualObjValue() > 1e-7) {

				// Keep track of how many were ill-posed
				if (statSol == mosek::fusion::SolutionStatus::Undefined || statSol == mosek::fusion::SolutionStatus::Unknown) {
					numIllPosed++;
				}

				// Write to a file if told to
				if (logFileName.size() > 0 && statSol != mosek::fusion::SolutionStatus::Undefined && statSol != mosek::fusion::SolutionStatus::Unknown && statProb != mosek::fusion::ProblemStatus::PrimalInfeasible) {
					auto sol = *(xM->level());
					std::vector<polyType> solVec(xM->getSize());
					for (int i=0; i<solVec.size(); i++) {
						solVec[i] = sol[i];
					}
					for (int i=0; i<toProcess[0].size(); i++) {
						logFile << solVec[firstMonomInds[i]] << ", ";
					}
					logFile << "no";
					logFile << std::endl;
				}

				// Update the total area count
				totalArea += areaCovered;

			// Otherwise, extract the result and figure out where to split
			} else {

				// Get the solution values
				auto sol = *(xM->level());
				polyType objPrimal = M->primalObjValue();
				polyType objDual = M->dualObjValue();

				// Output the relevant moments
				std::vector<polyType> solVec(xM->getSize());
				for (int i=0; i<solVec.size(); i++) {
					solVec[i] = sol[i];
				}

				// Check the resulting vector for a good place to split 
				std::vector<double> errors(maxVariables);
				for (int i=0; i<monoms.size(); i++) {
					if (monoms[i].size() == 2*digitsPerInd) {
						int ind1 = std::stoi(monoms[i].substr(0, digitsPerInd));
						int ind2 = std::stoi(monoms[i].substr(digitsPerInd, digitsPerInd));
						errors[ind1] += std::pow(solVec[i] - solVec[firstMonomInds[ind1]]*solVec[firstMonomInds[ind2]], 2);
						errors[ind2] += std::pow(solVec[i] - solVec[firstMonomInds[ind1]]*solVec[firstMonomInds[ind2]], 2);
					}
				}

				// Find the biggest error
				double biggestError = -10000;
				int bestInd = -1;
				for (int i=0; i<maxVariables; i++) {
					if (errors[i] > biggestError) {
						biggestError = errors[i];
						bestInd = i;
					}
				}

				// If we've converged
				if (biggestError < 1e-8) {

					// If logging, write to file and continue
					if (logFileName.size() > 0) {
						for (int i=0; i<toProcess[0].size(); i++) {
							logFile << solVec[firstMonomInds[i]] << ", ";
						}
						logFile << "yes";
						logFile << std::endl;

					// Otheriwse write to std::out and stop
					} else {
						std::cout << std::endl;
						std::vector<double> vals(maxVariables);
						for (int i=0; i<maxVariables; i++) {
							vals[i] = solVec[firstMonomInds[i]];
						}
						std::cout << "converged in " << iter << " iters to: " << vals << std::endl;
						break;
					}

				// If there's still space to split
				} else {

					// Split it 
					double minPoint = toProcess[0][bestInd].first;
					double maxPoint = toProcess[0][bestInd].second;
					double midPoint = (minPoint + maxPoint) / 2.0;
					double mostFeasiblePoint = solVec[firstMonomInds[bestInd]];
					double distanceBetween = mostFeasiblePoint - midPoint;
					double splitPoint = mostFeasiblePoint;
					if (verbosity >= 2) {
						std::cout << "    splitting var " << bestInd << " at " << splitPoint << std::endl;
					}
					auto copyLeft = toProcess[0];
					auto copyRight = toProcess[0];
					copyLeft[bestInd].second = splitPoint;
					copyRight[bestInd].first = splitPoint;

					// Add the new paths to the queue
					toProcess.insert(toProcess.begin()+1, copyLeft);
					toProcess.insert(toProcess.begin()+1, copyRight);

				}

			}

			// Time estimation
			end = std::chrono::steady_clock::now();
			secondsPerIter = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / (iter * 1.0e6);
			double areaPerIter = totalArea / iter;
			double itersRemaining = (maxArea - totalArea) / areaPerIter;
			double secondsRemaining = itersRemaining * secondsPerIter;

			// Per-iteration output
			std::cout << std::defaultfloat;
			if (verbosity >= 2) {
				std::cout << iter << "i  " << 100.0 * totalArea / maxArea << "%  " << 100.0 * areaPerIter / maxArea << "%/i  " << representTime(secondsPerIter) << "/i  " << numIllPosed << "  " << representTime(secondsRemaining) << "  " << 100.0 * areaCovered / maxArea << "%\n" << std::flush;
			} else if (verbosity >= 1) {
				std::cout << iter << "i  " << 100.0 * totalArea / maxArea << "%  " << 100.0 * areaPerIter / maxArea << "%/i  " << representTime(secondsPerIter) << "/i  " << numIllPosed << "  " << representTime(secondsRemaining) << "  " << 100.0 * areaCovered / maxArea << "%                  \r" << std::flush;
			}

			// Remove the one we just processed
			toProcess.erase(toProcess.begin());

			// Keep track of the iteration number
			iter++;
			if (maxIters >= 0 && iter > maxIters) {
				break;
			}	

		}
		if (verbosity >= 1) { 
			std::cout << std::endl;
			std::cout << representTime(iter * secondsPerIter) << std::endl;
		}

	}

};

// Generic overload for outputting vector of vector of Polynomials
template <typename type> 
std::ostream &operator<<(std::ostream &output, const std::vector<std::vector<Polynomial<type>>> &arr) {

	// Use fixed precision
	int precision = 7;
	output << std::showpos << std::fixed << std::setprecision(precision);

	// For each column, figure out the output widths
	int columns = 0;
	if (arr.size() > 0) {
		columns = arr[0].size();
	}
	std::vector<int> columnWidths(columns, 0);
	for (int i=0; i<arr.size(); i++) {
		for (int j=0; j<arr[i].size(); j++) {
			columnWidths[j] = std::max(columnWidths[j], arr[i][j].getOutputWidth()+1);
		}
	}

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
			int paddingNeeded = columnWidths[x]-arr[y][x].getOutputWidth();
			output << std::string(paddingNeeded, ' ');
			output << arr[y][x];
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
	
	// Reset formatting
	output << std::noshowpos << std::defaultfloat;
	
	return output;

}

#endif
	
