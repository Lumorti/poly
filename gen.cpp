#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <complex>
#include <iomanip>
#include <random>
#include <fstream>
#include <chrono>
#include <math.h>

// Print an unordered map
template<typename K, typename V>
void pretty(std::unordered_map<K, V> const &m) {
    for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}

// Print an equation
void pretty(std::unordered_map<int,std::string> invertedMap, std::vector<std::pair<int,int>> a) {
	std::cout << "[";
	for (int i=0; i<a.size(); i++) {
		std::cout << invertedMap[a[i].first] << ": " << a[i].second << ", ";
	}
	std::cout << "]" << std::endl;
}

// Multiply two systems e.g. x1^2 * y1^1 = 1^2,1^1
std::string multiply(std::string a, std::string b) {
	return a + "," + b;
}

// Standard cpp entry point
int main(int argc, char ** argv) {

	// Get the problem from the args
	int d = 2;
	int n = 2;
	int k = 1;
	int numVarsNonConj = n*d*d;
	int numVars = 2*n*d*d;
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}
	if (argc > 3) {
		k = std::stoi(argv[3]);
	}

	// Init the monomial ordering
	int conjDelta = numVarsNonConj;
	int nextInd = 1;
	std::unordered_map<std::string, int> map = {{"", 0}};
	std::vector<std::string> varNames;
	for (int i1=0; i1<numVars; i1++) {
		varNames.push_back(std::to_string(i1));
	}

	// Linear terms
	std::vector<std::vector<std::pair<int, int>>> eqns;
	for (int i1=0; i1<numVars; i1++) {
		map[varNames[i1]] = nextInd;
		nextInd++;
	}
	
	// Quadratic terms
	for (int i1=0; i1<numVars; i1++) {
		for (int i2=i1+1; i2<numVars; i2++) {
			map[multiply(varNames[i1], varNames[i2])] = nextInd;
			nextInd++;
		}
	}

	// Invert the monomial map
	std::unordered_map<int, std::string> invertedMap;
    for (auto const &pair: map) {
		invertedMap[pair.second] = pair.first;
    }

	// Generate normalisation equations
	for (int i=0; i<n; i++) {
		for (int j=0; j<d; j++) {
			std::vector<std::pair<int,int>> eqn = {{0, -1}};
			for (int k=0; k<d; k++) {
				int ind = i*d*d + j*d + k;
				eqn.push_back(std::pair<int,int>(map[multiply(varNames[ind], varNames[ind+conjDelta])], 1));
			}
			eqns.push_back(eqn);
		}
	}

	// Generate orthogonality equations
	for (int i=0; i<n; i++) {
		for (int j=0; j<d; j++) {
			for (int k=j+1; k<d; k++) {
				std::vector<std::pair<int,int>> eqn = {};
				std::vector<std::pair<int,int>> eqnConj = {};
				for (int l=0; l<d; l++) {
					int ind1 = i*d*d + j*d + l;
					int ind2 = i*d*d + k*d + l;
					eqn.push_back(std::pair<int,int>(map[multiply(varNames[ind1], varNames[ind2+conjDelta])], 1));
					eqnConj.push_back(std::pair<int,int>(map[multiply(varNames[ind2], varNames[ind1+conjDelta])], 1));
				}
				eqns.push_back(eqn);
				eqns.push_back(eqnConj);
			}
		}
	}

	// Generate MUB equations TODO
	
	// Check everything
	pretty(map);
	for (int i=0; i<eqns.size(); i++) {
		pretty(invertedMap, eqns[i]);
	}

	// Get the polynomial to multiply by

	// Multiply each equation by each element of this polynomial 

	// Write to file

}
	
