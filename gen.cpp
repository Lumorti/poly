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

// TODO
// https://www.math.ucdavis.edu/~deloera/RECENT_WORK/jsc09_issac08.pdf

// Print an unordered map
template<typename K, typename V>
void pretty(std::unordered_map<K, V> const &m) {
    for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}

// Print an equation
void pretty(std::unordered_map<int,std::string> invertedMap, std::vector<std::pair<int,char>> a) {
	std::cout << "[";
	for (int i=0; i<a.size(); i++) {
		std::cout << invertedMap[a[i].first] << ": " << a[i].second << ", ";
	}
	std::cout << "]" << std::endl;
}

// Multiply two terms together, sorting the result TODO
// e.g. 1,3,7, * 2,5, = 1,2,3,5,7,
// this should be O(n+m)
std::string multiplyFull(std::string a, std::string b) {

	// Start with a blank string, go through and add the other items
	std::string newString = "";

	// Loop through the first string
	int firstStart = 0;
	int secondStart = 0;
	int j = 0;
	for (int i=0; i<a.size(); i++) {

		// Stop at each comma, we have a new value
		if (a[i] == ',') {

			// Get the value
			std::string string1 = a.substr(firstStart, i-firstStart);
			int val1 = std::stoi(string1);
			firstStart = i+1;

			// Go through the other string
			while (j<b.size()) {

				// Stop at each comma, we have a new value
				if (b[j] == ',') {

					// Get the value
					std::string string2 = b.substr(secondStart, j-secondStart);
					int val2 = std::stoi(string2);
					secondStart = j+1;

					// If this is less, add it
					if (val2 < val1) {
						newString += string2 + ",";
					} else {
						j -= j-secondStart+1;
						secondStart = j-1;
						break;
					}

				}

				j++;

			}

			// Then finally add the first string
			newString += string1 + ",";

		}

	}

	return newString;
}

// Standard cpp entry point
int main(int argc, char ** argv) {

	// Get the problem from the args
	int d = 2;
	int n = 2;
	int k = 0;
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}
	if (argc > 3) {
		k = std::stoi(argv[3]);
	}
	int numVarsNonConj = n*d*d;
	int numVars = 2*n*d*d;

	// Init the monomial ordering
	int conjDelta = numVarsNonConj;
	int nextInd = 1;
	std::unordered_map<std::string, int> map = {{"", 0}};
	std::unordered_map<int, std::string> invertedMap;
	std::vector<std::vector<std::pair<int,char>>> eqns;
	std::vector<std::string> varNames;
	for (int i1=0; i1<numVars; i1++) {
		varNames.push_back(std::to_string(i1));
	}

	// Generate normalisation equations
	std::cout << "Generating normalisation equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<d; j++) {

			// Start with a -1 term
			std::vector<std::pair<int,char>> eqn = {{0, -1}};

			// Magnitude of this vector should sum to 1
			for (int k=0; k<d; k++) {

				// The variables to multiply
				int ind = i*d*d + j*d + k;

				// Add to the mapping if it's new
				std::string term = varNames[ind] + "," + varNames[ind+conjDelta] + ",";
				if (map.find(term) == map.end()) {
					map[term] = nextInd;
					invertedMap[nextInd] = term;
					nextInd++;
				}

				// Add this term to the equation
				eqn.push_back(std::pair<int,char>(map[term], 1));

			}

			// Add it to the list
			eqns.push_back(eqn);

		}
	}

	// Generate orthogonality equations
	std::cout << "Generating orthogonality equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<d; j++) {
			for (int k=j+1; k<d; k++) {

				// Start with two blank equations
				std::vector<std::pair<int,char>> eqn = {};
				std::vector<std::pair<int,char>> eqnConj = {};

				// The inner product of these two vectors should equal zero
				for (int l=0; l<d; l++) {

					// The variables to multiply
					int ind1 = i*d*d + j*d + l;
					int ind2 = i*d*d + k*d + l;

					// Add to the mapping if it's new
					std::string term1 = varNames[ind1] + "," + varNames[ind2+conjDelta] + ",";
					std::string term2 = varNames[ind1+conjDelta] + "," + varNames[ind2] + ",";
					if (map.find(term1) == map.end()) {
						map[term1] = nextInd;
						invertedMap[nextInd] = term1;
						nextInd++;
					}
					if (map.find(term2) == map.end()) {
						map[term2] = nextInd;
						invertedMap[nextInd] = term2;
						nextInd++;
					}

					// Add this term to the equation
					eqn.push_back(std::pair<int,char>(map[term1], 1));
					eqnConj.push_back(std::pair<int,char>(map[term2], 1));

				}

				// Add them to the list
				eqns.push_back(eqn);
				eqns.push_back(eqnConj);

			}
		}
	}

	// Generate MUB equations
	std::cout << "Generating MUB equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=i+1; j<n; j++) {

			// For each combination between these sets
			for (int k=0; k<d; k++) {
				for (int l=0; l<d; l++) {

					// Start with a -1 term (should be -1/d but we multiply by d for integers)
					std::vector<std::pair<int,char>> eqn = {{0, -1}};
					std::vector<std::pair<int,char>> eqnConj = {{0, -1}};

					// Multiply out the brackets
					for (int t1=0; t1<d; t1++) {
						for (int t2=0; t2<d; t2++) {

							// The variables to multiply
							int var1 = i*d*d + k*d + t1;
							int var2 = j*d*d + l*d + t1;
							int var3 = i*d*d + k*d + t2;
							int var4 = j*d*d + l*d + t2;

							// Add to the mapping if it's new
							std::string term1 = varNames[var3] + "," + varNames[var2]  + "," + varNames[var1+conjDelta] + "," + varNames[var4+conjDelta] + ",";
							std::string term2 = varNames[var1] + "," + varNames[var4]  + "," + varNames[var3+conjDelta] + "," + varNames[var2+conjDelta] + ",";
							if (map.find(term1) == map.end()) {
								map[term1] = nextInd;
								invertedMap[nextInd] = term1;
								nextInd++;
							}
							if (map.find(term2) == map.end()) {
								map[term2] = nextInd;
								invertedMap[nextInd] = term2;
								nextInd++;
							}

							// Add the term to the equation
							eqn.push_back(std::pair<int,char>(map[term1], d));
							eqnConj.push_back(std::pair<int,char>(map[term2], d));

						}
					}

					// Add them to the list
					eqns.push_back(eqn);
					eqns.push_back(eqnConj);

				}
			}

		}
	}

	// TODO combine
	
	// Get the polynomial to multiply by
	std::cout << "Generating moment vector..." << std::endl;
	std::vector<std::string> termsToMultiply = {""};
	if (k >= 1) {
		for (int i=0; i<numVars; i++) {
			termsToMultiply.push_back(std::to_string(i)+",");
		}
	}
	if (k >= 2) {
		for (int i=0; i<numVars; i++) {
			for (int j=0; j<numVars; j++) {
				termsToMultiply.push_back(std::to_string(i)+","+std::to_string(j)+",");
			}
		}
	}
	if (k >= 3) {
		for (int i=0; i<numVars; i++) {
			for (int j=0; j<numVars; j++) {
				for (int k=0; k<numVars; k++) {
					termsToMultiply.push_back(std::to_string(i)+","+std::to_string(j)+","+std::to_string(k)+",");
				}
			}
		}
	}
	if (k >= 4) {
		for (int i=0; i<numVars; i++) {
			for (int j=0; j<numVars; j++) {
				for (int k=0; k<numVars; k++) {
					for (int l=0; l<numVars; l++) {
						termsToMultiply.push_back(std::to_string(i)+","+std::to_string(j)+","+std::to_string(k)+","+std::to_string(l)+",");
					}
				}
			}
		}
	}

	// Multiply each original equation by each element of this polynomial 
	std::cout << "Multiplying equations..." << std::endl;
	int numOrig = eqns.size();
	std::vector<std::vector<std::pair<int,char>>> newEqns;
	for (int j=0; j<termsToMultiply.size(); j++) {
		std::cout << j << "/" << termsToMultiply.size() << std::endl;
		for (int i=0; i<numOrig; i++) {

			// Start with a copy
			std::vector<std::pair<int,char>> newEqn = eqns[i];

			// Multiply each term
			for (int k=0; k<newEqn.size(); k++) {

				// Convert back then multiply
				std::string term = multiplyFull(invertedMap[newEqn[k].first], termsToMultiply[j]);

				// Add to the mapping if it's new
				if (map.find(term) == map.end()) {
					map[term] = nextInd;
					invertedMap[nextInd] = term;
					nextInd++;
				}
				
				// Update this term
				newEqn[k].first = map[term];

			}

			// Add to the new array
			newEqns.push_back(newEqn);

		}
	}

	// Check everything TODO check for k>0
	//pretty(map);
	//std::cout << "before" << std::endl;
	//for (int i=0; i<eqns.size(); i++) {
		//pretty(invertedMap, eqns[i]);
	//}
	//std::cout << "after" << std::endl;
	//for (int i=0; i<newEqns.size(); i++) {
		//pretty(invertedMap, newEqns[i]);
	//}

	// Write to file
	std::cout << "Writing file..." << std::endl;
	std::string fileName = "matrices/d" + std::to_string(d) + "n" + std::to_string(n) + "k" + std::to_string(k) + ".csv"; 
	std::ofstream outFile;
	outFile.open(fileName);
	outFile << newEqns.size() << " " << nextInd << " " << eqns.size() << " " << numVars << std::endl;
	for (int i=0; i<newEqns.size(); i++) {
		std::cout << i << "/" << newEqns.size() << std::endl;
		for (int j=0; j<newEqns[i].size(); j++) {
			outFile << i << " " << newEqns[i][j].first << " " << int(newEqns[i][j].second) << std::endl;
		}
	}

}
	
