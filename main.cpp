#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include <math.h>

// Algorithm settings
int maxLevel = 25;
double zeroVal = 1e-13;
int maxTries = 10;
long int maxIters = 1e10;
int useTest = 0;
bool outputCert = false;
double instabilityThresh = 1e5;

// Generate/add the vector of r-combinations of 1...n
void genCombinations(std::vector<std::vector<int>> &combs, int n, int r) {

	// https://stackoverflow.com/questions/9430568/generating-combinations-in-c
	std::vector<bool> v(n);
	std::fill(v.begin(), v.begin() + r, true);
	do {
		std::vector<int> comb;
		for (int i=0; i<n; ++i) {
            if (v[i]) {
				comb.push_back(i);
            }
        }
		combs.push_back(comb);
	} while (std::prev_permutation(v.begin(), v.end()));

}

// Print an unordered map
template<typename K, typename V>
void pretty(std::unordered_map<K, V> const &m) {
	std::cout << "{ ";
    for (auto const &pair: m) {
        std::cout << pair.first << ": " << pair.second << ", ";
    }
	std::cout << "} " << std::endl;;
}

// Remove zero terms from an equation
void cleanEqn(std::unordered_map<std::string,double> &a) {
	std::vector<std::string> toRem;
	for (auto& it: a) {
		if (std::abs(it.second) < zeroVal) {
			toRem.push_back(it.first);
		}
	}
	for (unsigned long int i=0; i<toRem.size(); i++) {
		a.erase(toRem[i]);
	}
}

// Add equation b to a in-place, scaling by c
void addEqns(std::unordered_map<std::string,double> &a, std::unordered_map<std::string,double> b, double c) {
	for (auto& it: b) {
		if (a.find(it.first) != a.end()) {
			a[it.first] += it.second * c;
		} else {
			a[it.first] = it.second * c;
		}
	}
}

// Multiply two terms together, sorting the result
// e.g. 1,3,7, * 2,5, = 1,2,3,5,7,
// this should be O(n+m)
std::string multiplyFull(std::string a, std::string b, int padding) {

	// If adding to nothing, just return
	if (a.size() == 0 || b.size() == 0) {
		return a + b;
	}

	// Start with a blank string, go through and add the other items
	std::string newString = "";

	// Loop through the first string
	unsigned long int j = 0;
	for (unsigned long int i=0; i<a.size(); i+=padding) {

		// Get the value
		std::string string1 = a.substr(i, padding);
		int val1 = std::stoi(string1);

		// Go through the other string
		for (j=j; j < b.size(); j+=padding) {

			// Get the value
			std::string string2 = b.substr(j, padding);
			int val2 = std::stoi(string2);

			// If this is less, add it
			if (val2 < val1) {
				newString += string2;
			} else {
				//j -= padding-1;
				break;
			}

		}

		// Then finally add the first string
		newString += string1;

	}

	// Add anything in b that bigger than all of a
	if (j < b.size()) {
		newString += b.substr(j, b.size()-j);
	}

	return newString;
}

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Process args
	std::string filename = "example.eqn";
	std::string arg = "";
	for (int i=1; i<argc; i++) {
		arg = argv[i];

		// Output the help
		if (arg == "-h") {
			std::cout << "" << std::endl;
			std::cout << "--------------------------------" << std::endl;
			std::cout << "             LMN                " << std::endl;
			std::cout << "  (Low Memory Nullstellensatz)  " << std::endl;
			std::cout << "         By Lumorti" << std::endl;
			std::cout << "--------------------------------" << std::endl;
			std::cout << "  Usage:" << std::endl;
			std::cout << "     ./lmn [args?] [equation file]" << std::endl;
			std::cout << "  Arguments:" << std::endl;
			std::cout << "     -e          use a test example" << std::endl;
			std::cout << "     -o [int]    set the max order (def 25)" << std::endl;
			std::cout << "     -z [float]  set the tolerance for zero (def 1e-13)" << std::endl;
			std::cout << "     -i [int]    set the max iterations (def 1e10)" << std::endl;
			std::cout << "     -t [int]    set the max tries (def 1000)" << std::endl;
			std::cout << "     -s [float]  set the threshold for instability" << std::endl;
			std::cout << "     -c          write the certificate linear system to a file" << std::endl;
			std::cout << "  Equation file format:" << std::endl;
			std::cout << "     x^2 + xy + y + 2 = 0 is represented as:" << std::endl;
			std::cout << "     0,0:1 + 0,1:1 + 1:1 + :2" << std::endl;
			std::cout << "     such that 0 is x and 1 is y." << std::endl;
			std::cout << "     (see the example file in the repo)" << std::endl;
			std::cout << "" << std::endl;
			return 0;

		// Set the max order
		} else if (arg == "-o") {
			maxLevel = std::stoi(argv[i+1]);
			i += 1;

		// Set the zero tolerance
		} else if (arg == "-z") {
			zeroVal = std::stod(argv[i+1]);
			i += 1;

		// Set the instability threshold
		} else if (arg == "-s") {
			instabilityThresh = std::stod(argv[i+1]);
			i += 1;

		// Set the max num of iterations
		} else if (arg == "-i") {
			maxIters = std::stoi(argv[i+1]);
			i += 1;

		// Set the max num of tries per iter
		} else if (arg == "-t") {
			maxTries = std::stoi(argv[i+1]);
			i += 1;

		// If told to output the certificate
		} else if (arg == "-c") {
			outputCert = true;

		// Otherwise assume it's a filename
		} else {
			filename = arg;

		}

	}

	// Open the file
	std::ifstream infile(filename);
	if (!infile) {
		std::cout << "ERROR - file does not exist" << std::endl;
		return 0;
	}

	// Process the file
	std::vector<std::unordered_map<std::string,double>> rawEqns;
	int numVars = 0;
	std::string line;
	std::string currentFirst = "";
	std::string currentSecond = "";
	bool onFirst = true;
	while (std::getline(infile, line)) {

		// Remove spaces from the line
		line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
		std::cout << line << std::endl;

		// Loop over the line
		currentFirst = "";
		currentSecond = "";
		onFirst = true;
		std::unordered_map<std::string,double> eqn;
		for (unsigned long int i=0; i<line.size(); i++) {

			// We're done with the key, now the value
			if (line[i] == ':') {
				onFirst = false;

				// Check all the var nums in this
				std::stringstream ss(currentFirst);
				for (int j; ss >> j;) {
					if (j+1 > numVars) {
						numVars = j+1;
					}
					if (ss.peek() == ',') {
						ss.ignore();
					}
				}

			// If it's a new term, add the previous
			} else if (line[i] == '+') {
				eqn[currentFirst] = std::stod(currentSecond);
				currentFirst = "";
				currentSecond = "";
				onFirst = true;

			// Add to the key
			} else if (onFirst) {
				currentFirst += line[i];

			// Add to the value
			} else {
				currentSecond += line[i];

			}

		}

		// Add whatever's left
		if (currentSecond.size() > 0) {
			eqn[currentFirst] = std::stod(currentSecond);
		}

		// Add this equation to the list
		rawEqns.push_back(eqn);

	}

	// Calculate the padding width needed
	int paddingWidth = std::log(numVars-1);

	// Fix the padding in all the equations
	int maxOrd = 0;
	std::vector<std::unordered_map<std::string,double>> eqns;
	for (unsigned long int i=0; i<rawEqns.size(); i++) {
		std::unordered_map<std::string,double> eqn;

		// For each term
		for (auto& it: rawEqns[i]) {
			std::string newString = "";

			// Check all the var nums in this
			std::stringstream ss(it.first);
			for (int j; ss >> j;) {

				// Pad it correctly
				std::string term = std::to_string(j);
				term.insert(term.begin(), paddingWidth-term.size(), ' ');
				newString += term; 

				// Skip the commas
				if (ss.peek() == ',') {
					ss.ignore();
				}

			}

			// Set the new key
			eqn[newString] = it.second;

			// What's the max order of this equation?
			int termOrd = int(it.first.size() / paddingWidth);
			if (termOrd > maxOrd) {
				maxOrd = termOrd;
			}

		}

		// Add to the corrected equation list
		eqns.push_back(eqn);

	}

	// Save a bit of memory by clearing the raw data array
	rawEqns.clear();

	// Start with 1
	std::unordered_map<std::string,double> mainEqn;
	mainEqn[""] = 1;

	// Keep track of which equations have been used for the certificate
	std::vector<std::unordered_map<std::string,double>> eqnsUsed;

	// List the equations available
	std::cout << numVars << std::endl;
	std::cout << paddingWidth << std::endl;
	std::cout << "Equations available:" << std::endl;
	for (unsigned long int i=0; i<eqns.size(); i++) {
		pretty(eqns[i]);
	}
	std::cout << "Starting equation:" << std::endl;
	pretty(mainEqn);

	// The list of equations that can be used to remove terms
	std::unordered_map<std::string,std::unordered_map<std::string,double>> eqnToRemoveTerm;
	std::string termToRemove = "";
	int termOrder = 0;

	// String of x's used to remove section of a string later
	std::string xString = "";
	for (int i=0; i<paddingWidth; i++) {
		xString += "x";
	}

	// Try to remove the lowest order terms first
	int ord;
	for (ord=1; ord<maxLevel; ord++) {

		// Stop if we've reached zero
		if (mainEqn.size() == 0) {
			break;
		}

		// Keep repeating for a while
		for (int iter=0; iter<maxIters; iter++) {

			// Stop if we've reached zero
			if (mainEqn.size() == 0) {
				break;
			}

			// Find the largest term in the main equation of this order (or less)
			termToRemove = "x";
			double bestMag = 0;
			double mag = 0;
			for (auto& it: mainEqn) {
				mag = std::abs(it.second);
				if (int(it.first.size() / paddingWidth) - 1 <= ord && mag > bestMag) {
					termToRemove = it.first;
					bestMag = mag;
				}
			}
			if (termToRemove == "x") {
				break;
			}
			termOrder = termToRemove.size() / paddingWidth - 1;
			std::cout << ord << " " << iter << "  in main = " << mainEqn.size() << " rem ord = " << termOrder+1 << std::endl;

			// Try to find something that could be multiplied to make this
			std::unordered_map<std::string,double> foundEqn;
			double bestCoeff = 0;
			std::string lookingFor = "";
			std::string toMultiply = "";
			int triesLeft = maxTries;
			for (int k=termOrder+1; k>=0; k--) {

				// Generate the combinations to check
				std::vector<std::vector<int>> combs;
				genCombinations(combs, termOrder+1, k);

				// Loop over these combinations
				for (unsigned long int m=0; m<combs.size(); m++) {

					// Extract the strings to search for
					lookingFor = "";
					toMultiply = termToRemove;
					for (unsigned long int l=0; l<combs[m].size(); l++) {
						lookingFor += termToRemove.substr(combs[m][l]*paddingWidth, paddingWidth);
					}
					for (unsigned long int l=0; l<combs[m].size(); l++) {
						toMultiply.replace(combs[m][l]*paddingWidth, paddingWidth, xString);
					}
					toMultiply.erase(std::remove(toMultiply.begin(), toMultiply.end(), 'x'), toMultiply.end());

					// See if this term is in one of the equations
					for (unsigned long int i=0; i<eqns.size(); i++) {
						std::unordered_map<std::string,double> checkEqn = eqns[i];

						// Multiply this by the factor needed to reach the term to remove
						std::unordered_map<std::string,double> mulEqn;
						for (auto& it: checkEqn) {
							std::string t = multiplyFull(it.first, toMultiply, paddingWidth);
							mulEqn[t] = it.second;
						}
						std::unordered_map<std::string,double> mulEqnBackup = mulEqn;

						// Reduce this using the equations that have already reduced the main
						bool somethingRemoved = true;
						while (somethingRemoved) {
							somethingRemoved = false;
							for (auto& it: mulEqn) {
								if (std::abs(it.second) > zeroVal && eqnToRemoveTerm.find(it.first) != eqnToRemoveTerm.end()) {
									addEqns(mulEqn, eqnToRemoveTerm[it.first], -mulEqn[it.first] / eqnToRemoveTerm[it.first][it.first]);
									somethingRemoved = true;
								}
							}
						}
						cleanEqn(mulEqn);

						// See if it's still valid after reduction
						if (mulEqn.find(termToRemove) != mulEqn.end()) {

							// If this is the equation with the biggest coeff for the term
							if (std::abs(mulEqn[termToRemove]) > bestCoeff) {
								foundEqn = mulEqn;
								bestCoeff = std::abs(mulEqn[termToRemove]);
							}

							// Only try a certain number of times before we assume we've found something big enough
							if (bestCoeff > 1e-3) {
								triesLeft--;
								if (triesLeft <= 0) {
									break;
								}
							}

						}

					}

					// Stop searching once we've tried enough times
					if (triesLeft <= 0) {
						break;
					}

				}

				// Stop searching once we've tried enough times
				if (triesLeft <= 0) {
					break;
				}

			}

			// Process the chosen row
			if (foundEqn.size() > 0) {

				// Apply this to all of the reduction equations
				for (auto& it: eqnToRemoveTerm) {
					if (it.second.find(termToRemove) != it.second.end()) {
						addEqns(it.second, foundEqn, -it.second[termToRemove] / foundEqn[termToRemove]);
						cleanEqn(it.second);
					}
				}

				// This can now be used to reduce other equations
				eqnToRemoveTerm[termToRemove] = foundEqn;

				// Reduce the main equations
				double factor = -mainEqn[termToRemove] / foundEqn[termToRemove];
				addEqns(mainEqn, foundEqn, -mainEqn[termToRemove] / foundEqn[termToRemove]);
				cleanEqn(mainEqn);

				// Stop if it's becoming unstable
				if (factor > instabilityThresh) {
					std::cout << "ERROR - system became unstable, try a higher try count" << std::endl;
					return 0;
				}

				// Keep track of the sums for the certificate
				if (outputCert) {
					eqnsUsed.push_back(foundEqn);
				}

			}

		}

	}

	// Final output
	std::cout << "Final equation:" << std::endl;
	pretty(mainEqn);

	// If we reached zero
	if (mainEqn.size() == 0) {
		std::cout << "Reached 0, meaning original system is inconsistent" << std::endl;
		std::cout << "Needed order " << ord  << " meaning full linalg would be " << std::pow(numVars, maxOrd+ord) << " x " << eqns.size()*std::pow(numVars, ord) << std::endl;
		std::cout << "Meanwhile this used a hash table with at most " << eqnToRemoveTerm.size()  << " elements" << std::endl;

		// If told to output the certificate linear system
		if (outputCert) {

			// Create the linear system using these
			std::cout << "Creating monomial mapping..." << std::endl;
			int newInd = 1;
			std::unordered_map<std::string,int> testMap;
			testMap[""] = 0;
			for (unsigned long int i=0; i<eqnsUsed.size(); i++) {
				for (auto& it: eqnsUsed[i]) {
					if (testMap.find(it.first) == testMap.end()) {
						testMap[it.first] = newInd;
						newInd++;
					}
				}
			}

			// Save the certificate system to file
			std::cout << "Writing file..." << std::endl;
			std::string fileName = "cert.csv"; 
			std::ofstream outFile;
			outFile.open(fileName);
			for (unsigned long int i=0; i<eqnsUsed.size(); i++) {
				for (auto& it: eqnsUsed[i]) {
					outFile << testMap[it.first] << " " << i << " " << it.second << std::endl;
				}
			}

		}

	// If we didn't
	} else {
		std::cout << "Didn't reach 0, meaning original system may be consistent" << std::endl;

	}

	return 0;

}
	
