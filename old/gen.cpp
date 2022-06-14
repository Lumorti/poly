#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <complex>
#include <iomanip>
#include <random>
#include <fstream>
#include <chrono>
#include <math.h>

// Algorithm settings
int maxLevel = 25;
double zeroVal = 1e-13;
int maxTries = 1000;

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

// Sort and print an unordered map
template<typename K, typename V>
void prettySorted(std::unordered_map<K, V> &m) {
	std::vector<std::string> keys;
    for (auto const &pair: m) {
        keys.push_back(pair.first);
    }
	std::sort(keys.begin(), keys.end());
	std::cout << "{ " << std::endl;
    for (int i=0; i<keys.size(); i++) {
        std::cout << keys[i] << ": " << m[keys[i]] << ", " << std::endl;
    }
	std::cout << "} " << std::endl;;
}

// Print an equation
void pretty(std::unordered_map<int,std::string> invertedMap, std::vector<std::pair<int,char>> a) {
	std::cout << "[";
	for (int i=0; i<a.size(); i++) {
		std::cout << invertedMap[a[i].first] << ": " << int(a[i].second) << ", ";
	}
	std::cout << "]" << std::endl;
}

// Print an equation
std::string makePretty(std::unordered_map<std::string,double> eqn) {
	std::string ret = "[";
    for (auto const &pair: eqn) {
        ret += pair.first + ": " + std::to_string(pair.second) + ", ";
    }
	ret += "]";
	return ret;
}

// Print an equation
void pretty(std::vector<std::pair<std::string,double>> eqn) {
	std::cout << "[";
	for (int i=0; i<eqn.size(); i++) {
		std::cout << eqn[i].first << ": " << eqn[i].second << ", ";
	}
	std::cout << "]" << std::endl;
}

// Remove zero terms from an equation
void cleanEqn(std::unordered_map<std::string,double> &a) {
	std::vector<std::string> toRem;
	for (auto& it: a) {
		if (std::abs(it.second) < zeroVal) {
			toRem.push_back(it.first);
		}
	}
	for (int i=0; i<toRem.size(); i++) {
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

// Substitute and evaluate an equation
double eval(std::unordered_map<int,std::string> invertedMap, std::vector<double> ideals, std::vector<std::pair<int,char>> a, int padding) {

	// Add up each term in this equation
	double total = 0;
	for (int i=0; i<a.size(); i++) {

		// Start with the constant term
		double term = int(a[i].second);

		// Go through this 1,3,4
		std::string label = invertedMap[a[i].first];
		int firstStart = 0;
		for (int j=0; j<label.size(); j+=padding) {

			// Extract the index
			std::string string1 = label.substr(j, padding);
			int val1 = std::stoi(string1);

			// Multiply by this ideal value
			term *= ideals[val1];

		}

		// Add this term to the total
		total += term;

	}

	return total;

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
	int j = 0;
	for (int i=0; i<a.size(); i+=padding) {

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

	// Get the problem from the args
	int d = 2;
	int n = 4;
	int kVal = 0;
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}
	if (argc > 3) {
		kVal = std::stoi(argv[3]);
	}
	int numVarsNonConj = n*d*d;
	int numVars = 2*n*d*d;

	// The width of all of the digits in the string
	int paddingWidth = std::floor(std::log(numVars));

	// Init the monomial ordering
	int conjDelta = numVarsNonConj;
	int nextInd = 1;
	std::unordered_map<std::string, int> map = {{"", 0}};
	std::unordered_map<int, std::string> invertedMap;
	std::vector<std::vector<std::pair<int,char>>> eqns;
	std::vector<std::string> varNames;
	std::string name = "";
	for (int i1=0; i1<numVars; i1++) {
		name = std::to_string(i1);
		name.insert(name.begin(), paddingWidth - name.size(), ' ');
		varNames.push_back(name);
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
					std::string term1 = varNames[ind1] + varNames[ind2+conjDelta];
					std::string term2 = varNames[ind2] + varNames[ind1+conjDelta];
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
				std::string term = varNames[ind] + varNames[ind+conjDelta];
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
							std::string term1 = varNames[var3] + varNames[var2]  + varNames[var1+conjDelta] + varNames[var4+conjDelta];
							std::string term2 = varNames[var1] + varNames[var4]  + varNames[var3+conjDelta] + varNames[var2+conjDelta];
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

	//map.clear();
	//invertedMap.clear();
	//eqns.clear();
	//nextInd = 0;
	//numVars = 4;

	//std::string newTerm = "";

	//newTerm = "";
	//map[newTerm] = nextInd;
	//invertedMap[nextInd] = newTerm;
	//nextInd++;

	//newTerm = "  0";
	//map[newTerm] = nextInd;
	//invertedMap[nextInd] = newTerm;
	//nextInd++;

	//newTerm = "  1";
	//map[newTerm] = nextInd;
	//invertedMap[nextInd] = newTerm;
	//nextInd++;

	//newTerm = "  2";
	//map[newTerm] = nextInd;
	//invertedMap[nextInd] = newTerm;
	//nextInd++;

	//newTerm = "  3";
	//map[newTerm] = nextInd;
	//invertedMap[nextInd] = newTerm;
	//nextInd++;

	//newTerm = "  0  1";
	//map[newTerm] = nextInd;
	//invertedMap[nextInd] = newTerm;
	//nextInd++;

	//newTerm = "  2  3";
	//map[newTerm] = nextInd;
	//invertedMap[nextInd] = newTerm;
	//nextInd++;

	//eqns.push_back({{0, -1}, {1, 1}});
	//eqns.push_back({{0, -1}, {2, 1}});
	//eqns.push_back({{0, -1}, {3, 1}});
	//eqns.push_back({{0, -1}, {4, 1}});
	//eqns.push_back({{0, 1}, {5, -1}, {6, -1}});

	// Get the polynomial to multiply by
	//std::cout << "Generating moment vector..." << std::endl;
	//std::vector<std::string> termsToMultiply = {""};
	//if (kVal >= 1) {
		//for (int i=0; i<numVars; i++) {
			//termsToMultiply.push_back(varNames[i]);
		//}
	//}
	//if (kVal >= 2) {
		//for (int i=0; i<numVars; i++) {
			//for (int j=0; j<numVars; j++) {
				//termsToMultiply.push_back(varNames[i]+varNames[j]);
			//}
		//}
	//}
	//if (kVal >= 3) {
		//for (int i=0; i<numVars; i++) {
			//for (int j=0; j<numVars; j++) {
				//for (int k=0; k<numVars; k++) {
					//termsToMultiply.push_back(varNames[i]+varNames[j]+varNames[k]);
				//}
			//}
		//}
	//}
	//if (kVal >= 4) {
		//for (int i=0; i<numVars; i++) {
			//for (int j=0; j<numVars; j++) {
				//for (int k=0; k<numVars; k++) {
					//for (int l=0; l<numVars; l++) {
						//termsToMultiply.push_back(varNames[i]+varNames[j]+varNames[k]+varNames[l]);
					//}
				//}
			//}
		//}
	//}
	//if (kVal >= 5) {
		//for (int i=0; i<numVars; i++) {
			//for (int j=0; j<numVars; j++) {
				//for (int k=0; k<numVars; k++) {
					//for (int l=0; l<numVars; l++) {
						//for (int m=0; m<numVars; m++) {
							//termsToMultiply.push_back(varNames[i]+varNames[j]+varNames[k]+varNames[l]+varNames[m]);
						//}
					//}
				//}
			//}
		//}
	//}

	//// Multiply each original equation by each element of this polynomial 
	//std::cout << "Multiplying equations..." << std::endl;
	//int numOrig = eqns.size();
	//std::vector<std::vector<std::pair<int,char>>> newEqns;
	//for (int j=0; j<termsToMultiply.size(); j++) {
		////std::cout << "multiplying " << j << "/" << termsToMultiply.size() << std::endl;
		//for (int i=0; i<numOrig; i++) {

			//// Start with a copy
			//std::vector<std::pair<int,char>> newEqn = eqns[i];

			//// Multiply each term
			//for (int k=0; k<newEqn.size(); k++) {

				//// Convert back then multiply
				//std::string term = multiplyFull(invertedMap[newEqn[k].first], termsToMultiply[j], paddingWidth);

				//// Add to the mapping if it's new
				//if (map.find(term) == map.end()) {
					//map[term] = nextInd;
					//invertedMap[nextInd] = term;
					//nextInd++;
				//}
				
				//// Update this term
				//newEqn[k].first = map[term];

			//}

			//// Add to the new array
			//newEqns.push_back(newEqn);

		//}
	//}

	//// Check everything
	////pretty(map);
	////std::cout << "before" << std::endl;
	////for (int i=0; i<eqns.size(); i++) {
		////pretty(invertedMap, eqns[i]);
	////}
	////std::cout << "after" << std::endl;
	////double sqrt2 = 1.0 / std::sqrt(2.0);
	////std::vector<double> idealVals;
	////idealVals = {1, 0,   0, 1,   sqrt2, sqrt2,   sqrt2, -sqrt2, 1, 0,   0, 1,   sqrt2, sqrt2,   sqrt2, -sqrt2};
	////double maxEval = 0;
	////for (int i=0; i<newEqns.size(); i++) {
		////pretty(invertedMap, newEqns[i]);
		////double v = eval(invertedMap, idealVals, newEqns[i], paddingWidth);
		////std::cout << v << std::endl;
		////if (v > maxEval) {
			////maxEval = v;
		////}
	////}
	////std::cout << "max eval = " << maxEval << std::endl;

	//// Write to file
	//std::cout << "Writing file..." << std::endl;
	//std::string fileName = "matrices/d" + std::to_string(d) + "n" + std::to_string(n) + "k" + std::to_string(kVal) + ".csv"; 
	//std::ofstream outFile;
	//outFile.open(fileName);
	//for (int i=0; i<newEqns.size(); i++) {
		//std::cout << "writing " << i << "/" << newEqns.size() << std::endl;
		//for (int j=0; j<newEqns[i].size(); j++) {
			//outFile << i << " " << newEqns[i][j].first << " " << int(newEqns[i][j].second) << std::endl;
		//}
	//}

	//return 0;

	// A low-memory reduction method
	std::vector<std::unordered_map<std::string,double>> newnewEqns;
	std::unordered_map<std::string,double> mainEqn;

	// Test system
	if (argc == 1) {
		paddingWidth = 3;
		newnewEqns.push_back({{" 17", 1}, {"", -1}});
		newnewEqns.push_back({{"  0", 1}, {"", -1}});
		newnewEqns.push_back({{" 16", 1}, {"", -1}});
		newnewEqns.push_back({{"  1", 1}, {"", -1}});
		newnewEqns.push_back({{"  1 17", 1}, {"  0 16", 1}, {"", -1}});

	// Use the real deal
	} else {

		// Get the equations into a nice form
		for (int i=0; i<eqns.size(); i++) {
			std::unordered_map<std::string,double> newE;
			for (int j=0; j<eqns[i].size(); j++) {
				newE[invertedMap[eqns[i][j].first]] = double(eqns[i][j].second);
			}
			newnewEqns.push_back(newE);
		}

	}

	// Start with 1
	mainEqn[""] = 1;

	// Keep track of which equations have been used for the certificate
	std::vector<std::unordered_map<std::string,double>> eqnsUsed;

	// List the equations available
	std::cout << "equations available:" << std::endl;
	for (int i=0; i<newnewEqns.size(); i++) {
		pretty(newnewEqns[i]);
	}
	std::cout << "main equation:" << std::endl;
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
	for (int ord=1; ord<maxLevel; ord++) {

		// Keep repeating for a while
		for (int iter=0; iter<10000000000; iter++) {

			// Stop if we've reached zero
			if (mainEqn.size() == 0) {
				break;
			}

			// Find the largest term in the main equation of this order (or less)
			termToRemove = "x";
			double bestMag = 0;
			double maxVal = 0;
			double minVal = 10000;
			double mag = 0;
			for (auto& it: mainEqn) {
				mag = std::abs(it.second);
				if (int(it.first.size() / paddingWidth) - 1 <= ord && mag > bestMag) {
					termToRemove = it.first;
					bestMag = mag;
				}
				if (mag > maxVal) {
					maxVal = mag;
				}
				if (mag < minVal) {
					minVal = mag;
				}
			}
			if (termToRemove == "x") {
				break;
			}
			termOrder = termToRemove.size() / paddingWidth - 1;
			std::cout << ord << " " << iter << "  in main = " << mainEqn.size() << " rem ord = " << termOrder << " min/max = " << minVal << " " << maxVal << std::endl;

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
				for (int m=0; m<combs.size(); m++) {

					// Extract the strings to search for
					lookingFor = "";
					toMultiply = termToRemove;
					for (int l=0; l<combs[m].size(); l++) {
						lookingFor += termToRemove.substr(combs[m][l]*paddingWidth, paddingWidth);
					}
					for (int l=0; l<combs[m].size(); l++) {
						toMultiply.replace(combs[m][l]*paddingWidth, paddingWidth, xString);
					}
					toMultiply.erase(std::remove(toMultiply.begin(), toMultiply.end(), 'x'), toMultiply.end());
					//std::cout << "looking for '" << lookingFor  << "' (multiply by '" << toMultiply << "')" << std::endl;

					// See if this term is in one of the equations
					for (int i=0; i<newnewEqns.size(); i++) {
						std::unordered_map<std::string,double> checkEqn = newnewEqns[i];

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

							// If this is the equation with the biggest coeff for the term TODO
							if (std::abs(mulEqn[termToRemove]) > bestCoeff) {
								foundEqn = mulEqn;
								bestCoeff = std::abs(mulEqn[termToRemove]);
								if (bestCoeff > 1e-3) {
									triesLeft = 0;
								}
							}

							// Only try a certain number of times before we assume we've found something big enough
							triesLeft--;
							if (triesLeft <= 0) {
								break;
							}

						}

					}

					if (triesLeft <= 0) {
						break;
					}
					//if (foundEqn.size() > 0) {
						//break;
					//}

				}

				if (triesLeft <= 0) {
					break;
				}
				//if (foundEqn.size() > 0) {
					//break;
				//}

			}

			// Process the chosen row TODO
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

				//std::cout << "found with '" << lookingFor  << "' (multiply by '" << toMultiply << "')" << std::endl;
				//std::cout << "will now reduce '" << termToRemove  << "' with:" << std::endl;
				//pretty(mulEqn);

				// Reduce the main equations
				addEqns(mainEqn, foundEqn, -mainEqn[termToRemove] / foundEqn[termToRemove]);
				cleanEqn(mainEqn);

				// Keep track of the sums for the certificate
				eqnsUsed.push_back(foundEqn);

				// Output
				//std::cout << "new main equation:" << std::endl;
				//prettySorted(mainEqn);

			}

		}

	}

	std::cout << "final equation:" << std::endl;
	prettySorted(mainEqn);

	// If we reached zero
	if (mainEqn.size() == 0) {

		// Create the linear system using these
		std::cout << "Creating monomial mapping..." << std::endl;
		int newInd = 1;
		std::unordered_map<std::string,int> testMap;
		testMap[""] = 0;
		for (int i=0; i<eqnsUsed.size(); i++) {
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
		for (int i=0; i<eqnsUsed.size(); i++) {
			for (auto& it: eqnsUsed[i]) {
				outFile << testMap[it.first] << " " << i << " " << it.second << std::endl;
			}
		}

	}

	// New method:
	// 1) least squares optimise
	// 2) add equations based on residiuals
	// 3) repeat

	// Create the linear system using these
	//int newInd = 1;
	//std::unordered_map<std::string,int> symMap;
	//std::unordered_map<int,std::string> symMapInverted;
	//symMap[""] = 0;
	//std::vector<Eigen::Triplet<double>> coefficients;
	//for (int i=0; i<newnewEqns.size(); i++) {
		//for (auto& it: newnewEqns[i]) {
			//if (symMap.find(it.first) == symMap.end()) {
				//symMap[it.first] = newInd;
				//symMapInverted[newInd] = it.first;
				//newInd++;
			//}
			//coefficients.push_back(Eigen::Triplet<double>(symMap[it.first], i, it.second));
		//}
	//}
	//Eigen::SparseMatrix<double> A(newInd, newnewEqns.size());
	//A.setFromTriplets(coefficients.begin(), coefficients.end());
	//Eigen::SparseVector<double> b(A.rows());
	//b.coeffRef(0) = 1.0;

	//int iters = 20;
	//double tolerance = 1e-10;

	//std::string xString = "";
	//for (int i=0; i<paddingWidth; i++) {
		//xString += "x";
	//}

	//std::unordered_set<std::string> eqnsUsed;

	////Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
	//Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
	//solver.setMaxIterations(iters);
	//solver.setTolerance(tolerance);

	//double prevError = 10000;
	//double error = 1;
	//int prevCols = 1;
	//int prevRows = 1;
	//Eigen::VectorXd sol;
	//Eigen::VectorXd resids;
	//int origSize = b.size();
	//int prevResidZero = -10;
	//double prevScore = 10000000000000;

	//double zeroThresh = 1e-6;

	////for (int i=0; i<3; i++) {
	////for (int i=0; i<1000; i++) {
	//for (int i=0; i<10000000000000000; i++) {

		////Eigen::SparseMatrix<double> W(A.rows(), A.rows());
		////for (int j=0; j<W.rows(); j++) {
			////W.insert(j,j) = std::pow(W.rows()-j, 10);
			////W.insert(j,j) = std::pow(W.rows()-j, 10);
			////if (rand() > 0.1) {
				////W.insert(j,j) = 1;
			////} else {
				////W.insert(j,j) = 1000000;
			////}
		////}
		////W.makeCompressed();

		////std::cout << Eigen::MatrixXd(W) << std::endl;
		////std::cout << std::endl;
		////std::cout << Eigen::MatrixXd(A) << std::endl;
		////std::cout << std::endl;
		////std::cout << Eigen::MatrixXd(W*A) << std::endl;

		//// Minimise this system to get the min residuals
		//solver.compute(A);
		//Eigen::VectorXd sol = solver.solve(b);
		//Eigen::VectorXd resids = A*sol - b;
		//resids = resids.cwiseProduct(resids);
		//error = resids.norm();
		////std::cout << resids << std::endl;
		////std::cout << sol << std::endl;

		//int residZero = 0;
		//std::vector<int> resChoices;
		//for (int k=0; k<resids.size(); k++) {
			//if (std::abs(resids(k)) < zeroThresh) {
				//residZero += 1;
			//} else {
				//resChoices.push_back(k);
			//}
		//}

		//// Remove the latest column if nothing changed
		////if (error > prevError) {
		////double score = residZero - resids.size();
		//double score = error;
		//if (score > prevScore) {
			//for (int k=prevRows; k<b.size(); k++) {
				//std::string toRemove = symMapInverted[k];
				//symMap.erase(toRemove);
				//symMapInverted.erase(k);
			//}
			//A.conservativeResize(prevRows, prevCols);
			//b.conservativeResize(prevRows);
			//nextInd = prevRows;
		//} else {
			//prevScore = score;
		//}

		//std::cout << i << " err = " << error << " A!0 = " << A.nonZeros() << " mons = " << b.size() << " s = " << score << std::endl;

		//// Stop if the error becomes lower than some value
		//if (resids.norm() < zeroThresh) {
			//break;
		//}

		//// SCOREBOARD for 1000
		//// norm = 0.00605595 first = 0.0038096   iters = 20   nnz = 3195
		//// norm = 0.00541202 first = 0.00332837   iters = 20   nnz = 2851
		//// 294102 err = 0.00359397 A!0 = 4934 s = 0.00359397

		//// Get the biggest residuals
		////double smallestRes = -100;
		////for (int k=0; k<resids.size(); k++) {
			////if (-(smallestRes - resids(k)) > 1e-5) {
				////smallestRes = resids(k);
				////resChoices.clear();
				////resChoices.push_back(k);
			////} else if (std::abs(smallestRes-resids(k)) < 1e-5) {
				////resChoices.push_back(k);
			////}
		////}

		////for (int k=0; k<resChoices.size(); k++) {
			////std::cout << resChoices[k] << std::endl;
		////}

		////int ind = rand() % origSize;
		////if (i > 10000) {
			////ind = rand() % (origSize*origSize);
		////}
		////if (i > 10000) {
			////ind = rand() % (origSize*origSize);
		////}
		////int ind = rand() % newInd;
		//int ind1 = resChoices[rand()%resChoices.size()];
		//int ind2 = resChoices[rand()%resChoices.size()];
		////int ind = resChoices[0];
		////int ind = 3;
		////Eigen::Index t;
		////resids.maxCoeff(&t);
		////ind = t;
		//std::string biggestResid = symMapInverted[ind1];
		//std::string resid2 = symMapInverted[ind2];
		////std::cout << "ind chosen = " << ind << std::endl;
		////std::cout << "biggest resid = " << biggestResid << std::endl;

		//// Find an equation containing this variable
		//int termOrder = biggestResid.size() / paddingWidth - 1;
		//std::unordered_map<std::string,double> foundEqn;
		//std::string lookingFor = "";
		//std::string toMultiply = "";

		//// Try find the equation requiring the least amount of multiplication
		//for (int k=termOrder+1; k>=0; k--) {

			//// Generate the combinations to check
			//std::vector<std::vector<int>> combs;
			//genCombinations(combs, termOrder+1, k);

			//// Loop over these combinations
			//for (int m=0; m<combs.size(); m++) {

				//// Extract the strings to search for
				//lookingFor = "";
				//toMultiply = biggestResid;
				//for (int l=0; l<combs[m].size(); l++) {
					//lookingFor += biggestResid.substr(combs[m][l]*paddingWidth, paddingWidth);
				//}
				//for (int l=0; l<combs[m].size(); l++) {
					//toMultiply.replace(combs[m][l]*paddingWidth, paddingWidth, xString);
				//}
				//toMultiply.erase(std::remove(toMultiply.begin(), toMultiply.end(), 'x'), toMultiply.end());

				//// See if this term is in one of the equations
				//for (int i=0; i<newnewEqns.size(); i++) {
					//if (newnewEqns[i].find(lookingFor) != newnewEqns[i].end()) {

						//// Multiply this by the factor needed to reach the term to remove
						//std::unordered_map<std::string,double> checkEqn = newnewEqns[i];
						//std::unordered_map<std::string,double> mulEqn;
						//for (auto& it: checkEqn) {
							//std::string t = multiplyFull(it.first, toMultiply, paddingWidth);
							//mulEqn[t] = it.second;
						//}

						//int otherTerms = 0;
						//for (auto& it: mulEqn) {
							//if (symMap.find(it.first) != symMap.end()) {
								//if (std::abs(resids(symMap[it.first])) > zeroThresh) {
									//otherTerms++;
								//}
							//}
						//}
						//if (otherTerms == 0) {
							//continue;
						//}

						//// Ensure we don't use the same equation twice
						//std::string eqnName = toMultiply + "*" + std::to_string(i);
						//if (eqnsUsed.find(eqnName) == eqnsUsed.end()) {

							//// Stop this equation from being used again
							//eqnsUsed.insert(eqnName);

							//// Stop searching once we have something
							//foundEqn = mulEqn;
							//break;

						//}

					//}
				//}

				//// Stop searching once we have something
				//if (foundEqn.size() > 0) {
					//break;
				//}

			//}

			//// Stop searching once we have something
			//if (foundEqn.size() > 0) {
				//break;
			//}

		//}

		//// Add any new terms to the mappings
		//for (auto& it: foundEqn) {
			//if (symMap.find(it.first) == symMap.end()) {
				//symMap[it.first] = newInd;
				//symMapInverted[newInd] = it.first;
				//newInd++;
			//}
		//}

		//// Expand and add to A
		//prevRows = A.rows();
		//prevCols = A.cols();
		//A.conservativeResize(newInd, A.cols()+1);
		//for (auto& it: foundEqn) {
			//A.insert(symMap[it.first], A.cols()-1) = it.second;
		//}
		//A.makeCompressed();

		//// Expand b
		//b.conservativeResize(newInd);

	//}

	// Check what this would give if evaluated
	//std::vector<double> ids(numVars);
	//ids[0] = 1;
	//ids[1] = 0;
	//ids[16] = 1;
	//ids[17] = 0;
	//double tot = 0;
	//for (auto& it: mainEqn) {
		//std::string term = it.first;
		//double val = it.second;
		//for (int i=0; i<term.size(); i+=paddingWidth) {
			//val *= ids[std::stoi(term.substr(i, paddingWidth))];
		//}
		//tot += val;
	//}
	//std::cout << tot << std::endl;

	return 0;

}
	
