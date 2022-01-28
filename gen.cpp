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

// Print an equation
void pretty(std::unordered_map<int,std::string> invertedMap, std::vector<std::pair<int,char>> a) {
	std::cout << "[";
	for (int i=0; i<a.size(); i++) {
		std::cout << invertedMap[a[i].first] << ": " << int(a[i].second) << ", ";
	}
	std::cout << "]" << std::endl;
}

// Print an equation
void pretty(std::vector<std::pair<std::string,float>> eqn) {
	std::cout << "[";
	for (int i=0; i<eqn.size(); i++) {
		std::cout << eqn[i].first << ": " << eqn[i].second << ", ";
	}
	std::cout << "]" << std::endl;
}

// Add equation b to a in-place, scaling by c TODO
void addEqns(std::unordered_map<std::string,float> &a, std::unordered_map<std::string,float> b, float c) {
	for (auto& it: b) {
		if (a.find(it.first) != a.end()) {
			a[it.first] += it.second * c;
			if (std::abs(a[it.first]) < 1e-5) {
				a.erase(it.first);
			}
		} else {
			a[it.first] = it.second * c;
		}
	}
}

// Substitute and evaluate an equation
float eval(std::unordered_map<int,std::string> invertedMap, std::vector<float> ideals, std::vector<std::pair<int,char>> a, int padding) {

	// Add up each term in this equation
	float total = 0;
	for (int i=0; i<a.size(); i++) {

		// Start with the constant term
		float term = int(a[i].second);

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
	int firstStart = 0;
	int secondStart = 0;
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
	int n = 2;
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

	// Check everything
	//pretty(map);
	//std::cout << "before" << std::endl;
	//for (int i=0; i<eqns.size(); i++) {
		//pretty(invertedMap, eqns[i]);
	//}
	//std::cout << "after" << std::endl;
	//float sqrt2 = 1.0 / std::sqrt(2.0);
	//std::vector<float> idealVals;
	//idealVals = {1, 0,   0, 1,   sqrt2, sqrt2,   sqrt2, -sqrt2, 1, 0,   0, 1,   sqrt2, sqrt2,   sqrt2, -sqrt2};
	//float maxEval = 0;
	//for (int i=0; i<newEqns.size(); i++) {
		//pretty(invertedMap, newEqns[i]);
		//float v = eval(invertedMap, idealVals, newEqns[i], paddingWidth);
		//std::cout << v << std::endl;
		//if (v > maxEval) {
			//maxEval = v;
		//}
	//}
	//std::cout << "max eval = " << maxEval << std::endl;

	// TODO try low-memory reduction method
	std::vector<std::unordered_map<std::string,float>> newnewEqns;
	std::unordered_map<std::string,float> mainEqn;
	//for (int i=0; i<eqns.size(); i++) {
		//std::unordered_map<std::string,float> newE;
		//for (int j=0; j<eqns[i].size(); j++) {
			//newE[invertedMap[eqns[i][j].first]] = float(eqns[i][j].second);
		//}
		//newnewEqns.push_back(newE);
	//}

	newnewEqns.push_back({{"  0  1", 1}, {"  0", 1}, {"", 1}});
	newnewEqns.push_back({{"  0", 1}, {"", -1}});
	newnewEqns.push_back({{"  1", 1}, {"", -1}});

	// Use the negative of the last equation as a starting point
	mainEqn = newnewEqns[newnewEqns.size()-1];
	for (auto& it: mainEqn) {
		it.second = -it.second;
	}

	std::cout << "equations available:" << std::endl;
	for (int i=0; i<newnewEqns.size(); i++) {
		pretty(newnewEqns[i]);
	}
	std::cout << "main equation:" << std::endl;
	pretty(mainEqn);

	std::unordered_map<std::string,std::unordered_map<std::string,float>> eqnToRemoveTerm;
	eqnToRemoveTerm[""] = mainEqn;
	std::string termToRemove = "";
	int termOrder = 0;

	int maxLevel = 3;
	std::vector<std::vector<std::vector<int>>> combs(maxLevel);
	std::vector<int> blank;
	for (int i=0; i<maxLevel; i++) {
		for (int j=i+1; j>0; j--) {
			genCombinations(combs[i], i+1, j);
		}
		combs[i].push_back(blank);
	}

	std::vector<std::vector<std::vector<int>>> negCombs(combs.size());
	for (int l=0; l<combs.size(); l++) {
		std::vector<int> full;
		for (int i=0; i<l+1; i++) {
			full.push_back(i);
		}
		for (int i=0; i<combs[l].size(); i++) {
			std::vector<int> fullCopy = full;
			for (int j=0; j<combs[l][i].size(); j++) {
				fullCopy.erase(std::find(fullCopy.begin(), fullCopy.end(), combs[l][i][j]));
			}
			negCombs[l].push_back(fullCopy);
		}
	}

	for (int ord=0; ord<maxLevel; ord++) {

		for (int iter=0; iter<3; iter++) {

			termToRemove = "";
			for (auto& it: mainEqn) {
				if (it.first.size() > 0) {
					termToRemove = it.first;
					termOrder = termToRemove.size() / paddingWidth - 1;
					if (termOrder == ord) {
						break;
					}
				}
			}
			if (termToRemove.size() == 0) {
				break;
			}
			std::cout << "trying to remove " << termToRemove << " (order " << termOrder << ")" << std::endl;
			bool termRemoved = false;

			std::unordered_map<std::string,float> foundEqn;
			std::string lookingFor = "";
			std::string toMultiply = "";
			for (int k=0; k<combs[termOrder].size(); k++) {

				lookingFor = "";
				toMultiply = "";
				for (int l=0; l<combs[termOrder][k].size(); l++) {
					lookingFor += termToRemove.substr(combs[termOrder][k][l]*paddingWidth, paddingWidth);
				}
				for (int l=0; l<negCombs[termOrder][k].size(); l++) {
					toMultiply += termToRemove.substr(negCombs[termOrder][k][l]*paddingWidth, paddingWidth);
				}
				std::cout << "looking for '" << lookingFor  << "' (multiply by '" << toMultiply << "')" << std::endl;

				// See if this term is in one of the equations
				for (int i=0; i<newnewEqns.size(); i++) {
					if (newnewEqns[i].find(lookingFor) != newnewEqns[i].end()) {
						std::unordered_map<std::string,float> checkEqn = newnewEqns[i];

						std::cout << "found:" << std::endl;
						pretty(checkEqn);

						std::unordered_map<std::string,float> mulEqn;
						for (auto& it: checkEqn) {
							std::string t = multiplyFull(it.first, toMultiply, paddingWidth);
							mulEqn[t] = it.second;
						}

						// Reduce this using the equations that have already reduced the main
						for (auto& it: mulEqn) {
							if (eqnToRemoveTerm.find(it.first) != eqnToRemoveTerm.end()) {
								addEqns(mulEqn, eqnToRemoveTerm[it.first], -mulEqn[it.first] / eqnToRemoveTerm[it.first][it.first]);
							}
						}

						std::cout << "multiplied and reduced:" << std::endl;
						pretty(mulEqn);

						// See if it's still valid after reduction
						if (mulEqn.find(termToRemove) != mulEqn.end()) {

							eqnToRemoveTerm[termToRemove] = mulEqn;
							addEqns(mainEqn, mulEqn, -mainEqn[termToRemove] / mulEqn[termToRemove]);

							std::cout << "new main equation:" << std::endl;
							pretty(mainEqn);
							
							termRemoved = true;
							break;

						}

					}

				}

				if (termRemoved) {
					break;
				}

			}

			if (!termRemoved) {
				std::cout << "couldn't solve system" << std::endl;
				return 0;
			}

		}

	}

	// Write to file
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

	return 0;

}
	
