#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Get the problem from the args
	int numSpins = 10;
	int verbosity = 1;
	int level = 1;
	int maxIters = -1;
	bool bruteForce = false;
	std::string filename = "";
	for (int i=0; i<argc; i++) {
		std::string arg = argv[i];
		if (arg == "-h") {
			std::cout << " -n [int]    set the number of spins" << std::endl;
			std::cout << " -v [int]    set the verbosity level" << std::endl;
			std::cout << " -l [int]    set the level of the branch and bound" << std::endl;
			std::cout << " -i [int]    set the maximum number of iterations" << std::endl;
			std::cout << " -f [str]    load the coefficients from a file" << std::endl;
			std::cout << " -b          brute force to find the optimum" << std::endl;
			std::cout << " -r          use a random seed" << std::endl;
			return 0;
		} else if (arg == "-n" && i+1 < argc) {
			numSpins = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-f" && i+1 < argc) {
			filename = argv[i+1];
			i++;
		} else if (arg == "-b") {
			bruteForce = true;
		} else if (arg == "-v" && i+1 < argc) {
			verbosity = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-l" && i+1 < argc) {
			level = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-i" && i+1 < argc) {
			maxIters = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-r") {
			std::srand(time(0));
		}

	}

	// If a filename is given, load the data from it
	std::vector<int> is;
	std::vector<int> js;
	std::vector<double> coeffs;
	if (filename != "") {

		// Open the file and check it's valid
		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cout << "Could not open file " << filename << std::endl;
			return 1;
		}

		// Get each line
		std::string line;
		while (std::getline(file, line)) {

			// Split into the i, j, and coefficient
			std::stringstream ss(line);
			int i, j;
			double coeff;
			ss >> i >> j >> coeff;

			// Determine the total number of spins
			numSpins = std::max(numSpins, std::max(i,j)+1);

			// Add to the vectors
			is.push_back(i);
			js.push_back(j);
			coeffs.push_back(coeff);

		}
	}

	// Create an Ising spin system with that many spins
	PolynomialProblem<double> prob(numSpins);
	for (int i=0; i<numSpins; i++) {
		prob.varIsBinary[i] = true;
		prob.varBounds[i] = {-1,1};
	}

	// If no file is given, generate a random problem
	if (filename == "") {
		for (int i=0; i<numSpins; i++) {
			for (int j=i+1; j<numSpins; j++) {
				prob.obj.addTerm(2*(double(rand())/RAND_MAX)-1, {i,j});
			}
		}

	// Otherwise use the coefficients from the file
	} else {
		for (int i=0; i<is.size(); i++) {
			prob.obj.addTerm(coeffs[i], {is[i], js[i]});
		}

	}

	// Brute force to find the optimum for comparison
	if (bruteForce) {
		auto res = prob.bruteForce();
		std::cout << "From brute forcing:" << std::endl;
		std::cout << res.first << " " << res.second << std::endl;
		std::cout << std::endl;
	}

	// Optimize using branch and bound
	auto res2 = prob.minimize(level, verbosity, maxIters);
	std::cout << "From branch and bound:" << std::endl;
	std::cout << res2.first << " " << res2.second << std::endl;
	std::cout << std::endl;

	return 0;

}
