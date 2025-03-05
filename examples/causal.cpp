#include "../poly.h"

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Get the problem from the args
    int d = 2;
	int verbosity = 1;
	int level = 1;
	int maxIters = -1;
	bool bruteForce = false;
	bool benchmark = false;
	std::string filename = "";
	for (int i=0; i<argc; i++) {
		std::string arg = argv[i];
		if (arg == "-h") {
			std::cout << " -d [int]    set the dimension" << std::endl;
			std::cout << " -v [int]    set the verbosity level" << std::endl;
			std::cout << " -l [int]    set the level of the branch and bound" << std::endl;
			std::cout << " -i [int]    set the maximum number of iterations" << std::endl;
			std::cout << " -r          use a random seed" << std::endl;
			return 0;
		} else if (arg == "-d" && i+1 < argc) {
			d = std::stoi(argv[i+1]);
			i++;
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

	// Generate the problem for GYNI TODO
    int widthW = std::pow(d, 4);
    int widthA = std::pow(d, 2);
    int widthB = std::pow(d, 2);
    int numVars = widthW*(widthW+1)/2 + widthA*(widthA+1)/2 + widthB*(widthB+1)/2;
	PolynomialProblem<double> prob(numVars);
    prob.obj = Polynomial<double>(numVars, 1, {1});

    // If verbose, output the problem
    if (verbosity > 0) {
        std::cout << "Problem:" << std::endl;
        std::cout << prob << std::endl;
    }

	// Optimize using branch and bound
	auto res2 = prob.optimize(level, verbosity, maxIters);
	if (!benchmark) {
		std::cout << std::endl;
		std::cout << "From branch and bound:" << std::endl;
		std::cout << res2.first << " " << res2.second << std::endl;
		std::cout << std::endl;
	}

	return 0;

}
