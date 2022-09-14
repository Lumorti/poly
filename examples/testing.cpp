#include "../poly.h"

std::vector<bool> intToBool(int val, int bits) {
	std::vector<bool> toReturn(bits, false);
	for (int i=0; i<toReturn.size(); i++) {
		toReturn[toReturn.size()-i-1] = ((val >> i) & 1);
	}
	return toReturn;
}

int boolToInt(std::vector<bool> val) {
	int toReturn = 0;
	for (int i=0; i<val.size(); i++) {
		toReturn += val[val.size()-1-i]*std::pow(2, i);
	}
	return toReturn;
}

bool f(std::vector<bool> x) {

	// Search O(numBits)
	//for (int i=0; i<x.size(); i++) {
		//if (x[i]) {
			//return true;
		//}
	//}
	//return false;

	// Primality testing O(numBits^6)
	int asInt = boolToInt(x);
	for (int i=2; i<std::ceil(std::sqrt(asInt)); i++) {
		if (asInt % i == 0) {
			return false;
		}
	}
	return true;

	// TODO addition or multiplication

}

// Standard cpp entry point 
int main(int argc, char ** argv) {

	int numBits = std::stoi(argv[1]);
	Polynomial<int> g(numBits);
	for (int i=0; i<std::pow(2, numBits); i++) {
		auto vec = intToBool(i, numBits);
		std::cout << vec << " " << f(vec) << std::endl;
		Polynomial<int> term(numBits, f(vec), {});
		for (int j=0; j<numBits; j++) {
			Polynomial<int> temp(numBits);
			if (vec[j]) {
				temp.addTerm(1, {numBits-j-1});
			} else {
				temp.addTerm(1, {});
				temp.addTerm(-1, {numBits-j-1});
			}
			term *= temp;
		}
		g += term;
	}
	g = g.prune();

	std::cout << g.asMathematica() << std::endl;
	std::cout << g.minimalHorner() << std::endl;

}
