#include "../poly.h"

#include <vector>
#include <stdexcept>
#include <iostream>

// Addition of two matrices of polynomials
template <typename T>
std::vector<std::vector<Polynomial<T>>> operator+(const std::vector<std::vector<Polynomial<T>>>& A, const std::vector<std::vector<Polynomial<T>>>& B) {
    std::vector<std::vector<Polynomial<T>>> C(A.size(), std::vector<Polynomial<T>>(A[0].size(), Polynomial<T>(A[0][0].maxVariables)));
    for (int i=0; i<A.size(); i++) {
        for (int j=0; j<A[0].size(); j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

// Generic addition of two matrices
template <typename T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
    std::vector<std::vector<T>> C(A.size(), std::vector<T>(A[0].size()));
    for (int i=0; i<A.size(); i++) {
        for (int j=0; j<A[0].size(); j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

// Subtraction of two matrices of polynomials
template <typename T>
std::vector<std::vector<Polynomial<T>>> operator-(const std::vector<std::vector<Polynomial<T>>>& A, const std::vector<std::vector<Polynomial<T>>>& B) {
    std::vector<std::vector<Polynomial<T>>> C(A.size(), std::vector<Polynomial<T>>(A[0].size(), Polynomial<T>(A[0][0].maxVariables)));
    for (int i=0; i<A.size(); i++) {
        for (int j=0; j<A[0].size(); j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

// Trace of a polynomial matrix
template <typename T>
Polynomial<T> trace(const std::vector<std::vector<Polynomial<T>>>& A) {
    Polynomial<T> tr(A[0][0].maxVariables);
    for (int i=0; i<A.size(); i++) {
        tr += A[i][i];
    }
    return tr;
}

// Generic trace of a matrix
template <typename T>
T trace(const std::vector<std::vector<T>>& A) {
    T tr = T(0);
    for (int i=0; i<A.size(); i++) {
        tr += A[i][i];
    }
    return tr;
}

// Matrix product of two matrices of polynomials
template <typename T>
std::vector<std::vector<Polynomial<T>>> operator*(const std::vector<std::vector<Polynomial<T>>>& A, const std::vector<std::vector<Polynomial<T>>>& B) {
    std::vector<std::vector<Polynomial<T>>> C(A.size(), std::vector<Polynomial<T>>(B[0].size(), Polynomial<T>(A[0][0].maxVariables)));
    for (int i=0; i<A.size(); i++) {
        for (int j=0; j<B[0].size(); j++) {
            for (int k=0; k<A[0].size(); k++) {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return C;
}

// Generic matrix product of two matrices
template <typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
    std::vector<std::vector<T>> C(A.size(), std::vector<T>(B[0].size()));
    for (int i=0; i<A.size(); i++) {
        for (int j=0; j<B[0].size(); j++) {
            for (int k=0; k<A[0].size(); k++) {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return C;
}

// Hadamard product of two matrices of polynomials
template <typename T>
std::vector<std::vector<Polynomial<T>>> hadamard(const std::vector<std::vector<Polynomial<T>>>& A, const std::vector<std::vector<Polynomial<T>>>& B) {
    std::vector<std::vector<Polynomial<T>>> C(A.size(), std::vector<Polynomial<T>>(A[0].size(), Polynomial<T>(A[0][0].maxVariables)));
    for (int i=0; i<A.size(); i++) {
        for (int j=0; j<A[0].size(); j++) {
            C[i][j] = A[i][j]*B[i][j];
        }
    }
    return C;
}

// Tensor product of two matrices of polynomials
template <typename T>
std::vector<std::vector<Polynomial<T>>> tensor(const std::vector<std::vector<Polynomial<T>>>& A, const std::vector<std::vector<Polynomial<T>>>& B) {
    std::vector<std::vector<Polynomial<T>>> C(A.size()*B.size(), std::vector<Polynomial<T>>(A[0].size()*B[0].size(), Polynomial<T>(A[0][0].maxVariables)));
    for (int i=0; i<A.size(); i++) {
        for (int j=0; j<A[0].size(); j++) {
            for (int k=0; k<B.size(); k++) {
                for (int l=0; l<B[0].size(); l++) {
                    C[i*B.size()+k][j*B[0].size()+l] = A[i][j]*B[k][l];
                }
            }
        }
    }
    return C;
}

// Generic tensor product of two matrices
template <typename T>
std::vector<std::vector<T>> tensor(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
    std::vector<std::vector<T>> C(A.size()*B.size(), std::vector<T>(A[0].size()*B[0].size()));
    for (int i=0; i<A.size(); i++) {
        for (int j=0; j<A[0].size(); j++) {
            for (int k=0; k<B.size(); k++) {
                for (int l=0; l<B[0].size(); l++) {
                    C[i*B.size()+k][j*B[0].size()+l] = A[i][j]*B[k][l];
                }
            }
        }
    }
    return C;
}

// Given a complex matrix that should be positive, get the real version
// i.e. M >= 0
// =>
// ( M_r -M_i)
// ( M_i M_r ) >= 0
std::vector<std::vector<Polynomial<double>>> complexSDPToReal(const std::vector<std::vector<Polynomial<std::complex<double>>>>& M) {

    // The result is a 2n x 2n matrix
    int n = M.size();
    std::vector<std::vector<Polynomial<double>>> result(2*n, std::vector<Polynomial<double>>(2*n, Polynomial<double>(M[0][0].maxVariables)));

    // Fill the result matrix
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            result[i][j] = real(M[i][j]);
            result[i+n][j+n] = real(M[i][j]);
            result[i+n][j] = imag(M[i][j]);
            result[i][j+n] = -imag(M[i][j]);
        }
    }

    return result;
}

/**
 * Convert an integer idx in [0, d^n - 1] to its base-d representation
 * as an n-element array. The least significant "digit" is stored
 * at result[n-1], i.e. result is in "big-endian" from left to right.
 */
std::vector<int> indexToMultiIndex(int idx, int n, int d) {
    std::vector<int> result(n);
    for(int i = n - 1; i >= 0; --i) {
        result[i] = idx % d;
        idx /= d;
    }
    return result;
}

/**
 * Convert an n-element array of base-d "digits" to a single integer
 * in [0, d^n - 1]. The input is assumed to be a valid representation
 * in base d. 
 */
int multiIndexToIndex(const std::vector<int>& digits, int d) {
    int idx = 0;
    for(int val : digits) {
        idx = idx * d + val;
    }
    return idx;
}

/**
 * Insert an integer x in [0, d-1] at position 'subsystem' 
 * into an (n-1)-element base-d array, producing an n-element base-d array.
 */
std::vector<int> insertSubsystemIndex(const std::vector<int>& shortIndex, 
                                      int subsystem, int x) {
    // shortIndex has length (n-1), new array has length n
    // We insert x at position 'subsystem'.
    std::vector<int> result;
    result.reserve(shortIndex.size() + 1);

    // copy up to 'subsystem'
    for(int i = 0; i < subsystem; ++i) {
        result.push_back(shortIndex[i]);
    }

    // insert x
    result.push_back(x);

    // copy the remainder
    for(int i = subsystem; i < (int)shortIndex.size(); ++i) {
        result.push_back(shortIndex[i]);
    }

    return result;
}


/**
 * Convert an integer idx in [0, d^(n-1) - 1] to its base-d representation
 * of length (n-1), then insert a digit x at position 'subsystem'
 * to get an n-element representation. Finally convert that to an 
 * integer in [0, d^n - 1].
 *
 * This effectively gives the row/column index in the full matrix
 * if the partial-trace row/column index is idx, after inserting
 * the subsystem index x.
 */
int embedIndexWithSubsystem(int idx, int n, int d, int subsystem, int x) {
    // n -> dimension of the full system is d^n
    // We are dealing with (n-1)-subsystem indices => dimension is d^(n-1)
    // Convert idx -> (n-1)-digit base-d
    std::vector<int> shortIndex(n-1);
    {
        int tmp = idx;
        for(int i = (n-1) - 1; i >= 0; --i) {
            shortIndex[i] = tmp % d;
            tmp /= d;
        }
    }

    // Insert x at the 'subsystem'-th position
    std::vector<int> fullIndex = insertSubsystemIndex(shortIndex, subsystem, x);

    // Convert back to integer
    return multiIndexToIndex(fullIndex, d);
}


/**
 * Compute the partial trace of a d^n x d^n matrix over subsystem 'subsystemToTraceOut'.
 *
 * Template parameter T is the type of the matrix elements, e.g. double, std::complex<double>, etc.
 *
 * \param mat              The input matrix: size d^n x d^n.
 * \param n                Number of subsystems.
 * \param d                Dimension of each subsystem.
 * \param subsystemToTraceOut  Which subsystem (0-based) to trace out.
 * \return                 A matrix (d^(n-1) x d^(n-1)) that is the partial trace.
 *
 * Example usage:
 *     std::vector<std::vector<double>> M = ...  // d^n x d^n
 *     auto partial = partialTrace(M, 3, 2, 1);  // trace out subsystem 1 in a 2^3 system
 */
template <typename T>
std::vector<std::vector<Polynomial<T>>> partialTrace(
    const std::vector<std::vector<Polynomial<T>>>& mat,
    int n,
    int d,
    int subsystemToTraceOut
) {
    // Check the input size
    const int dimFull = 1;
    for(int k = 0; k < n; ++k) {
        // This naive approach calculates d^n by repeated multiplication
        // but in practice you might want to check for overflow, etc.
    }
    // Actually compute d^n
    int dimCheck = 1;
    for(int k = 0; k < n; ++k) {
        dimCheck *= d;
    }
    if((int)mat.size() != dimCheck || (int)mat[0].size() != dimCheck) {
        throw std::invalid_argument("Matrix must be d^n x d^n in partialTrace");
    }

    // The result is a d^(n-1) x d^(n-1) matrix
    int dimPartial = 1;
    for(int k = 0; k < n - 1; ++k) {
        dimPartial *= d;
    }
    std::vector<std::vector<Polynomial<T>>> result(dimPartial, std::vector<Polynomial<T>>(dimPartial, Polynomial<T>(mat[0][0].maxVariables)));

    // For each row/col in the d^(n-1) x d^(n-1) partial trace,
    // we sum over x in [0, d-1], i.e. the index that was "traced out".
    for(int rowP = 0; rowP < dimPartial; ++rowP) {
        for(int colP = 0; colP < dimPartial; ++colP) {

            Polynomial<T> sumVal = Polynomial<T>(mat[0][0].maxVariables);

            // Sum over all x (the basis states of the traced-out subsystem)
            for(int x = 0; x < d; ++x) {
                // Convert rowP to a row index in the full space, embedding x 
                int fullRow = embedIndexWithSubsystem(rowP, n, d, subsystemToTraceOut, x);
                // Convert colP to a col index in the full space, embedding x 
                int fullCol = embedIndexWithSubsystem(colP, n, d, subsystemToTraceOut, x);

                sumVal += mat[fullRow][fullCol];
            }

            result[rowP][colP] = sumVal;
        }
    }

    return result;
}

/**
 * Reorder the subsystems of a d^n x d^n matrix according to a permutation p.
 *
 * - n: number of subsystems
 * - d: dimension of each subsystem
 * - p: a permutation of {0, 1, ..., n-1} describing how to reorder.
 *
 * The new matrix M' is such that M'(row', col') = M(row, col),
 * where row'/col' are the base-d digits (subsystem indices) rearranged
 * by p from row/col.
 *
 * Example:
 *   // Original system: A (x) B (x) C => p = {0,1,2} identity)
 *   // Suppose we want B (x) A (x) C => p = {1,0,2}
 *   // We reorder a 2^3 x 2^3 matrix from ABC-basis to BAC-basis.
 */
template <typename T>
std::vector<std::vector<Polynomial<T>>> reorderSubsystems(
    const std::vector<std::vector<Polynomial<T>>>& mat,
    int n,
    int d,
    const std::vector<int>& p
) {
    // 1) Check that p is a valid permutation of size n
    if((int)p.size() != n) {
        throw std::invalid_argument("Permutation p must have size n.");
    }
    // Check that each integer in 0..n-1 appears exactly once
    {
        std::vector<int> check = p;
        std::sort(check.begin(), check.end());
        for(int i = 0; i < n; ++i) {
            if(check[i] != i) {
                throw std::invalid_argument("Permutation p must contain 0..n-1 exactly once.");
            }
        }
    }

    // 2) Check the input matrix dimension: must be d^n x d^n
    int dimFull = 1;
    for(int k = 0; k < n; ++k) {
        dimFull *= d;
    }
    if((int)mat.size() != dimFull || (int)mat[0].size() != dimFull) {
        throw std::invalid_argument("Matrix must be d^n x d^n in reorderSubsystems.");
    }

    // 3) Allocate the new (reordered) matrix
    std::vector<std::vector<Polynomial<T>>> newMat(dimFull, std::vector<Polynomial<T>>(dimFull, Polynomial<T>(mat[0][0].maxVariables)));

    // 4) For each row/col in the old matrix, determine where it goes in the new matrix
    for(int oldRow = 0; oldRow < dimFull; ++oldRow) {
        // decode the base-d digits for oldRow
        std::vector<int> oldRowDigits = indexToMultiIndex(oldRow, n, d);

        // reorder the digits: newRowDigits[i] = oldRowDigits[p[i]]
        std::vector<int> newRowDigits(n);
        for(int i = 0; i < n; ++i) {
            newRowDigits[i] = oldRowDigits[p[i]];
            //newRowDigits[p[i]] = oldRowDigits[i];
        }
        // convert back to integer
        int newRow = multiIndexToIndex(newRowDigits, d);

        // Do the same for columns
        for(int oldCol = 0; oldCol < dimFull; ++oldCol) {
            std::vector<int> oldColDigits = indexToMultiIndex(oldCol, n, d);
            std::vector<int> newColDigits(n);
            for(int i = 0; i < n; ++i) {
                newColDigits[i] = oldColDigits[p[i]];
            }
            int newCol = multiIndexToIndex(newColDigits, d);

            // 5) The reordered matrix’s (newRow, newCol) is the old matrix’s (oldRow, oldCol)
            newMat[newRow][newCol] = mat[oldRow][oldCol];
        }
    }

    return newMat;
}

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Get the problem from the args
    int d = 2;
	int verbosity = 1;
	int level = 1;
	int maxIters = -1;
    std::string problemName = "gyni";
    std::string ineqString = "";
    bool outputFinal = false;
    bool useKnown = false;
    bool useMatProds = false;
    bool useAnsatz = false;
    bool bEqualsA = false;
    bool fixedAB = false;
    bool seesaw = false;
    bool justPrintHessian = false;
	for (int i=0; i<argc; i++) {
		std::string arg = argv[i];
		if (arg == "-h" || arg == "--help") {
			std::cout << " -d [int]    set the dimension" << std::endl;
			std::cout << " -v [int]    set the verbosity level" << std::endl;
			std::cout << " -l [int]    set the level of the branch and bound" << std::endl;
			std::cout << " -i [int]    set the maximum number of iterations" << std::endl;
			std::cout << " --gyni      use the GYNI problem" << std::endl;
			std::cout << " --lgyni     use the LGYNI problem" << std::endl;
			std::cout << " --rlgyni    use the RLGYNI problem" << std::endl;
			std::cout << " -m          include products of matrices" << std::endl;
			std::cout << " -o          output the final inequality" << std::endl;
			std::cout << " -s          optimize using seesaw" << std::endl;
			std::cout << " -H          just print the Hessian of the objective" << std::endl;
			std::cout << " -a          use ansatz approach" << std::endl;
			std::cout << " -B          set all A and B matrices to be equal" << std::endl;
			std::cout << " -f          use fixed A and B matrices" << std::endl;
			std::cout << " -k          test with a known solution" << std::endl;
			std::cout << " -r          use a random seed" << std::endl;
			return 0;
		} else if (arg == "-d" && i+1 < argc) {
			d = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-v" && i+1 < argc) {
			verbosity = std::stoi(argv[i+1]);
			i++;
        } else if (arg == "-H") {
            justPrintHessian = true;
        } else if (arg == "-o") {
            outputFinal = true;
        } else if (arg == "--gyni") {
            problemName = "gyni";
        } else if (arg == "--lgyni") {
            problemName = "lgyni";
        } else if (arg == "--rlgyni") {
            problemName = "rlgyni";
        } else if (arg == "-m") {
            useMatProds = true;
        } else if (arg == "-s") {
            seesaw = true;
        } else if (arg == "-B") {
            bEqualsA = true;
        } else if (arg == "-f") {
            fixedAB = true;
        } else if (arg == "-a") {
            useAnsatz = true;
		} else if (arg == "-l" && i+1 < argc) {
			level = std::stoi(argv[i+1]);
			i++;
        } else if (arg == "-k") {
            useKnown = true;
		} else if (arg == "-i" && i+1 < argc) {
			maxIters = std::stoi(argv[i+1]);
			i++;
		} else if (arg == "-r") {
			std::srand(time(0));
		}

	}

	// Generate the problem for GYNI
    int d2 = d*d;
    int d3 = d*d*d;
    int d4 = d*d*d*d;
    int widthW = d4;
    int numAs = 4;
    int numBs = 4;
    int widthA = d2;
    int widthB = d2;
    int fullWidth = widthW + numAs*widthA + numBs*widthB;
    int numVars = 2*(widthW*(widthW+1)/2 + numAs*widthA*(widthA+1)/2 + numBs*widthB*(widthB+1)/2);
	PolynomialProblem<double> prob(numVars);

    // Output some sizes
    if (verbosity > 0) {
        std::cout << "Full (complex) matrix size: " << fullWidth << std::endl;
        std::cout << "Number of variables: " << numVars << std::endl;
        std::cout << "Width of W: " << widthW << std::endl;
        std::cout << "Width of A: " << widthA << std::endl;
        std::cout << "Width of B: " << widthB << std::endl;
    }

    // Known valid solution
    if (verbosity > 0) {
        std::cout << "Constructing known solution..." << std::endl;
    }
    std::vector<std::vector<std::complex<double>>> knownW(widthW, std::vector<std::complex<double>>(widthW, 0));
    std::vector<std::vector<std::vector<std::complex<double>>>> knownAs;
    for (int i=0; i<numAs; i++) {
        std::vector<std::vector<std::complex<double>>> A(widthA, std::vector<std::complex<double>>(widthA, 0));
        knownAs.push_back(A);
    }
    std::vector<std::vector<std::vector<std::complex<double>>>> knownBs;
    for (int i=0; i<numBs; i++) {
        std::vector<std::vector<std::complex<double>>> B(widthB, std::vector<std::complex<double>>(widthB, 0));
        knownBs.push_back(B);
    }

    // Define the matrices
    std::vector<std::vector<std::complex<double>>> Z(2, std::vector<std::complex<double>>(2, 0));
    Z[0][0] = 1;
    Z[1][1] = -1;
    std::vector<std::vector<std::complex<double>>> Y(2, std::vector<std::complex<double>>(2, 0));
    Y[0][1] = -1i;
    Y[1][0] = 1i;
    std::vector<std::vector<std::complex<double>>> X(2, std::vector<std::complex<double>>(2, 0));
    X[0][1] = 1;
    X[1][0] = 1;
    std::vector<std::vector<std::complex<double>>> I2(2, std::vector<std::complex<double>>(2, 0));
    I2[0][0] = 1;
    I2[1][1] = 1;
    std::vector<std::vector<std::complex<double>>> proj0(2, std::vector<std::complex<double>>(2, 0));
    proj0[0][0] = 1;
    std::vector<std::vector<std::complex<double>>> proj1(2, std::vector<std::complex<double>>(2, 0));
    proj1[1][1] = 1;
    std::vector<std::vector<std::complex<double>>> ZZZI = tensor(Z, tensor(Z, tensor(Z, I2)));
    std::vector<std::vector<std::complex<double>>> ZIXX = tensor(Z, tensor(I2, tensor(X, X)));

    // W = 0.25 * (I_8 + (1/sqrt(2)) * (ZZZI + ZIXX))
    if (d == 2) {
        for (int i=0; i<widthW; i++) {
            knownW[i][i] = 1;
        }
        for (int i=0; i<widthW; i++) {
            for (int j=i; j<widthW; j++) {
                knownW[i][j] += (1/sqrt(2)) * (ZZZI[i][j] + ZIXX[i][j]);
            }
        }
        for (int i=0; i<widthW; i++) {
            for (int j=0; j<widthW; j++) {
                knownW[i][j] *= 0.25;
            }
        }
        for (int i=0; i<widthW; i++) {
            for (int j=i; j<widthW; j++) {
                knownW[j][i] = std::conj(knownW[i][j]);
            }
        }

        // A_0_0 = 0
        // A_1_0 = |00><00| + |00><11| + |11><00| + |11><11|
        knownAs[2][0][0] = 1;
        knownAs[2][0][3] = 1;
        knownAs[2][3][0] = 1;
        knownAs[2][3][3] = 1;
        
        // A_0_1 = |0><0| (x) |0><0|
        // A_1_1 = |1><1| (x) |0><0|
        knownAs[1][0][0] = 1;
        knownAs[3][2][2] = 1;

        // Bs the same as the As
        for (int i=0; i<numBs; i++) {
            knownBs[i] = knownAs[i];
        }

    }

    // Ansats approach
    std::vector<std::vector<Polynomial<std::complex<double>>>> W(widthW, std::vector<Polynomial<std::complex<double>>>(widthW, Polynomial<std::complex<double>>(numVars, 0i)));
    std::vector<std::vector<std::vector<Polynomial<std::complex<double>>>>> As;
    std::vector<std::vector<std::vector<Polynomial<std::complex<double>>>>> Bs;
    if (useAnsatz) {
        std::vector<std::vector<std::complex<double>>> IIII = tensor(I2, tensor(I2, tensor(I2, I2)));
        std::vector<std::vector<std::complex<double>>> ZIZI = tensor(Z, tensor(I2, tensor(Z, I2)));
        std::vector<std::vector<std::complex<double>>> ZIII = tensor(Z, tensor(I2, tensor(I2, I2)));
        std::vector<std::vector<std::complex<double>>> IIZI = tensor(I2, tensor(I2, tensor(Z, I2)));
        std::vector<std::vector<std::complex<double>>> ZIIZ = tensor(Z, tensor(I2, tensor(I2, Z)));
        std::vector<std::vector<std::complex<double>>> IZZI = tensor(I2, tensor(Z, tensor(Z, I2)));
        std::vector<std::vector<std::complex<double>>> ZIZZ = tensor(Z, tensor(I2, tensor(Z, Z)));
        std::vector<std::vector<std::complex<double>>> ZZZI = tensor(Z, tensor(Z, tensor(Z, I2)));
        std::vector<std::vector<std::complex<double>>> ZIXX = tensor(Z, tensor(I2, tensor(X, X)));
        std::vector<std::vector<std::complex<double>>> ZIYY = tensor(Z, tensor(I2, tensor(Y, Y)));
        std::vector<std::vector<std::complex<double>>> XXZI = tensor(X, tensor(X, tensor(Z, I2)));
        std::vector<std::vector<std::complex<double>>> YYZI = tensor(Y, tensor(Y, tensor(Z, I2)));
        std::vector<std::vector<std::complex<double>>> II = tensor(I2, I2);
        std::vector<std::vector<std::complex<double>>> IX = tensor(I2, X);
        std::vector<std::vector<std::complex<double>>> IY = tensor(I2, Y);
        std::vector<std::vector<std::complex<double>>> IZ = tensor(I2, Z);
        std::vector<std::vector<std::complex<double>>> XI = tensor(X, I2);
        std::vector<std::vector<std::complex<double>>> XX = tensor(X, X);
        std::vector<std::vector<std::complex<double>>> XY = tensor(X, Y);
        std::vector<std::vector<std::complex<double>>> XZ = tensor(X, Z);
        std::vector<std::vector<std::complex<double>>> YI = tensor(Y, I2);
        std::vector<std::vector<std::complex<double>>> YX = tensor(Y, X);
        std::vector<std::vector<std::complex<double>>> YY = tensor(Y, Y);
        std::vector<std::vector<std::complex<double>>> YZ = tensor(Y, Z);
        std::vector<std::vector<std::complex<double>>> ZI = tensor(Z, I2);
        std::vector<std::vector<std::complex<double>>> ZX = tensor(Z, X);
        std::vector<std::vector<std::complex<double>>> ZY = tensor(Z, Y);
        std::vector<std::vector<std::complex<double>>> ZZ = tensor(Z, Z);
        std::vector<std::vector<std::complex<double>>> proj01 = tensor(proj0, proj1);
        std::vector<std::vector<std::complex<double>>> proj10 = tensor(proj1, proj0);
        std::vector<std::vector<std::complex<double>>> proj00 = tensor(proj0, proj0);
        std::vector<std::vector<std::complex<double>>> psi(d2, std::vector<std::complex<double>>(d2, 0));
        psi[0][0] = 1;
        psi[d2-1][d2-1] = 1;
        psi[0][d2-1] = 1;
        psi[d2-1][0] = 1;

        // Which to use
        std::vector<std::vector<std::vector<std::complex<double>>>> termsW = {

            // TODO ./causal -a -i 50
            // 37s -1.17 < -0.41
            // 52s -1.10 < -0.46
            // 48s -0.83 < -0.39
            // 1m -0.77 < -0.39
            //IIII,
            //ZZZI,
            //ZIXX,

            // 0.5694
            IIII,
            ZIZI,
            ZIII,
            IIZI,
            ZIIZ,
            IZZI,
            ZIZZ,
            ZZZI,
            ZIXX,
            ZIYY,
            XXZI,
            YYZI,

        };
        std::vector<std::vector<std::vector<std::complex<double>>>> termsA = {
            proj10,
            proj00,
            psi,
            //II,
            //IX,
            //IY,
            //IZ,
            //XI,
            //XX,
            //XY,
            //XZ,
            //YI,
            //YX,
            //YY,
            //YZ,
            //ZI,
            //ZX,
            //ZY,
            //ZZ,
        };
        std::vector<std::vector<std::vector<std::complex<double>>>> termsB = {
            proj10,
            proj00,
            psi,
            //II,
            //IX,
            //IY,
            //IZ,
            //XI,
            //XX,
            //XY,
            //XZ,
            //YI,
            //YX,
            //YY,
            //YZ,
            //ZI,
            //ZX,
            //ZY,
            //ZZ,
        };

        // W
        int nextVar = 0;
        for (int i=0; i<termsW.size(); i++) {
            for (int j=0; j<widthW; j++) {
                for (int k=j; k<widthW; k++) {
                    if (std::abs(termsW[i][j][k]) > 1e-8) {
                        W[j][k] += Polynomial<std::complex<double>>(numVars, 0.25, {nextVar}) * termsW[i][j][k];
                    }
                }
            }
            nextVar++;
        }
        for (int j=0; j<widthW; j++) {
            for (int k=j; k<widthW; k++) {
                W[j][k] = W[j][k].prune();
                W[k][j] = std::conj(W[j][k]);
            }
        }

        // The As
        for (int i=0; i<numAs; i++) {
            std::vector<std::vector<Polynomial<std::complex<double>>>> A(widthA, std::vector<Polynomial<std::complex<double>>>(widthA, Polynomial<std::complex<double>>(numVars, 0i)));
            for (int j=0; j<termsA.size(); j++) {
                for (int k=0; k<widthA; k++) {
                    for (int l=k; l<widthA; l++) {
                        if (std::abs(termsA[j][k][l]) > 1e-8) {
                            A[k][l] += Polynomial<std::complex<double>>(numVars, 0.5, {nextVar}) * termsA[j][k][l];
                        }
                    }
                }
                nextVar++;
            }
            for (int j=0; j<widthA; j++) {
                for (int k=j; k<widthA; k++) {
                    A[j][k] = A[j][k].prune();
                    A[k][j] = std::conj(A[j][k]);
                }
            }
            As.push_back(A);
        }

        // The Bs
        for (int i=0; i<numBs; i++) {
            std::vector<std::vector<Polynomial<std::complex<double>>>> B(widthB, std::vector<Polynomial<std::complex<double>>>(widthB, Polynomial<std::complex<double>>(numVars, 0i)));
            for (int j=0; j<termsB.size(); j++) {
                for (int k=0; k<widthB; k++) {
                    for (int l=k; l<widthB; l++) {
                        if (std::abs(termsB[j][k][l]) > 1e-8) {
                            B[k][l] += Polynomial<std::complex<double>>(numVars, 0.5, {nextVar}) * termsB[j][k][l];
                        }
                    }
                }
                nextVar++;
            }
            for (int j=0; j<widthB; j++) {
                for (int k=j; k<widthB; k++) {
                    B[j][k] = B[j][k].prune();
                    B[k][j] = std::conj(B[j][k]);
                }
            }
            Bs.push_back(B);
        }

    // Otherwise construct more generally
    } else {

        // Construct W
        if (verbosity > 0) {
            std::cout << "Constructing variable matrices..." << std::endl;
        }
        int nextVar = 0;
        for (int i=0; i<widthW; i++) {
            for (int j=i; j<widthW; j++) {
                if (i == j) {
                    W[i][j] = Polynomial<std::complex<double>>(numVars, 1, {nextVar});
                    nextVar++;
                } else {
                    W[i][j] = Polynomial<std::complex<double>>(numVars, 1, {nextVar}) + Polynomial<std::complex<double>>(numVars, 1i, {nextVar+1});
                    W[j][i] = std::conj(W[i][j]);
                    nextVar += 2;
                }
            }
        }

        // Construct As
        for (int i=0; i<numAs; i++) {
            std::vector<std::vector<Polynomial<std::complex<double>>>> A(widthA, std::vector<Polynomial<std::complex<double>>>(widthA, Polynomial<std::complex<double>>(numVars, 0i)));
            for (int j=0; j<widthA; j++) {
                for (int k=j; k<widthA; k++) {
                    if (j == k) {
                        A[j][k] = Polynomial<std::complex<double>>(numVars, 1, {nextVar});
                        nextVar++;
                    } else {
                        A[j][k] = Polynomial<std::complex<double>>(numVars, 1, {nextVar}) + Polynomial<std::complex<double>>(numVars, 1i, {nextVar+1});
                        nextVar+=2;
                    }
                    A[k][j] = std::conj(A[j][k]);
                }
            }
            As.push_back(A);
        }

        // Construct Bs
        for (int i=0; i<numBs; i++) {
            std::vector<std::vector<Polynomial<std::complex<double>>>> B(widthB, std::vector<Polynomial<std::complex<double>>>(widthB, Polynomial<std::complex<double>>(numVars, 0i)));
            for (int j=0; j<widthB; j++) {
                for (int k=j; k<widthB; k++) {
                    if (j == k) {
                        B[j][k] = Polynomial<std::complex<double>>(numVars, 1, {nextVar});
                        nextVar++;
                    } else {
                        B[j][k] = Polynomial<std::complex<double>>(numVars, 1, {nextVar}) + Polynomial<std::complex<double>>(numVars, 1i, {nextVar+1});
                        nextVar += 2;
                    }
                    B[k][j] = std::conj(B[j][k]);
                }
            }
            Bs.push_back(B);
        }

    }

    // If we should use known optimum for A and B
    if (fixedAB) {
        for (int i=0; i<numAs; i++) {
            for (int j=0; j<widthA; j++) {
                for (int k=0; k<widthA; k++) {
                    As[i][j][k] = knownAs[i][j][k];
                }
            }
        }
        for (int i=0; i<numBs; i++) {
            for (int j=0; j<widthB; j++) {
                for (int k=0; k<widthB; k++) {
                    Bs[i][j][k] = knownBs[i][j][k];
                }
            }
        }
    }

    // If each B should equal each A
    if (bEqualsA) {
        for (int i=0; i<numAs; i++) {
            for (int j=0; j<widthA; j++) {
                for (int k=0; k<widthA; k++) {
                    Bs[i][j][k] = As[i][j][k];
                }
            }
        }
    }

    // Print everything
    if (verbosity > 1) {
        std::cout << "W:" << std::endl;
        std::cout << W << std::endl;
        for (int i=0; i<numAs; i++) {
            std::cout << "A[" << i << "]:" << std::endl;
            std::cout << As[i] << std::endl;
        }
        for (int i=0; i<numBs; i++) {
            std::cout << "B[" << i << "]:" << std::endl;
            std::cout << Bs[i] << std::endl;
        }
    }

    // Positivity of everything
    if (verbosity >= 1) {
        std::cout << "Adding constraints..." << std::endl;
    }
    prob.conPSD.push_back(complexSDPToReal(W));
    for (int i=0; i<numAs; i++) {
        prob.conPSD.push_back(complexSDPToReal(As[i]));
    }
    for (int i=0; i<numBs; i++) {
        prob.conPSD.push_back(complexSDPToReal(Bs[i]));
    }

    // TODO add products of matrices
    if (useMatProds) {

        // all 2nd order hadamard products
        for (int i=0; i<numAs; i++) {
            for (int j=0; j<numBs; j++) {
                prob.conPSD.push_back(complexSDPToReal(hadamard(As[i], Bs[j])));
            }
        }
        for (int i=0; i<numAs; i++) {
            for (int j=i; j<numAs; j++) {
                prob.conPSD.push_back(complexSDPToReal(hadamard(As[i], As[j])));
            }
        }
        for (int i=0; i<numBs; i++) {
            for (int j=i; j<numBs; j++) {
                prob.conPSD.push_back(complexSDPToReal(hadamard(Bs[i], Bs[j])));
            }
        }
        for (int i=0; i<numAs; i++) {
            prob.conPSD.push_back(complexSDPToReal(hadamard(As[i], W)));
        }
        for (int i=0; i<numBs; i++) {
            prob.conPSD.push_back(complexSDPToReal(hadamard(Bs[i], W)));
        }

        // all 3rd order hadamard products
        for (int i=0; i<numAs; i++) {
            for (int j=0; j<numBs; j++) {
                prob.conPSD.push_back(complexSDPToReal(hadamard(As[i], hadamard(Bs[j], W))));
            }
        }

        // since A and B commute, AB is positive
        for (int i=0; i<numAs; i++) {
            for (int j=0; j<numBs; j++) {
                prob.conPSD.push_back(complexSDPToReal(As[i] * Bs[j]));
            }
        }

    }

    // Identities
    std::vector<std::vector<Polynomial<std::complex<double>>>> I_2(d, std::vector<Polynomial<std::complex<double>>>(d, Polynomial<std::complex<double>>(numVars, 0i)));
    for (int i=0; i<d; i++) {
        I_2[i][i] = Polynomial<std::complex<double>>(numVars, 1.0/d);
    }
    std::vector<std::vector<Polynomial<std::complex<double>>>> I_4(d2, std::vector<Polynomial<std::complex<double>>>(d2, Polynomial<std::complex<double>>(numVars, 0i)));
    for (int i=0; i<d2; i++) {
        I_4[i][i] = Polynomial<std::complex<double>>(numVars, 1.0/d2);
    }
    std::vector<std::vector<Polynomial<std::complex<double>>>> I_8(d3, std::vector<Polynomial<std::complex<double>>>(d3, Polynomial<std::complex<double>>(numVars, 0i)));
    for (int i=0; i<d3; i++) {
        I_8[i][i] = Polynomial<std::complex<double>>(numVars, 1.0/d3);
    }

    // Constraints on the W
    // The order of W is AI, AO, BI, BO
    // W_something means just the W_{something} part (the rest traced out and replaced with identity)
    Polynomial<std::complex<double>> trace_W = trace(W);
    std::vector<std::vector<Polynomial<std::complex<double>>>> W_AI_AO = tensor(tensor(partialTrace(partialTrace(W, 4, d, 3), 3, d, 2), I_2), I_2);
    std::vector<std::vector<Polynomial<std::complex<double>>>> W_BI_BO = tensor(I_2, tensor(I_2, partialTrace(partialTrace(W, 4, d, 0), 3, d, 0)));
    std::vector<std::vector<Polynomial<std::complex<double>>>> W_AI = tensor(partialTrace(partialTrace(partialTrace(W, 4, d, 3), 3, d, 2), 2, d, 1), I_8);
    std::vector<std::vector<Polynomial<std::complex<double>>>> W_BI = tensor(I_4, tensor(partialTrace(partialTrace(partialTrace(W, 4, d, 3), 3, d, 0), 2, d, 0), I_2));
    std::vector<std::vector<Polynomial<std::complex<double>>>> W_AI_AO_BI = tensor(partialTrace(W, 4, d, 3), I_2);
    std::vector<std::vector<Polynomial<std::complex<double>>>> W_AI_BI_BO = reorderSubsystems(tensor(I_2, partialTrace(W, 4, d, 1)), 4, d, {1, 0, 2, 3});
    std::vector<std::vector<Polynomial<std::complex<double>>>> W_AI_BI = reorderSubsystems(tensor(I_4, partialTrace(partialTrace(W, 4, d, 3), 3, d, 1)), 4, d, {2, 1, 3, 0});

    // tr_W  = d^2 
    Polynomial<std::complex<double>> traceCon = trace_W - Polynomial<std::complex<double>>(numVars, d*d);
    prob.conZero.push_back(real(traceCon));
    prob.conZero.push_back(imag(traceCon));

    // W_{AI,AO} = W_{AI}
    for (int i=0; i<widthW; i++) {
        for (int j=i; j<widthW; j++) {
            Polynomial<std::complex<double>> con = W_AI_AO[i][j] - W_AI[i][j];
            prob.conZero.push_back(real(con));
            prob.conZero.push_back(imag(con));
        }
    }

    // W_{BI,BO} = W_{BI}
    for (int i=0; i<widthW; i++) {
        for (int j=i; j<widthW; j++) {
            Polynomial<std::complex<double>> con = W_BI_BO[i][j] - W_BI[i][j];
            prob.conZero.push_back(real(con));
            prob.conZero.push_back(imag(con));
        }
    }

    // W = W_{AI,AO,BI} + W_{AI,BI,BO} - W_{AI,BI}
    for (int i=0; i<widthW; i++) {
        for (int j=i; j<widthW; j++) {
            Polynomial<std::complex<double>> con = W_AI_AO_BI[i][j] + W_AI_BI_BO[i][j] - W_AI_BI[i][j] - W[i][j];
            prob.conZero.push_back(real(con));
            prob.conZero.push_back(imag(con));
        }
    }

    // Constraints on the As
    // tr_AO \sum_a A^{AO,AI}_{a,x} = I
    std::vector<std::vector<Polynomial<std::complex<double>>>> A_0_0_AI = partialTrace(As[0], 2, d, 1);
    std::vector<std::vector<Polynomial<std::complex<double>>>> A_0_1_AI = partialTrace(As[1], 2, d, 1);
    std::vector<std::vector<Polynomial<std::complex<double>>>> A_1_0_AI = partialTrace(As[2], 2, d, 1);
    std::vector<std::vector<Polynomial<std::complex<double>>>> A_1_1_AI = partialTrace(As[3], 2, d, 1);
    std::vector<std::vector<Polynomial<std::complex<double>>>> A_0_sum = A_0_0_AI + A_1_0_AI;
    std::vector<std::vector<Polynomial<std::complex<double>>>> A_1_sum = A_0_1_AI + A_1_1_AI;
    for (int i=0; i<A_0_sum.size(); i++) {
        for (int j=i; j<A_0_sum.size(); j++) {
            if (i == j) {
                Polynomial<std::complex<double>> con1 = A_0_sum[i][j] - Polynomial<std::complex<double>>(numVars, 1);
                Polynomial<std::complex<double>> con2 = A_1_sum[i][j] - Polynomial<std::complex<double>>(numVars, 1);
                prob.conZero.push_back(real(con1));
                prob.conZero.push_back(imag(con1));
                prob.conZero.push_back(real(con2));
                prob.conZero.push_back(imag(con2));
            } else {
                prob.conZero.push_back(real(A_0_sum[i][j]));
                prob.conZero.push_back(imag(A_0_sum[i][j]));
                prob.conZero.push_back(real(A_1_sum[i][j]));
                prob.conZero.push_back(imag(A_1_sum[i][j]));
            }
        }
    }

    // Constraints on the Bs
    // tr_BO \sum_b B^{BO,BI}_{b,x} = I
    std::vector<std::vector<Polynomial<std::complex<double>>>> B_0_0_BI = partialTrace(Bs[0], 2, d, 1);
    std::vector<std::vector<Polynomial<std::complex<double>>>> B_0_1_BI = partialTrace(Bs[1], 2, d, 1);
    std::vector<std::vector<Polynomial<std::complex<double>>>> B_1_0_BI = partialTrace(Bs[2], 2, d, 1);
    std::vector<std::vector<Polynomial<std::complex<double>>>> B_1_1_BI = partialTrace(Bs[3], 2, d, 1);
    std::vector<std::vector<Polynomial<std::complex<double>>>> B_0_sum = B_0_0_BI + B_1_0_BI;
    std::vector<std::vector<Polynomial<std::complex<double>>>> B_1_sum = B_0_1_BI + B_1_1_BI;
    for (int i=0; i<B_0_sum.size(); i++) {
        for (int j=i; j<B_0_sum.size(); j++) {
            if (i == j) {
                Polynomial<std::complex<double>> con1 = B_0_sum[i][j] - Polynomial<std::complex<double>>(numVars, 1);
                Polynomial<std::complex<double>> con2 = B_1_sum[i][j] - Polynomial<std::complex<double>>(numVars, 1);
                prob.conZero.push_back(real(con1));
                prob.conZero.push_back(imag(con1));
                prob.conZero.push_back(real(con2));
                prob.conZero.push_back(imag(con2));
            } else {
                prob.conZero.push_back(real(B_0_sum[i][j]));
                prob.conZero.push_back(imag(B_0_sum[i][j]));
                prob.conZero.push_back(real(B_1_sum[i][j]));
                prob.conZero.push_back(imag(B_1_sum[i][j]));
            }
        }
    }

    // Objective
    if (verbosity >= 1) {
        std::cout << "Constructing objective..." << std::endl;
    }
    prob.obj = Polynomial<double>(numVars);

    // GYNI = p(a = y, b = x)
    //      = 0.25*(p_0_0_0_0 + p_0_1_1_0 + p_1_0_0_1 + p_1_1_1_1);
    if (problemName == "gyni") {
        for (int y=0; y<2; y++) {
            for (int x=0; x<2; x++) {
                int indA = x*2 + y;
                int indB = y*2 + x;
                Polynomial<std::complex<double>> p = trace(tensor(As[indA], Bs[indB]) * W).prune();
                if (verbosity >= 3) {
                    std::cout << "Prob " << x << " " << y << " : " << p << std::endl;
                }
                prob.obj += real(p);
            }
        }
        prob.obj = prob.obj * 0.25;
    } else if (problemName == "lgyni") {
        for (int x=0; x<2; x++) {
            for (int y=0; y<2; y++) {
                for (int a=0; a<2; a++) {
                    for (int b=0; b<2; b++) {
                        if (x*(a ^ y) == 0 && y*(b ^ x) == 0) {
                            int indA = a*2 + x;
                            int indB = b*2 + y;
                            Polynomial<std::complex<double>> p = trace(tensor(As[indA], Bs[indB]) * W).prune();
                            if (verbosity >= 3) {
                                std::cout << "Prob " << x << " " << y << " " << a << " " << b << " : " << p << std::endl;
                            }
                            prob.obj += real(p);
                        }
                    }
                }
            }
        }
        prob.obj = prob.obj * 0.25;
    } else if (problemName == "rlgyni") { // TODO
        for (int x=0; x<2; x++) {
            for (int y=0; y<2; y++) {
                for (int a=0; a<2; a++) {
                    for (int b=0; b<2; b++) {
                        if (x*(a ^ y) == 0 && y*(b ^ x) == 0) {
                            int indA = a*2 + x;
                            int indB = b*2 + y;
                            Polynomial<std::complex<double>> p = trace(tensor(As[indA], Bs[indB]) * W).prune();
                            if (verbosity >= 3) {
                                std::cout << "Prob " << x << " " << y << " " << a << " " << b << " : " << p << std::endl;
                            }
                            double randCoeff = ((double) rand() / (RAND_MAX)) * 2 - 1;
                            prob.obj += randCoeff * real(p);
                            ineqString += std::to_string(-randCoeff) + "<A" + std::to_string(indA) + "B" + std::to_string(indB) + "> +";
                        }
                    }
                }
            }
        }
    }

    // Tidy the inequality string
    // remove all spaces
    std::string newIneqString;
    for (size_t i=0; i<ineqString.size(); i++) {
        if (ineqString[i] != ' ') {
            newIneqString += ineqString[i];
        }
    }
    ineqString = newIneqString;
    // replace "+-" with "-"
    size_t pos = 0;
    while ((pos = ineqString.find("+-", pos)) != std::string::npos) {
        ineqString.replace(pos, 2, "-");
        pos += 1;
    }

    // Custom moment matrix
    //std::vector<std::vector<std::vector<Polynomial<std::complex<double>>>>> matsInTopRow;
    //std::vector<std::vector<Polynomial<std::complex<double>>>> fullI(widthW, std::vector<Polynomial<std::complex<double>>>(widthW, Polynomial<std::complex<double>>(numVars, 0i)));
    //std::vector<std::vector<Polynomial<std::complex<double>>>> partialI(widthA, std::vector<Polynomial<std::complex<double>>>(widthA, Polynomial<std::complex<double>>(numVars, 0i)));
    //for (int i=0; i<widthW; i++) {
        //fullI[i][i] = Polynomial<std::complex<double>>(numVars, 1);
    //}
    //for (int i=0; i<widthA; i++) {
        //partialI[i][i] = Polynomial<std::complex<double>>(numVars, 1);
    //}
    //matsInTopRow.push_back(fullI);
    //matsInTopRow.push_back(W);
    //std::cout << "W done" << std::endl;
    //for (int i=0; i<numAs; i++) {
        //matsInTopRow.push_back(tensor(As[i], partialI));
        //std::cout << "A[" << i << "] done" << std::endl;
    //}
    //for (int i=0; i<numBs; i++) {
        //matsInTopRow.push_back(tensor(partialI, Bs[i]));
        //std::cout << "B[" << i << "] done" << std::endl;
    //}
    //for (int i=0; i<numBs; i++) {
        //matsInTopRow.push_back(tensor(partialI, Bs[i]) * W);
        //std::cout << "B[" << i << "]W done" << std::endl;
    //}

    // Form the custom matrix based on the products of the above matrices
    //int customMatWidth = matsInTopRow.size() * widthW;
    //std::vector<std::vector<Polynomial<double>>> customMomentMat(customMatWidth, std::vector<Polynomial<double>>(customMatWidth, Polynomial<double>(numVars, 0)));
    //for (int i=0; i<matsInTopRow.size(); i++) {
        //for (int j=0; j<matsInTopRow.size(); j++) {
            //std::cout << "Multiplying " << i << " " << j << std::endl;
            //std::vector<std::vector<Polynomial<std::complex<double>>>> mat = matsInTopRow[i] * matsInTopRow[j];
            //for (int k=0; k<mat.size(); k++) {
                //for (int l=0; l<mat.size(); l++) {
                    //customMomentMat[widthW*i+k][widthW*j+l] += real(mat[k][l]);
                //}
            //}
        //}
    //}
    //level = 0;
    //prob.conPSD.push_back(customMomentMat);

    // Try restricting to a subset of the objective
    //Polynomial<double> objSubset = Polynomial<double>(numVars, 0);
    //int numTerms = 5;
    //for (int i=0; i<numTerms; i++) {
        //int randInd = rand() % prob.obj.size();
        //int count = 0;
        //for (auto it=prob.obj.coeffs.begin(); it!=prob.obj.coeffs.end(); it++) {
            //if (count == randInd) {
                //objSubset.coeffs[it->first] = it->second;
                //break;
            //}
            //count++;
        //}
    //}
    //prob.obj = objSubset;

    // Clean all the constraints
    if (verbosity > 0) {
        std::cout << "Cleaning constraints..." << std::endl;
    }
    std::vector<Polynomial<double>> conZero;
    for (int i=0; i<prob.conZero.size(); i++) {
        Polynomial<double> p = prob.conZero[i].prune();
        if (p.size() > 0) {
            conZero.push_back(p);
        }
    }
    prob.conZero = conZero;

    // If verbose, output the problem
    if (verbosity >= 2) {
        std::cout << "Problem:" << std::endl;
        std::cout << prob << std::endl;
        std::cout << "Num linear cons: " << prob.conZero.size() << std::endl;
        std::cout << "Size of obj: " << prob.obj.size() << std::endl;
        std::cout << "Size of moment mats: ";
        for (int i=0; i<prob.conPSD.size(); i++) {
            std::cout << prob.conPSD[i].size() << " x " << prob.conPSD[i][0].size() << ", ";
        }
        std::cout << std::endl;
    }

    // Try reducing the problem
    if (verbosity > 0) {
        std::cout << "Reducing problem..." << std::endl;
    }
    prob = prob.removeTrivial();
    std::vector<Polynomial<double>> conZeroNew;
    for (int i=0; i<prob.conZero.size(); i++) {
        Polynomial<double> p = prob.conZero[i].prune();
        if (p.size() > 0) {
            conZeroNew.push_back(p);
        }
    }
    prob.conZero = conZeroNew;

    // If verbose, output the problem
    if (verbosity > 1) {
        std::cout << "Problem after reduction:" << std::endl;
        std::cout << prob << std::endl;
        std::cout << "Num linear cons: " << prob.conZero.size() << std::endl;
        std::cout << "Size of obj: " << prob.obj.size() << std::endl;
        std::cout << "Size of moment mat: " << prob.conPSD.size() << " x " << prob.conPSD[0].size() << std::endl;
    }

    // Put everything in
    std::vector<double> valueVec;
    if (useKnown) {
        valueVec = std::vector<double>(numVars, 0);
        std::unordered_map<int, double> valueMap;
        int nextVarInd = 0;
        for (int i=0; i<widthW; i++) {
            for (int j=i; j<widthW; j++) {
                if (i == j) {
                    valueMap[nextVarInd] = std::real(knownW[i][j]);
                    nextVarInd++;
                } else {
                    valueMap[nextVarInd] = std::real(knownW[i][j]);
                    valueMap[nextVarInd+1] = std::imag(knownW[i][j]);
                    nextVarInd+=2;
                }
            }
        }
        for (int i=0; i<numAs; i++) {
            for (int j=0; j<widthA; j++) {
                for (int k=j; k<widthA; k++) {
                    if (j == k) {
                        valueMap[nextVarInd] = std::real(knownAs[i][j][k]);
                        nextVarInd++;
                    } else {
                        valueMap[nextVarInd] = std::real(knownAs[i][j][k]);
                        valueMap[nextVarInd+1] = std::imag(knownAs[i][j][k]);
                        nextVarInd+=2;
                    }
                }
            }
        }
        for (int i=0; i<numBs; i++) {
            for (int j=0; j<widthB; j++) {
                for (int k=j; k<widthB; k++) {
                    if (j == k) {
                        valueMap[nextVarInd] = std::real(knownBs[i][j][k]);
                        nextVarInd++;
                    } else {
                        valueMap[nextVarInd] = std::real(knownBs[i][j][k]);
                        valueMap[nextVarInd+1] = std::imag(knownBs[i][j][k]);
                        nextVarInd+=2;
                    }
                }
            }
        }
        for (int i=0; i<numVars; i++) {
            valueVec[i] = valueMap[i];
        }
        std::cout << "Checking known solution..." << std::endl;
        bool feasible = prob.isFeasible(valueVec);
        std::cout << "Feasible: " << feasible << std::endl;
        std::cout << "Value: " << prob.obj.eval(valueVec) << std::endl;
        if (!justPrintHessian && !seesaw) {
            return 0;
        }
    }

    // Switch the minimal variable mapping
    if (useAnsatz) {
        std::cout << "Switching to minimal variable mapping..." << std::endl;
        std::unordered_map<int, int> minimalMap = prob.getMinimalMap();
        prob = prob.replaceWithVariable(minimalMap);
        std::cout << "Num variables after reduction: " << prob.maxVariables << std::endl;
    }

    // Add bounds and make it a maximization problem
    if (useAnsatz) {
        for (int i=0; i<prob.maxVariables; i++) {
            prob.varBounds[i] = std::make_pair(-1, 1);
        }
    } else {
        for (int i=0; i<prob.maxVariables; i++) {
            prob.varBounds[i] = std::make_pair(-d, d);
        }
    }
    prob.obj = -prob.obj;

    // Generate two random points and plot the objective between them
    //std::vector<std::vector<double>> points = prob.manyFeasible(10, verbosity);
    //points.push_back(valueVec);
    //int steps = 100;
    //std::ofstream outfile("temp2");

    // 2D
    //for (int i=0; i<steps; i++) {

        //// The mix
        //double t1 = double(i)/double(steps);

        //// For each combination of points
        //std::vector<double> results;
        //for (int j1=0; j1<points.size(); j1++) {
            //for (int j2=j1+1; j2<points.size(); j2++) {

                //// Mixing the two
                //std::vector<double> point(prob.maxVariables, 0);
                //for (int l=0; l<prob.maxVariables; l++) {
                    //point[l] = (1-t1)*points[j1][l] + t1*points[j2][l];
                //}

                //// Evaluate the objective
                //double obj = prob.obj.eval(point);
                //results.push_back(obj);

            //}
        //}

        //// Write it to file
        //outfile << t1 << ", ";
        //for (int j=0; j<results.size(); j++) {
            //outfile << results[j];
            //if (j < results.size()-1) {
                //outfile << ", ";
            //}
        //}
        //outfile << std::endl;

    //}

    // 3D
    //for (int i=0; i<steps; i++) {
        //for (int k=0; k<steps; k++) {

            //// The mix
            //double t1 = double(i)/double(steps);
            //double t2 = double(k)/double(steps);

            //// For each combination of points
            //std::vector<double> results;
            //for (int j1=0; j1<points.size(); j1++) {
                //for (int j2=j1+1; j2<points.size(); j2++) {
                    //for (int j3=j2+1; j3<points.size(); j3++) {

                        //// Mixing all three
                        //std::vector<double> point(prob.maxVariables, 0);
                        //for (int l=0; l<prob.maxVariables; l++) {
                            //point[l] = t1*points[j1][l] + t2*points[j2][l] + (1-t1-t2)*points[j3][l];
                        //}

                        //// Evaluate the objective
                        //double obj = prob.obj.eval(point);
                        //results.push_back(obj);

                    //}
                //}
            //}

            //// Write it to file
            //outfile << t1 << ", ";
            //outfile << t2 << ", ";
            //for (int j=0; j<results.size(); j++) {
                //outfile << results[j];
                //if (j < results.size()-1) {
                    //outfile << ", ";
                //}
            //}
            //outfile << std::endl;

        //}
    //}

    //return 0;

    // See how convex it is using TSNE, writing to a file
    //std::ofstream outfile("temp2");
    //std::ofstream outfile2("temp3");
    //auto points = prob.manyRandom(100, 1000);
    //for (int i=0; i<points.size(); i++) {
        //for (int j=0; j<points[i].size(); j++) {
            //outfile << points[i][j];
            //if (j < points[i].size()-1) {
                //outfile << "\t";
            //}
        //}
        //outfile << std::endl;
        //outfile2 << prob.obj.eval(points[i]) << std::endl;
    //}
    //return 0;

    // Get approximate linear objective
    //std::vector<double> direction(prob.maxVariables);
    //direction = valueVec;
    //for (int j=0; j<prob.maxVariables; j++) {

        //// Looking for this var
        //int digitsPerInd = prob.obj.digitsPerInd;
        //std::string indToFind = std::to_string(j);
        //indToFind.insert(0, digitsPerInd-indToFind.size(), ' ');

        //// Get the sum of all terms that contain this variable
        //double coeffSum = 0;
        //for (auto const &pair: prob.obj.coeffs) {
            
            //// Check if this term contains our variable
            //for (int j=0; j<pair.first.size(); j+=digitsPerInd) {
                //if (pair.first.substr(j, digitsPerInd) == indToFind) {
                    //coeffSum += pair.second / double(pair.first.size()/digitsPerInd);
                    //break;
                //}

            //}

        //}

        //// For now the direction is just the sum
        //direction[j] = -coeffSum;

    //}

    // This is the linear objective
    //Polynomial<double> linearObj = Polynomial<double>(numVars, 0);
    //for (int j=0; j<prob.maxVariables; j++) {
        //linearObj += (-direction[j])*Polynomial<double>(numVars, 1, {j});
    //}
    //Polynomial<double> trueObj = prob.obj;
    //prob.obj = linearObj;

    // If told to just print the Hessian
    if (justPrintHessian) {

        // Mess around with the objective
        std::cout << "Modifying the objective..." << std::endl;
        prob.obj = prob.obj * prob.obj;

        // Output the terminal
        std::cout << "Calculating Hessian..." << std::endl;
        std::vector<std::vector<Polynomial<double>>> hessian = prob.obj.hessian(8);

        // If using known, check the Hessian with the known optimum
        if (useKnown) {

            // Evaluate the Hessian at the known optimum
            std::cout << "Checking Hessian with known optimum..." << std::endl;
            Eigen::MatrixXd hessianMat(hessian.size(), hessian[0].size());
            for (int i=0; i<hessian.size(); i++) {
                for (int j=0; j<hessian[i].size(); j++) {
                    hessianMat(i, j) = hessian[i][j].eval(valueVec);
                }
            }

            // Get the eigenvalues
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(hessianMat);
            Eigen::VectorXd eigenvalues = es.eigenvalues();
            double minEigenvalue = eigenvalues.minCoeff();
            double maxEigenvalue = eigenvalues.maxCoeff();
            std::cout << "Eigenvalues from " << minEigenvalue << " to " << maxEigenvalue << std::endl;

        }

        // Write it as a csv
        std::ofstream outfile("hessian.csv");
        for (int i=0; i<hessian.size(); i++) {
            for (int j=0; j<hessian[i].size(); j++) {
                if (hessian[i][j].size() > 0) {
                    outfile << 1;
                } else {
                    outfile << 0;
                }
                if (j < hessian[i].size()-1) {
                    outfile << ", ";
                }
            }
            outfile << std::endl;
        }
        outfile.close();

        return 0;

    }

    // If seesawing
    std::pair<double, std::vector<double>> res;
    if (seesaw) {

        // Determine the variable sets
        std::set<int> varSetA;
        std::set<int> varSetB;
        std::set<int> varSetW;
        for (int i=0; i<numAs; i++) {
            for (int j=0; j<As[i].size(); j++) {
                for (int k=j; k<As[i][j].size(); k++) {
                    std::vector<int> vars = As[i][j][k].getVariables();
                    for (auto it=vars.begin(); it!=vars.end(); it++) {
                        varSetA.insert(*it);
                    }
                }
            }
        }
        for (int i=0; i<numBs; i++) {
            for (int j=0; j<Bs[i].size(); j++) {
                for (int k=j; k<Bs[i][j].size(); k++) {
                    std::vector<int> vars = Bs[i][j][k].getVariables();
                    for (auto it=vars.begin(); it!=vars.end(); it++) {
                        varSetB.insert(*it);
                    }
                }
            }
        }
        for (int i=0; i<W.size(); i++) {
            for (int j=i; j<W[i].size(); j++) {
                std::vector<int> vars = W[i][j].getVariables();
                for (auto it=vars.begin(); it!=vars.end(); it++) {
                    varSetW.insert(*it);
                }
            }
        }

        // Optimize using seesaw
        if (verbosity >= 1) {
            std::cout << "Optimizing using seesaw..." << std::endl;
        }
        res = prob.seesaw(verbosity, maxIters, {varSetW, varSetA, varSetB}, valueVec);
        if (verbosity >= 1) {
            std::cout << "Result: " << -res.first << std::endl;
        }

    // If branch and bounding
    } else {

        // Optimize using branch and bound
        if (verbosity >= 1) {
            std::cout << "Optimizing..." << std::endl;
        }
        res = prob.optimize(level, verbosity, maxIters);
        if (verbosity >= 1) {
            std::cout << "Result: " << -res.first << std::endl;
        }
        bool isFeas = prob.isFeasible(res.second);
        if (verbosity >= 1) {
            std::cout << "Feasible: " << isFeas << std::endl;
        }
        
    }

    // Output the optimal W
    //int nextVarInd = 0;
    //Eigen::MatrixXcd Wopt = Eigen::MatrixXcd::Zero(widthW, widthW);
    //for (int i=0; i<widthW; i++) {
        //for (int j=i; j<widthW; j++) {
            //if (i == j) {
                //Wopt(i, j) = res2.second[nextVarInd];
                //nextVarInd++;
            //} else {
                //Wopt(i, j) = res2.second[nextVarInd] + 1i*res2.second[nextVarInd+1];
                //Wopt(j, i) = std::conj(Wopt(i, j));
                //nextVarInd+=2;
            //}
        //}
    //}
    //// set near 0 elements to 0
    //for (int i=0; i<widthW; i++) {
        //for (int j=0; j<widthW; j++) {
            //if (std::abs(Wopt(i, j)) < 1e-6) {
                //Wopt(i, j) = 0;
            //}
            //if (std::abs(std::imag(Wopt(i, j))) < 1e-6) {
                //Wopt(i, j) = std::real(Wopt(i, j));
            //}
            //if (std::abs(std::real(Wopt(i, j))) < 1e-6) {
                //Wopt(i, j) = 1i*std::imag(Wopt(i, j));
            //}
        //}
    //}
    //std::cout << "Wopt:" << std::endl;
    //std::cout << Wopt << std::endl;

    // Output the ineqString if it is not empty TODO
    if (outputFinal) {
        std::cout << -res.first << ineqString << std::endl;
    }

	return 0;

}
