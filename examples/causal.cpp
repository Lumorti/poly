#include "../poly.h"

#include <vector>
#include <stdexcept>
#include <iostream>

// Generic addition of two vector of vectors
template <typename T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
    std::vector<std::vector<T>> C(A.size(), std::vector<T>(A[0].size(), T(0)));
    for (int i=0; i<A.size(); i++) {
        for (int j=0; j<A[0].size(); j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

// Trace of a matrix
template <typename T>
T trace(const std::vector<std::vector<T>>& A) {
    T tr = T(0);
    for (int i=0; i<A.size(); i++) {
        tr += A[i][i];
    }
    return tr;
}

// Matrix product of two matrices
template <typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
    std::vector<std::vector<T>> C(A.size(), std::vector<T>(B[0].size(), T(0)));
    for (int i=0; i<A.size(); i++) {
        for (int j=0; j<B[0].size(); j++) {
            for (int k=0; k<A[0].size(); k++) {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return C;
}

// Tensor product of two matrices
template <typename T>
std::vector<std::vector<T>> tensor(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
    std::vector<std::vector<T>> C(A.size()*B.size(), std::vector<T>(A[0].size()*B[0].size(), T(0)));
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
std::vector<std::vector<T>> partialTrace(
    const std::vector<std::vector<T>>& mat,
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
    std::vector<std::vector<T>> result(dimPartial, std::vector<T>(dimPartial, T(0)));

    // For each row/col in the d^(n-1) x d^(n-1) partial trace,
    // we sum over x in [0, d-1], i.e. the index that was "traced out".
    for(int rowP = 0; rowP < dimPartial; ++rowP) {
        for(int colP = 0; colP < dimPartial; ++colP) {

            T sumVal = T(0);

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
 *   // Original system: A (x) B (x) C => p = {0,1,2} (identity)
 *   // Suppose we want B (x) A (x) C => p = {1,0,2}
 *   // We reorder a 2^3 x 2^3 matrix from ABC-basis to BAC-basis.
 */
template <typename T>
std::vector<std::vector<T>> reorderSubsystems(
    const std::vector<std::vector<T>>& mat,
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
    std::vector<std::vector<T>> newMat(dimFull, std::vector<T>(dimFull, T(0)));

    // 4) For each row/col in the old matrix, determine where it goes in the new matrix
    for(int oldRow = 0; oldRow < dimFull; ++oldRow) {
        // decode the base-d digits for oldRow
        std::vector<int> oldRowDigits = indexToMultiIndex(oldRow, n, d);

        // reorder the digits: newRowDigits[i] = oldRowDigits[p[i]]
        std::vector<int> newRowDigits(n);
        for(int i = 0; i < n; ++i) {
            newRowDigits[i] = oldRowDigits[p[i]];
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

	// Generate the problem for GYNI
    int widthW = std::pow(d, 4);
    int numAs = 4;
    int numBs = 4;
    int widthA = std::pow(d, 2);
    int widthB = std::pow(d, 2);
    int fullWidth = widthW + numAs*widthA + numBs*widthB;
    int numVars = widthW*(widthW+1)/2 + numAs*widthA*(widthA+1)/2 + numBs*widthB*(widthB+1)/2;
	PolynomialProblem<double> prob(numVars);

    // Output some sizes
    std::cout << "Full matrix size: " << fullWidth << std::endl;
    std::cout << "Number of variables: " << numVars << std::endl;
    std::cout << "Width of W: " << widthW << std::endl;
    std::cout << "Width of A: " << widthA << std::endl;
    std::cout << "Width of B: " << widthB << std::endl;

    // Construct W
    std::vector<std::vector<Polynomial<double>>> W(widthW, std::vector<Polynomial<double>>(widthW, Polynomial<double>(numVars, 0)));
    int nextVar = 0;
    for (int i=0; i<widthW; i++) {
        for (int j=i; j<widthW; j++) {
            W[i][j] = Polynomial<double>(numVars, 1, {nextVar});
            W[j][i] = W[i][j];
            nextVar++;
        }
    }

    // Construct As
    std::vector<std::vector<std::vector<Polynomial<double>>>> As;
    for (int i=0; i<numAs; i++) {
        std::vector<std::vector<Polynomial<double>>> A(widthA, std::vector<Polynomial<double>>(widthA, Polynomial<double>(numVars, 0)));
        for (int j=0; j<widthA; j++) {
            for (int k=j; k<widthA; k++) {
                A[j][k] = Polynomial<double>(numVars, 1, {nextVar});
                A[k][j] = A[j][k];
                nextVar++;
            }
        }
        As.push_back(A);
    }

    // Construct Bs
    std::vector<std::vector<std::vector<Polynomial<double>>>> Bs;
    for (int i=0; i<numBs; i++) {
        std::vector<std::vector<Polynomial<double>>> B(widthB, std::vector<Polynomial<double>>(widthB, Polynomial<double>(numVars, 0)));
        for (int j=0; j<widthB; j++) {
            for (int k=j; k<widthB; k++) {
                B[j][k] = Polynomial<double>(numVars, 1, {nextVar});
                B[k][j] = B[j][k];
                nextVar++;
            }
        }
        Bs.push_back(B);
    }

    // Positivity of everything
    std::vector<std::vector<Polynomial<double>>> psdCon(fullWidth, std::vector<Polynomial<double>>(fullWidth, Polynomial<double>(numVars, 0)));
    int delta = 0;
    for (int i=0; i<widthW; i++) {
        for (int j=i; j<widthW; j++) {
            psdCon[i+delta][j+delta] = W[i][j];
            psdCon[j+delta][i+delta] = W[j][i];
        }
    }
    delta += widthW;
    for (int i=0; i<numAs; i++) {
        for (int j=0; j<widthA; j++) {
            for (int k=j; k<widthA; k++) {
                psdCon[j+delta][k+delta] = As[i][j][k];
                psdCon[k+delta][j+delta] = As[i][k][j];
            }
        }
        delta += widthA;
    }
    for (int i=0; i<numBs; i++) {
        for (int j=0; j<widthB; j++) {
            for (int k=j; k<widthB; k++) {
                psdCon[j+delta][k+delta] = As[i][j][k];
                psdCon[k+delta][j+delta] = As[i][k][j];
            }
        }
        delta += widthB;
    }
    prob.conPSD = psdCon;

    // Identities
    std::vector<std::vector<Polynomial<double>>> I_2(2, std::vector<Polynomial<double>>(2, Polynomial<double>(numVars, 0)));
    for (int i=0; i<2; i++) {
        I_2[i][i] = Polynomial<double>(numVars, 1);
    }
    std::vector<std::vector<Polynomial<double>>> I_4 = tensor(I_2, I_2);
    std::vector<std::vector<Polynomial<double>>> I_8 = tensor(I_4, I_2);

    // Constraints on the W
    // The order of W is AI, AO, BI, BO
    // W_something is tracing out that something then tensoring with the identity
    Polynomial<double> trace_W = trace(W);
    std::vector<std::vector<Polynomial<double>>> W_AI_AO = tensor(tensor(partialTrace(partialTrace(W, 4, d, 3), 3, d, 2), I_2), I_2);
    std::vector<std::vector<Polynomial<double>>> W_BI_BO = tensor(I_2, tensor(I_2, partialTrace(partialTrace(W, 4, d, 0), 3, d, 0)));
    std::vector<std::vector<Polynomial<double>>> W_AI = tensor(partialTrace(partialTrace(partialTrace(W, 4, d, 3), 3, d, 2), 2, d, 1), I_8);
    std::vector<std::vector<Polynomial<double>>> W_BI = tensor(I_4, tensor(partialTrace(partialTrace(partialTrace(W, 4, d, 3), 3, d, 0), 2, d, 0), I_2));
    std::vector<std::vector<Polynomial<double>>> W_AI_AO_BI = tensor(partialTrace(W, 4, d, 3), I_2);
    std::vector<std::vector<Polynomial<double>>> W_AI_BI_BO = reorderSubsystems(tensor(I_2, partialTrace(W, 4, d, 1)), 4, d, {1, 0, 2, 3});
    std::vector<std::vector<Polynomial<double>>> W_AI_BI = reorderSubsystems(tensor(I_4, partialTrace(partialTrace(W, 4, d, 3), 3, d, 1)), 4, d, {2, 0, 3, 1});
    // tr_W  = d^2 
    prob.conZero.push_back(trace_W - Polynomial<double>(numVars, 1, {d*d}));
    // W_{AI,AO} = W_{AI}
    for (int i=0; i<widthW; i++) {
        for (int j=i; j<widthW; j++) {
            prob.conZero.push_back(W_AI_AO[i][j] - W_AI[i][j]);
            prob.conZero.push_back(W_AI_AO[j][i] - W_AI[j][i]);
        }
    }
    // W_{BI,BO} = W_{BI}
    for (int i=0; i<widthW; i++) {
        for (int j=i; j<widthW; j++) {
            prob.conZero.push_back(W_BI_BO[i][j] - W_BI[i][j]);
            prob.conZero.push_back(W_BI_BO[j][i] - W_BI[j][i]);
        }
    }
    // W = W_{AI,AO,BI} + W_{AI,BI,BO} - W_{AI,BI}
    for (int i=0; i<widthW; i++) {
        for (int j=i; j<widthW; j++) {
            prob.conZero.push_back(W_AI_AO_BI[i][j] + W_AI_BI_BO[i][j] - W_AI_BI[i][j]);
            prob.conZero.push_back(W_AI_AO_BI[j][i] + W_AI_BI_BO[j][i] - W_AI_BI[j][i]);
        }
    }

    // Constraints on the As
    // tr_AO \sum_a A^{AO,AI}_{a,x} = I
    std::vector<std::vector<Polynomial<double>>> A_0_0_AI = partialTrace(As[0], 2, d, 1);
    std::vector<std::vector<Polynomial<double>>> A_0_1_AI = partialTrace(As[1], 2, d, 1);
    std::vector<std::vector<Polynomial<double>>> A_1_0_AI = partialTrace(As[2], 2, d, 1);
    std::vector<std::vector<Polynomial<double>>> A_1_1_AI = partialTrace(As[3], 2, d, 1);
    std::vector<std::vector<Polynomial<double>>> A_0_sum = A_0_0_AI + A_0_1_AI;
    std::vector<std::vector<Polynomial<double>>> A_1_sum = A_1_0_AI + A_1_1_AI;
    for (int i=0; i<widthA; i++) {
        for (int j=i; j<widthA; j++) {
            if (i == j) {
                prob.conZero.push_back(A_0_sum[i][j] - Polynomial<double>(numVars, 1));
                prob.conZero.push_back(A_0_sum[j][i] - Polynomial<double>(numVars, 1));
                prob.conZero.push_back(A_1_sum[i][j] - Polynomial<double>(numVars, 1));
                prob.conZero.push_back(A_1_sum[j][i] - Polynomial<double>(numVars, 1));
            } else {
                prob.conZero.push_back(A_0_sum[i][j]);
                prob.conZero.push_back(A_0_sum[j][i]);
                prob.conZero.push_back(A_1_sum[i][j]);
                prob.conZero.push_back(A_1_sum[j][i]);
            }
        }
    }

    // Constraints on the Bs
    // tr_BO \sum_b B^{BO,BI}_{b,x} = I
    std::vector<std::vector<Polynomial<double>>> B_0_0_BI = partialTrace(Bs[0], 2, d, 1);
    std::vector<std::vector<Polynomial<double>>> B_0_1_BI = partialTrace(Bs[1], 2, d, 1);
    std::vector<std::vector<Polynomial<double>>> B_1_0_BI = partialTrace(Bs[2], 2, d, 1);
    std::vector<std::vector<Polynomial<double>>> B_1_1_BI = partialTrace(Bs[3], 2, d, 1);
    std::vector<std::vector<Polynomial<double>>> B_0_sum = B_0_0_BI + B_0_1_BI;
    std::vector<std::vector<Polynomial<double>>> B_1_sum = B_1_0_BI + B_1_1_BI;
    for (int i=0; i<widthB; i++) {
        for (int j=i; j<widthB; j++) {
            if (i == j) {
                prob.conZero.push_back(B_0_sum[i][j] - Polynomial<double>(numVars, 1));
                prob.conZero.push_back(B_0_sum[j][i] - Polynomial<double>(numVars, 1));
                prob.conZero.push_back(B_1_sum[i][j] - Polynomial<double>(numVars, 1));
                prob.conZero.push_back(B_1_sum[j][i] - Polynomial<double>(numVars, 1));
            } else {
                prob.conZero.push_back(B_0_sum[i][j]);
                prob.conZero.push_back(B_0_sum[j][i]);
                prob.conZero.push_back(B_1_sum[i][j]);
                prob.conZero.push_back(B_1_sum[j][i]);
            }
        }
    }

    // Objective TODO
    Polynomial<double> p_0_0_0_0 = trace(tensor(As[0], Bs[0]) * W);
    Polynomial<double> p_0_0_0_1 = trace(tensor(As[0], Bs[1]) * W);
    Polynomial<double> p_0_0_1_0 = trace(tensor(As[0], Bs[2]) * W);
    Polynomial<double> p_0_0_1_1 = trace(tensor(As[0], Bs[3]) * W);
    Polynomial<double> p_0_1_0_0 = trace(tensor(As[1], Bs[0]) * W);
    Polynomial<double> p_0_1_0_1 = trace(tensor(As[1], Bs[1]) * W);
    Polynomial<double> p_0_1_1_0 = trace(tensor(As[1], Bs[2]) * W);
    Polynomial<double> p_0_1_1_1 = trace(tensor(As[1], Bs[3]) * W);
    Polynomial<double> p_1_0_0_0 = trace(tensor(As[2], Bs[0]) * W);
    Polynomial<double> p_1_0_0_1 = trace(tensor(As[2], Bs[1]) * W);
    Polynomial<double> p_1_0_1_0 = trace(tensor(As[2], Bs[2]) * W);
    Polynomial<double> p_1_0_1_1 = trace(tensor(As[2], Bs[3]) * W);
    Polynomial<double> p_1_1_0_0 = trace(tensor(As[3], Bs[0]) * W);
    Polynomial<double> p_1_1_0_1 = trace(tensor(As[3], Bs[1]) * W);
    Polynomial<double> p_1_1_1_0 = trace(tensor(As[3], Bs[2]) * W);
    Polynomial<double> p_1_1_1_1 = trace(tensor(As[3], Bs[3]) * W);
    // GYNI = p(a = y, b = x)
    prob.obj = p_0_0_0_0 + p_0_1_1_0 + p_1_0_0_1 + p_1_1_1_1;

    // If verbose, output the problem
    if (verbosity > 0) {
        std::cout << "Problem:" << std::endl;
        std::cout << prob << std::endl;
    }

	// Optimize using branch and bound
	//auto res2 = prob.optimize(level, verbosity, maxIters);
	//if (!benchmark) {
		//std::cout << std::endl;
		//std::cout << "From branch and bound:" << std::endl;
		//std::cout << res2.first << " " << res2.second << std::endl;
		//std::cout << std::endl;
	//}

	return 0;

}
