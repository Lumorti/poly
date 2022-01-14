#!/usr/bin/python3
import math
from sympy import nonlinsolve, nsolve, solve, groebner, Symbol, sympify, poly, I, simplify, expand
import cvxpy as cp
import numpy as np

# System params
d = 4
n = 4
numVars = 2*n*d*d

# The list of variables
def nameFromInd(i):
    if i % 2 == 0:
        return "a" + str(int(i/2))
    else:
        return "b" + str(int((i-1)/2))
xi = [Symbol(nameFromInd(i), real=True) for i in range(numVars)]
rt2 = 1.0 / math.sqrt(2.0)
ideal = [
            1.0, 0.0,  0.0, 0.0,
            0.0, 0.0,  1.0, 0.0,
            rt2, 0.0,  rt2, 0.0,
            rt2, 0.0, -rt2, 0.0
        ]
guess = [1 for i in range(numVars)]

# The list of polynomials
F = []

# For each set
for i in range(n):
    for j in range(d):
        for k in range(d):

            # Normalisation
            if j == k:
                start = i*d*d*2 + j*d*2
                expr = sympify(-1)
                for l in range(0, 2*d, 2):
                    expr += xi[start+l]*xi[start+l] + xi[start+l+1]*xi[start+l+1]
                F.append(expr)

# For each set
for i in range(n):
    for j in range(d):
        for k in range(d):

            # Orthogonality
            if j != k:
                exprReal = sympify(0)
                exprImag = sympify(0)
                startj = i*d*d*2 + j*d*2
                startk = i*d*d*2 + k*d*2
                for l in range(0, 2*d, 2):
                    exprReal += xi[startj+l]*xi[startk+l] + xi[startj+l+1]*xi[startk+l+1]
                    exprImag += xi[startj+l]*xi[startk+l+1] - xi[startj+l+1]*xi[startk+l]
                F.append(exprReal)
                F.append(exprImag)

# For each pair of sets
for i in range(n):
    for j in range(i+1, n):

        # For each combination of vectors
        for k in range(d):
            for l in range(d):

                startk = i*d*d*2 + k*d*2
                startl = j*d*d*2 + l*d*2
                expr = sympify(0)
                exprConj = sympify(0)
                for m in range(0, 2*d, 2):
                    expr += (xi[startk+m]-I*xi[startk+m+1])*(xi[startl+m]+I*xi[startl+m+1])
                    exprConj += (xi[startk+m]+I*xi[startk+m+1])*(xi[startl+m]-I*xi[startl+m+1])
                full = exprConj*expr + sympify(-1/d)
                F.append(expand(full))

# Combine into a single polynomial
toOpt = sympify(0)
for i in range(len(F)):
    print(F[i])
    toOpt += F[i]*F[i]
toOpt = poly(toOpt)

# Extract the coefficients of each term
coeffs = toOpt.coeffs()
monoms = toOpt.monoms()
print("num unique = ", len(coeffs))

# Generate sparse matrix for the objective TODO
sparseVals = []
sparseX = []
sparseY = []

print(sparseVals)
print(sparseX)
print(sparseY)

# Run the SDP
# C = cp.spmatrix(sparseVals, sparseX, sparseY)
# X = cp.Variable((n, n), symmetric=True)
# constraints = [X >> 0, X[0] == coeffs[-1]]
# prob = cp.Problem(cp.Minimize(cp.trace(C @ X)), constraints)
# prob.solve()


