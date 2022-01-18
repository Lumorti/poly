#!/usr/bin/python3
import math
from sympy import Symbol, sympify, poly, I, simplify, expand, Matrix
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
print("Generating equations...")
for i in range(n):
    for j in range(d):
        for k in range(d):

            # Normalisation
            if j == k:
                start = i*d*d*2 + j*d*2
                expr = sympify(-1)
                for l in range(0, 2*d, 2):
                    expr += xi[start+l]*xi[start+l] + xi[start+l+1]*xi[start+l+1]
                F.append(poly(expr))

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
                F.append(poly(exprReal))
                F.append(poly(exprImag))

# For each pair of sets
for i in range(n):
    for j in range(i+1, n):

        # For each combination of vectors
        for k in range(d):
            for l in range(d):

                # Their inner product should equal 1/d
                startk = i*d*d*2 + k*d*2
                startl = j*d*d*2 + l*d*2
                expr = sympify(0)
                exprConj = sympify(0)
                for m in range(0, 2*d, 2):
                    expr += (xi[startk+m]-I*xi[startk+m+1])*(xi[startl+m]+I*xi[startl+m+1])
                    exprConj += (xi[startk+m]+I*xi[startk+m+1])*(xi[startl+m]-I*xi[startl+m+1])
                full = exprConj*expr + sympify(-1/d)
                F.append(poly(expand(full)))

# Extract the coefficients for each equation
print("Extracting coefficients...")
coeffs = []
monoms = []
gens = []
for eqn in F:
    coeffs.append([int(a) for a in eqn.coeffs()])
    monoms.append([list(a) for a in eqn.monoms()])
    gens.append([a for a in eqn.gens])

# Get an ordering for the monomials
print("Getting monomial ordering...")
ind = 0
uniqueMonoms = {}
fullMonoms = []
for i in range(len(monoms)):
    newMonoms = []
    for j in range(len(monoms[i])):
        newMonom = ""
        for k in range(len(monoms[i][j])):
            if monoms[i][j][k] > 0:
                newMonom += str(monoms[i][j][k]) + str(gens[i][k]) + ","
        newMonoms.append(newMonom)
        if newMonom not in uniqueMonoms.keys():
            uniqueMonoms[newMonom] = ind
            ind += 1
    fullMonoms.append(newMonoms)

# Generate the coefficient matrix
print("Generating coefficient matrix...")
fullCoeffs = np.zeros((len(coeffs), ind))
for i in range(len(coeffs)):
    for j in range(len(coeffs[i])):
        fullCoeffs[i,uniqueMonoms[fullMonoms[i][j]]] = coeffs[i][j]

# Row reduce
print("Calculating matrix rank...")
rank = np.linalg.matrix_rank(fullCoeffs)
print("rank = ", rank)
print("num vars = ", numVars)


