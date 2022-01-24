#!/usr/bin/python3
import math
from sympy import conjugate, zeros, Symbol, sympify, poly, I, simplify, expand, Matrix
import numpy as np
import sys

# System params
d = 2
n = 2
kIndep = 1
if len(sys.argv) > 1:
    d = int(sys.argv[1])
if len(sys.argv) > 2:
    n = int(sys.argv[2])
if len(sys.argv) > 3:
    kIndep = int(sys.argv[3])
numVars = n*d*d

# The list of variables
def nameFromInd(i):
    return "v" + str(i)
xi = [Symbol(nameFromInd(i)) for i in range(numVars)]
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
                start = i*d*d + j*d
                expr = sympify(-1)
                for l in range(d):
                    expr += conjugate(xi[start+l])*xi[start+l]
                F.append(expr)

# For each set
for i in range(n):

    # For each combination in this set
    for j in range(d):
        for k in range(j+1, d):

            # Orthogonality
            if j != k:
                expr = sympify(0)
                startj = i*d*d + j*d
                startk = i*d*d + k*d
                for l in range(d):
                    expr += conjugate(xi[startj+l])*xi[startk+l]
                F.append(expr)

# For each pair of sets
for i in range(n):
    for j in range(i+1, n):

        # For each combination of vectors
        for k in range(d):
            for l in range(d):

                # Their inner product should equal 1/d
                startk = i*d*d + k*d
                startl = j*d*d + l*d
                expr = sympify(0)
                exprConj = sympify(0)
                for m in range(d):
                    expr += conjugate(xi[startk+m])*xi[startl+m]
                full = conjugate(expr)*expr + sympify(-1/d)
                F.append(expand(full))

numOrig = len(F)

# for eq in F:
    # print(eq)

# For higher k
print("Generating higher k extension...")
kPoly = []
if kIndep == 2:
    for x in xi:
        kPoly.append(x)
elif kIndep == 3:
    for x in xi:
        kPoly.append(x)
    for x in xi:
        for y in xi:
            kPoly.append(x*y)
elif kIndep == 4:
    for x in xi:
        kPoly.append(x)
    for x in xi:
        for y in xi:
            kPoly.append(x*y)
    for x in xi:
        for y in xi:
            for z in xi:
                kPoly.append(x*y*z)
newEqns = []
for pol in kPoly:
    for eqn in F:
        newEqns.append(eqn*pol)
F.extend(newEqns)

# Add the complex conjugates of each
newEqns2 = []
for eq in F:
    newEqns2.append(conjugate(eq))
F.extend(newEqns2)

# Extract the coefficients for each equation
print("Extracting coefficients for {} equations...".format(len(F)))
coeffs = []
monoms = []
gens = []
for i, eqn in enumerate(F):
    print("{} / {}".format(i, len(F)))
    eqn = poly(eqn)
    coeffs.append(eqn.coeffs())
    monoms.append([list(a) for a in eqn.monoms()])
    gens.append(eqn.gens)

# Get an ordering for the monomials
print("Getting monomial ordering...")
ind = 1
uniqueMonoms = {"": 0}
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
print("Matrix is {} by {}".format(len(coeffs), ind))
print("Writing coefficient matrix...")
totalWrites = 0
with open("matrices/d{}n{}k{}.csv".format(d,n,kIndep), "w") as f:
    f.write("{:d} {:d} {:d} {:d}\n".format(len(coeffs), ind, numOrig, numVars))
    for i in range(len(coeffs)):
        for j in range(len(coeffs[i])):
            f.write("{:d} {:d} {:.2f}\n".format(i, uniqueMonoms[fullMonoms[i][j]], float(coeffs[i][j])))
            totalWrites += 1
print("{} unique elements in matrix".format(totalWrites))

# Row reduce TODO
# print("Reducing matrix...")
# for i in range(1):

    # # Keep track of which rows have already been used (no swaps allowed)
    # rowUsed = [False for i in range(matY)]
    # rowUsed[i] = True

    # # For each element of the row, try to make it zero
    # for j in range(matX):

        # # Find a pivot row
        # piv = -1
        # for k in range(matY):
            # if fullCoeffs[k,j] != 0 and not rowUsed[k]:
                # piv = k
                # rowUsed[k] = True
                # break

        # # If there's no valid pivot, stop
        # if piv == -1:
            # continue

        # # Make this have the only non-zero in the column
        # for k in range(matY):
            # if k != piv:
                # factor = -fullCoeffs[k,j] / fullCoeffs[piv,j]
                # fullCoeffs[k,:] = fullCoeffs[k,:] + factor*fullCoeffs[piv,:]

# zeroRows = np.where(~fullCoeffs[0:numOrig+1,:].any(axis=1))[0]





