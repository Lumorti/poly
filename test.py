#!/usr/bin/python3
import math
from sympy import nonlinsolve, nsolve, solve, groebner, Symbol, sympify, poly, I, simplify, expand

# System params
d = 2
n = 2
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

for i in range(len(F)):
    print(F[i])
    # for j in range(numVars):
        # F[i] = F[i].subs(xi[j], ideal[j])
    # print(F[i])

# res = groebner(F, *xi, order='lex', method='f5b')
# print(res)

res = nsolve(F, xi, guess)
# res = nonlinsolve(F, xi)
print(res)

# for i in range(numVars):
    # if len(F[i].free_symbols) == 0:
        # print("no vars left in: ", F[i])
        # continue
    # varToRemove = next(iter(F[i].free_symbols))
    # print("removing var: ", varToRemove)
    # F[i] = expand(F[i])
    # print(F[i])
    # varSolved = solve(F[i], varToRemove)
    # if len(varSolved) > 0:
        # for j in range(len(F)):
            # F[j] = F[j].subs(varToRemove, varSolved[0])

# x4 = solve(F[-1], xi[4])
# print(F[-1])
# print()
# print(x4[0])
# print()
# print(expand(x4[0]))

# print()
# print(F)

