import numpy as np

# Output with less precision
np.set_printoptions(precision=3, suppress=True)

# Define the dimension d
d = 3

# Given a list of vectors and eigenvalues, find the operator that has them as eigenvectors
def findOperator(eigenvecs, eigenvals):
    n = len(eigenvecs)
    M = np.zeros((d, d), dtype=complex)
    for i in range(n):
        M += eigenvals[i] * np.outer(eigenvecs[i], eigenvecs[i].conj())
    return M

# Given an operator, find the eigenvectors
def findEigenvectors(M):
    eigenvals, eigenvecs = np.linalg.eig(M)
    return eigenvecs.T

# Create the identity operator
def genI(d):
    return np.eye(d, dtype=complex)
I = genI(d)

# Create the X operator
def genX(d):
    X = np.zeros((d, d), dtype=complex)
    for i in range(d - 1):
        X[i + 1, i] = 1
    X[0, d - 1] = 1
    return X
X = genX(d)

# Create the Z operator
def genZ(d):
    Z = np.zeros((d, d), dtype=complex)
    omega = np.exp(2.0 * np.pi * 1j / float(d))
    for i in range(d):
        Z[i, i] = np.power(omega, i)
    return Z
Z = genZ(d)

# Function that checks if a list of matrices are mutually unbiased
def checkMutualUnbiasedness(Ms):
    for i in range(len(Ms)):
        for j in range(i+1, len(Ms)):
            for k in range(len(Ms[i])):
                for l in range(len(Ms[j])):
                    innerProd = np.abs(np.dot(Ms[i][k], Ms[j][l].conj()))
                    if np.abs(innerProd - 1.0/np.sqrt(d)) > 1e-8:
                        return False
    return True

# Dim 3
opsToFind = []
if d == 3:
    opsToFind.append(Z)
    opsToFind.append(X)
    opsToFind.append(X @ Z)
    opsToFind.append(X @ (Z @ Z))

    # Print each set of eigenvalues
    for i in range(len(opsToFind)):
        print("Operator " + str(i) + ":")
        eigenvals, eigenvecs = np.linalg.eig(opsToFind[i])
        print(eigenvals)
        print()

# Dim 4
elif d == 4:
    vals = [1, -1, 1j, -1j]
    M_0 = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
    M_1 = np.array([
        [0.5, 0.5, 0.5, 0.5],
        [0.5, 0.5, -0.5, -0.5],
        [0.5, -0.5, -0.5, 0.5],
        [0.5, -0.5, 0.5, -0.5]
    ])
    M_2 = 0.5 * np.array([
        [1, -1, -1j, -1j],
        [1, -1, 1j, 1j],
        [1, 1, 1j, -1j],
        [1, 1, -1j, 1j]
    ])
    M_3 = 0.5 * np.array([
        [1, -1j, -1j, -1],
        [1, -1j, 1j, 1],
        [1, 1j, 1j, -1],
        [1, 1j, -1j, 1]
    ])
    M_4 = 0.5 * np.array([
        [1, -1j, -1, -1j],
        [1, -1j, 1, 1j],
        [1, 1j, -1, 1j],
        [1, 1j, 1, -1j]
    ])
    opsToFind.append(findOperator(M_0, vals))
    opsToFind.append(findOperator(M_1, vals))
    opsToFind.append(findOperator(M_2, vals))
    opsToFind.append(findOperator(M_3, vals))
    opsToFind.append(findOperator(M_4, vals))

# Print each operator
for i in range(len(opsToFind)):
    print("Operator " + str(i) + ":")
    print(opsToFind[i])
    print()

# Find the eigenvectors of each operator
vecs = [findEigenvectors(op) for op in opsToFind]

# See if the eigenvectors are mutually unbiased
print("Unbiased: ", checkMutualUnbiasedness(vecs))

exit()

# Iterate over all combinations of sums of paulis and phases given
def sumOfPaulisAndPhases(num, paulis, phases):

    # Generate the combinations of paulis and phases
    combs = []
    for i in range(len(paulis)):
        for k in range(len(paulis)):
            for j in range(len(phases)):
                combs.append(phases[j] * np.kron(paulis[i], paulis[k]))

    # Start with a blank matrix
    matDim = len(paulis[0]) ** 2
    M = np.zeros((matDim, matDim), dtype=np.complex128)

    # If the number is 0, return the first
    if num == 0:
        return M + combs[0]

    # Iterate over the number in ternary representation
    while num > 0:
        rem = num % len(combs)
        M = M + combs[rem]
        num = num // len(combs)

    return M

# Same as the above function but with strings
def sumOfPaulisAndPhasesString(num, paulis, phases):
    
    # Generate the combinations of paulis and phases
    combs = []
    for i in range(len(paulis)):
        for k in range(len(paulis)):
            for j in range(len(phases)):
                combs.append(phases[j] + " " + paulis[i] + "(x)" + paulis[k])

    # Start with a blank string
    M = ""

    # If the number is 0, return the first
    if num == 0:
        return M + combs[0]

    # Iterate over the number in ternary representation
    while num > 0:
        rem = num % len(combs)
        M = M + combs[rem]
        num = num // len(combs)
        if num > 0:
            M = M + " + "

    return M

# Mat1 = 1108 (0.5+0.5j) X(x)X + (0.5-0.5j) X(x)I
# Mat2 = 1774 (-0.5+0.5j) Y(x)Z + (-0.5-0.5j) X(x)Y
# Mat3 =
# Mat4 =

mat1Inds = []
mat2Inds = []
mat3Inds = []
mat4Inds = []
X2 = genX(2)
for matNum in range(0, 3000):
    M = sumOfPaulisAndPhases(matNum, [genI(2), genX(2), genY(2), genZ(2)], [0.5+0.5j, 0.5-0.5j, -0.5+0.5j, -0.5-0.5j])
    MString = sumOfPaulisAndPhasesString(matNum, ["I", "X", "Y", "Z"], ["(0.5+0.5j)", "(0.5-0.5j)", "(-0.5+0.5j)", "(-0.5-0.5j)"])
    expanded = M
    if len(mat1Inds) == 0:
        if np.allclose(expanded, opToFind1):
            mat1Inds.append(matNum)
            print("Found mat1: ", matNum, MString)
    if len(mat2Inds) == 0:
        if np.allclose(expanded, opToFind2):
            mat2Inds.append(matNum)
            print("Found mat2: ", matNum, MString)
    if len(mat3Inds) == 0:
        if np.allclose(expanded, opToFind3):
            mat3Inds.append(matNum)
            print("Found mat3: ", matNum, MString)
    if len(mat4Inds) == 0:
        if np.allclose(expanded, opToFind4):
            mat4Inds.append(matNum)
            print("Found mat4: ", matNum, MString)

print()
print("Mat1: ", mat1Inds)
print("Mat2: ", mat2Inds)
print("Mat3: ", mat3Inds)
print("Mat4: ", mat4Inds)

