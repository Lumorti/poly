import numpy as np

# Output with less precision
np.set_printoptions(precision=3, suppress=True)

# Define the dimension d
d = 4

# Create the X operator
def genX():
    X = np.zeros((d, d), dtype=complex)
    for i in range(d - 1):
        X[i + 1, i] = 1
    X[0, d - 1] = 1
    return X
X = genX()

# Create the Y operator
def genY(theta):
    Y = np.zeros((d, d), dtype=complex)
    omega = np.exp(theta * 2.0 * np.pi * 1j / float(d))
    for i in range(d - 1):
        Y[i + 1, i] = 1j * np.power(omega, i)
    Y[0, d - 1] = 1j * np.power(omega, d - 1)
    return Y
Y = genY(1)

# Create the Z operator
def genZ(theta):
    Z = np.zeros((d, d), dtype=complex)
    omega = np.exp(theta * 2.0 * np.pi * 1j / float(d))
    for i in range(d):
        Z[i, i] = np.power(omega, i)
    return Z
Z = genZ(1)

# The operator pool
operatorPool = [genX(), genZ(1), genZ(2)]
stringPool = ["X", "Z1", "Z2"]
numOps = len(operatorPool)

# The identity
I = np.eye(d, dtype=complex)

# Given a number, using ternary representation to generate the matrix
def matrixFromNumber(num):
    M = I
    string = ""
    while num > 0:
        rem = num % numOps
        M = operatorPool[rem] @ M
        string = stringPool[rem] + string
        num = num // numOps
    print(string)
    return M

minOverall = 1000
maxNum = 30
# for matNum1 in range(1, maxNum):
for matNum1 in range(1):
    # for matNum2 in range(1, maxNum):
    for matNum2 in range(1):
        # for matNum3 in range(1, maxNum):
        for matNum3 in range(1):
            # for matNum4 in range(1, maxNum):
            for matNum4 in range(1):

                # Define the M matrix list
                M = []
                M.append(Z)
                M.append(X)
                M.append(Z@X)
                M.append(Z@Z@X)
                # M.append(Z@(Z@(Z@X)))
                # M.append(Z@(Z@(Z@(X@(Z@(Z@Z))))))
                # M.append(matrixFromNumber(matNum1))
                # M.append(matrixFromNumber(matNum2))
                # M.append(matrixFromNumber(matNum3))
                # M.append(matrixFromNumber(matNum4))
                n = len(M)

                # Print the eigenvectors
                # print()
                eigenvecs = [None] * n
                for i in range(n):
                    eigenvecs[i] = np.linalg.eig(M[i])[1]
                    # for j in range(d):
                        # print(eigenvecs[i][:, j])

                # Calculate the inner products
                maxError = 0
                for i in range(n):
                    for j in range(i+1, n):
                        for k in range(d):
                            for l in range(d):
                                innerProduct = np.dot(eigenvecs[i][:, k], eigenvecs[j][:, l])
                                innerProduct = np.abs(innerProduct) - 1.0 / np.sqrt(d)
                                if innerProduct > maxError:
                                    maxError = innerProduct

                print("Max error: ", maxError)
                minOverall = min(minOverall, maxError)

print("Min overall: ", minOverall)

