import numpy as np
import matplotlib.pyplot as plt

# Load the points from the csv
points = np.loadtxt('temp2.csv', delimiter=',')

# Indices to use 
ind1 = 10
ind2 = 203
ind3 = 305

# Plot these points in 3d
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(points[:, ind1], points[:, ind2], points[:, ind3], c='r', marker='o')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
plt.show()


