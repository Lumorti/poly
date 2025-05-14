import numpy as np
import matplotlib.pyplot as plt

# Load the array from the csv
data = np.loadtxt('hessian.csv', delimiter=',')

# Show it using imshow
plt.imshow(data, cmap='gray', interpolation='nearest')
plt.colorbar()
plt.title('Data Visualization')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show()

