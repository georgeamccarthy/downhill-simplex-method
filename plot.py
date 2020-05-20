import matplotlib.pyplot as plt
import numpy as np
data = np.loadtxt("C:\\Users\\George\\Documents\\CCourse\\Coursework\\output.txt", delimiter = ',')
x = np.array(data[:,0])
F = np.array(data[:,1])
plt.plot(x, F)
plt.ylabel('y')
plt.xlabel(r'$x_{0}$')
# Changing the axis bounds zooms the graph.
plt.axis([-2, 2, 0, 909])
plt.draw()
plt.show()