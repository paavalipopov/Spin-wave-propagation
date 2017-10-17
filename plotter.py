from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
f = open('trash', 'w')
X, Y, Z= np.loadtxt("all results", delimiter='\t', usecols=(0, 1, 2), unpack=True, skiprows=1)
f.write(X, Y, Z)
X = np.arange(-5, 5, 0.05)
Y = np.arange(-5, 5, 0.05)
X, Y = np.meshgrid(X, Y)
Z = np.cos(np.sqrt(X**2 + Y**2))


# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

plt.show()