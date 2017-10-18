from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')


X, Y, Z= np.loadtxt("all results", delimiter='\t', usecols=(0, 1, 2), unpack=True, skiprows=1)

# Plot the surface.
surf = ax.plot_trisurf(X, Y, Z, linewidth=0.0005, antialiased=True, cmap=cm.autumn)

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()