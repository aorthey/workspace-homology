import pickle
import numpy as np

from pylab import *
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

X = pickle.load( open( "data/cspaceX.dat", "rb" ) )
Y = pickle.load( open( "data/cspaceY.dat", "rb" ) )
Z = pickle.load( open( "data/cspaceZ.dat", "rb" ) )

fig=figure(1)
fig.clf()
ax = fig.gca(projection='3d')
xi = linspace(min(X), max(X))
yi = linspace(min(Y), max(Y))
zi = griddata(X, Y, Z, xi, yi)
xim, yim = meshgrid(xi, yi)

#ax.plot_surface(X,Y,Z,alpha=0.5)
#ax.scatter(X,Y,Z,marker='o',c='r')
ax.plot_surface(xim,yim,zi,rstride=1, cstride=1, cmap=cm.jet,\
                linewidth=0.2,vmax=2.)
#ax.contourf3D(xi, yi, zi)
#ax.contour3D(xi, yi, zi)
#ax.plot_wireframe(X,Y,Z)

plt.show()
