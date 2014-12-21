import pickle
import numpy as np

from pylab import *
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

XX = pickle.load( open( "data/cspaceX.dat", "rb" ) )
YY = pickle.load( open( "data/cspaceY.dat", "rb" ) )
ZZ = pickle.load( open( "data/cspaceZ.dat", "rb" ) )

fig=figure(1)
fig.clf()

N=len(XX)/4

t=pi/3
R = array([(cos(t), -sin(t),0), (sin(t), cos(t),0),(0,0,1)])
R = array([(cos(t), -sin(t),0), (sin(t), cos(t),0),(0,0,1)])
def rota(x):
        #return np.dot(R,x)
        return x

for i in range(0,4):
        X = XX[1+i*N:N+i*N]
        Y = YY[1+i*N:N+i*N]
        Z = ZZ[1+i*N:N+i*N]
        D = np.array([X,Y,Z])
        Drot = np.apply_along_axis(rota, axis=1,arr=D.T)
        X=Drot[:,0]
        Y=Drot[:,1]
        Z=Drot[:,2]

        X = X+np.random.normal((X))/100000
        Y = Y+np.random.normal((X))/100000
        Z = Z+np.random.normal((X))/100000

        ax = fig.gca(projection='3d')

        xi = linspace(min(X), max(X),80)
        yi = linspace(min(Y), max(Y),20)
        zi = griddata(X, Y, Z, xi, yi)

        xim, yim = meshgrid(xi, yi)

        #ax.scatter(X,Y,Z,marker='o',c='r',s=5)
        ax.plot_surface(xim,yim,zi,rstride=1, cstride=1, cmap=cm.jet,\
                        linewidth=0.001)
        #ax.contourf3D(xi, yi, zi)
        #ax.contour3D(xi, yi, zi)

plt.show()
