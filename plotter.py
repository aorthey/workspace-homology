from pylab import *
from numpy import *
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
class Plotter:
        def __init__(self):
                self.fig=figure()
                self.ax = self.fig.gca(projection='3d')
        def polytope(self,A,b):
                #NIY
                return
        def points(self,V):
                self.hull = ConvexHull(V)
                faces = []
                for ia, ib, ic in self.hull.simplices:
                     faces.append(V[[ia, ib, ic]])
                items = Poly3DCollection(faces, facecolors=[(0, 0, 0, 0.1)])
                self.ax.add_collection(items)
                self.ax.scatter(V[:,0], V[:,1], V[:,2], 'o')

        def show(self):
                show()


