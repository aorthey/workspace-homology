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

        def polytopeFromPolygonVertices(self,Vin,thickness=0.1, fcolor=(0,0,0,0.1), ecolor=(1,1,1,0.4)):
                Vz = np.zeros((len(Vin),1))
                Vz.fill(thickness/2)
                Vup = column_stack((Vin,Vz))
                Vdown = column_stack((Vin,-Vz))
                V = np.vstack((Vdown, Vup))

                self.polytopeFromVertices(V, fcolor, ecolor)

        def polytopeFromVertices(self,V,fcolor=(0,0,0,0.1), ecolor=(1,1,1,0.4)):
                self.hull = ConvexHull(V)
                faces = []
                for ia, ib, ic in self.hull.simplices:
                     faces.append(V[[ia, ib, ic]])
                items = Poly3DCollection(faces, facecolors=[fcolor],
                                edgecolors=[ecolor])
                self.ax.add_collection(items)
                self.ax.scatter(V[:,0], V[:,1], V[:,2], 'r*')

        def show(self):
                show()

