from pylab import *
from numpy import *
import networkx as nx
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
class Plotter:
        def __init__(self):
                self.fig=figure(1)
                self.graph=figure(2)
                self.ax = self.fig.gca(projection='3d')

        def graphLayout(self, G):
                pos=nx.spring_layout(G)
                nx.draw_networkx_nodes(G, pos, node_color='r')
                nx.draw_networkx_edges(G, pos, edge_color='b', width=1.0,alpha=0.5)

        def polytopeFromPolygonVertices3D(self,Vin,thickness=0.1, fcolor=(0,0,0,0.1), ecolor=(1,1,1,0.4)):
                Vz = np.zeros((len(Vin),1))
                Vz.fill(thickness/2)
                Vup = column_stack((Vin[:,0:2],thickness/2+Vin[:,2]))
                Vdown = column_stack((Vin[:,0:2],-thickness/2+Vin[:,2]))
                V = np.vstack((Vdown, Vup))

                self.polytopeFromVertices(V, fcolor, ecolor)

#        def polytopeFromPolygonVertices(self,Vin,thickness=0.1, fcolor=(0,0,0,0.1), ecolor=(1,1,1,0.4)):
#                Vz = np.zeros((len(Vin),1))
#                Vz.fill(thickness/2)
#                Vup = column_stack((Vin,Vz))
#                Vdown = column_stack((Vin,-Vz))
#                V = np.vstack((Vdown, Vup))
#
#                self.polytopeFromVertices(V, fcolor, ecolor)

        def polytopeFromVertices2D(self,Vin,thickness=0.1, fcolor=(0,0,0,0.1), ecolor=(1,1,1,0.4)):
                Vz = np.zeros((len(Vin),1))
                Vz.fill(thickness/2)
                Vup = column_stack((Vin[:,0:2],Vz))
                Vdown = column_stack((Vin[:,0:2],-Vz))
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
                plt.show()

