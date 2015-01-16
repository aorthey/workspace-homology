from pylab import *
from numpy import *
import networkx as nx
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt

### counterclock wise rotation
def rotFromRPY(tx,ty,tz):
        Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
        Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
        Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])
        return np.dot(Rx, np.dot(Ry, Rz))

def axisFromRPY(roll,pitch,yaw):
        x = cos(yaw)*cos(pitch)
        y = sin(yaw)*cos(pitch)
        z = sin(pitch)

        return np.array((x,y,z))

DEFAULT_EDGE_COLOR = (0,0,0,0.4)
class Plotter:
        def __init__(self):
                self.fig=figure(1)
                #self.graph=figure(2)
                #self.robot_fig = figure(3)
                #self.rax = self.robot_fig.gca(projection='3d')
                self.ax = self.fig.gca(projection='3d')

        def graphLayout(self, G):
                pos=nx.spring_layout(G)
                nx.draw_networkx_nodes(G, pos, node_color='r')
                nx.draw_networkx_edges(G, pos, edge_color='b', width=1.0,alpha=0.5)

        def walkableSurface(self ,Vin,thickness=0.1, fcolor=(0,0,0,0.1), ecolor=DEFAULT_EDGE_COLOR):
                self.polytopeFromPolygonVertices3D(Vin,thickness=thickness, \
                                fcolor=fcolor,ecolor=ecolor)

        def polytopeFromPolygonVertices3D(self,Vin,thickness=0.1, fcolor=(0,0,0,0.1), ecolor=DEFAULT_EDGE_COLOR):
                Vz = np.zeros((len(Vin),1))
                Vz.fill(thickness/2)
                Vup = column_stack((Vin[:,0:2],thickness/2+Vin[:,2]))
                Vdown = column_stack((Vin[:,0:2],-thickness/2+Vin[:,2]))
                V = np.vstack((Vdown, Vup))

                self.polytopeFromVertices(V, fcolor, ecolor)

        def polytopeFromVertices2D(self,Vin,thickness=0.1, fcolor=(0,0,0,0.1), ecolor=DEFAULT_EDGE_COLOR):
                Vz = np.zeros((len(Vin),1))
                Vz.fill(thickness/2)
                Vup = column_stack((Vin[:,0:2],Vz))
                Vdown = column_stack((Vin[:,0:2],-Vz))
                V = np.vstack((Vdown, Vup))

                self.polytopeFromVertices(V, fcolor, ecolor)

        def polytopeFromVertices(self,V,fcolor=(0,0,0,0.1), ecolor=DEFAULT_EDGE_COLOR):
                self.hull = ConvexHull(V)
                faces = []
                for ia, ib, ic in self.hull.simplices:
                     faces.append(V[[ia, ib, ic]])
                items = Poly3DCollection(faces, facecolors=[fcolor],
                                edgecolors=[ecolor])
                self.ax.add_collection(items)
                #self.ax.scatter(V[:,0], V[:,1], V[:,2], 'r*')

        def capsule(self, radius, length, xyz, rpy):
                x = xyz[0]
                y = xyz[1]
                z = xyz[2]
                rr = rpy[0]
                rp = rpy[1]
                ry = rpy[2]
                R = rotFromRPY(rr,rp,ry)

                u, v = np.mgrid[0:2*np.pi:10j, 0:2*np.pi:20j]
                xx=radius*np.sin(u)*cos(v)+x
                yy=radius*np.sin(u)*sin(v)+y
                zz=radius*np.cos(u)+z
                print xx
                self.rax.plot_wireframe(xx,yy,zz, color="b")

        def point(self, x, size=100, color=(1,0,0,1)):
                self.ax.scatter(x[0],x[1],x[2], s=size, c=color )

        def showRobot(self):
                self.rax.set_xlabel('x-axis')
                self.rax.set_ylabel('y-axis')
                self.rax.set_zlabel('z-axis')
                self.robot_fig.show()

        def showEnvironment(self):
                self.ax.set_xlim(-2, 2)
                self.ax.set_ylim(-2, 2)
                self.ax.set_zlim(0, 1.5)
                self.fig.show()

        def pause(self,t):
                plt.pause(t)

        def clear(self):
                self.fig.clf()
        def show(self):
                plt.show()

        def set_view(self,azim,elev):
                self.ax.view_init(elev=elev, azim=azim)

        def line(self, L,style='-r',lw=5.0):
                x=[]
                y=[]
                z=[]
                for i in range(0,len(L)):
                        x.append(L[i][0])
                        y.append(L[i][1])
                        z.append(L[i][2])

                x = np.array(x).flatten()
                y = np.array(y).flatten()
                z = np.array(z).flatten()

                self.ax.plot(x, y, z, style,linewidth=lw)

        def lines(self, L,style='o-r'):
                N = L.shape[0]
                M = L.shape[1]
                for i in range(0,N):
                        l = L[i].T
                        x = l[:][0]
                        y = l[:][1]
                        z = l[:][2]
                        self.ax.plot(x, y, z, style,linewidth=8.0)


        def display_capsule2(self, Xc, rpy, length, r, color='r'):
                ## modified code from https://github.com/roboptim/roboptim-analysis/blob/master/bin/capsule_display
                #################################################################
                ## compute P1 from P0, the length and the euler angles rpy
                #################################################################
                V = axisFromRPY(rpy[0],rpy[1],rpy[2])

                P0 = Xc - (length/2)*V/np.linalg.norm(V)
                P1 = Xc + (length/2)*V/np.linalg.norm(V)

                print P0,P1

                mid = 0.5 * (P0 + P1)
                dir = np.array([0., 0., 1.])
                norm = np.linalg.norm(P1-P0)
                if norm > 1e-3:
                    dir = (P1-P0)/norm

                x0 = [1, 0, 0]
                y0 = [0, 1, 0]
                z0 = [0, 0, 1]

                # Capsule is prepared in the following basis: B1 = (x1, z1, dir)
                # Let's take x1 = z0 x dir (or y0 x dir)
                x1 = np.cross(z0, dir)
                if np.linalg.norm(x1) < 1e-6:
                    x1 = np.cross(y0, dir)
                x1 /= np.linalg.norm(x1)
                y1 = np.cross(dir, x1)
                y1 /= np.linalg.norm(y1)

                # We note R the rotation matrix from B = (x, y, z) to B'
                R = np.identity(3)
                R[0, :] = x1
                R[1, :] = y1
                R[2, :] = dir
                assert abs(np.linalg.det(R) - 1) < 1e-6

                def rotate_xyz(rot, x, y, z):
                    l1 = x.shape[0]
                    l2 = x.shape[1]
                    for i in xrange(l1):
                            for j in xrange(l2):
                                    v = np.array([x[i, j], y[i, j], z[i, j]])
                                    x[i, j] = R[:, 0].dot(v)
                                    y[i, j] = R[:, 1].dot(v)
                                    z[i, j] = R[:, 2].dot(v)
                    return x, y, z

                def translate_xyz(trans, x, y, z):
                    l1 = x.shape[0]
                    l2 = x.shape[1]
                    for i in xrange(l1):
                            for j in xrange(l2):
                                    x[i, j] = trans[0] + x[i, j]
                                    y[i, j] = trans[1] + y[i, j]
                                    z[i, j] = trans[2] + z[i, j]
                    return x, y, z

                n = 100
                dtheta = 2*pi/n
                dz = norm/n
                [theta, z] = np.mgrid[-pi:pi+dtheta:dtheta, -norm/2-dz:norm/2+dz:dz]
                x = r*cos(theta)
                y = r*sin(theta)
                x, y, z = rotate_xyz(R.T, x, y, z)
                x, y, z = translate_xyz(mid, x, y, z)
                cylinder = self.rax.plot_surface(x, y, z, color=color)


                [theta, phi] = np.mgrid[-pi:pi+dtheta:dtheta, -dtheta:pi/2+dtheta:dtheta]
                x = r*cos(theta)*sin(phi)
                y = r*sin(theta)*sin(phi)
                z = -norm/2 - r*cos(phi)
                x, y, z = rotate_xyz(R.T, x, y, z)
                x, y, z = translate_xyz(mid, x, y, z)
                halfsphere0 = self.rax.plot_surface(x, y, z, color=color)

                x = r*cos(theta)*sin(phi)
                y = r*sin(theta)*sin(phi)
                z = norm/2 + r*cos(phi)
                x, y, z = rotate_xyz(R.T, x, y, z)
                x, y, z = translate_xyz(mid, x, y, z)
                halfsphere1 = self.rax.plot_surface(x, y, z, color=color)

                #axes = zip (P0,P1)
                #endPoints = mlab.points3d(axes[0],
                #                          axes[1],
                #                          axes[2],
                #                          scale_factor=0.01,
                #                          color=color)
