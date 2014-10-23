import numpy as np
import networkx as nx
from networkx import graphviz_layout

from src.polytope import Polytope
from src.walkable import WalkableSurface, WalkableSurfacesFromPolytopes, ProjectPolytopesDownInsideBox
from src.urdfparser import URDFtoPolytopes
from src.linalg import intersection
from src.linalg import distanceWalkableSurfaceWalkableSurface
from src.robotspecifications import ROBOT_FOOT_RADIUS, ROBOT_MAX_SLOPE, ROBOT_SPHERE_RADIUS
from src.robotspecifications import ROBOT_FOOT_HEIGHT
from src.robotspecifications import ROBOT_VOLUME_MIN_HEIGHT, ROBOT_VOLUME_MAX_HEIGHT
from plotter import Plotter
plot = Plotter()

pobjects = URDFtoPolytopes("wall.urdf")
wsurfaces = WalkableSurfacesFromPolytopes(pobjects)

wboxesLow=[]
for i in range(0,len(wsurfaces)):
        wboxesLow.append(wsurfaces[i].createBox(0, ROBOT_FOOT_HEIGHT))

###############################################################################
# Visualize Complete Scene
###############################################################################
#for i in range(0,len(pobjects)):
        #plot.polytopeFromVertices(pobjects[i].getVertexRepresentation(),(0,0,0,0.05))

###############################################################################
# Project Objects in WboxesLow down
###############################################################################

Wsurfaces_decomposed = []
for i in range(0,len(wsurfaces)):

        ap = wsurfaces[i].ap
        bp = wsurfaces[i].bp
        iObject = wsurfaces[i].iObject
        p = ProjectPolytopesDownInsideBox(pobjects, wsurfaces[i], wboxesLow[i])

        ###############################################################################
        # Project Objects in Upper Box above splitted surfaces Down
        ###############################################################################

        for j in range(0,len(p)):
                Wsplit = WalkableSurface.fromVertices(ap,bp,p[j],iObject)
                #plot.polytopeFromPolygonVertices3D(Wsplit.getVertexRepresentation(),0.1,(1,0,0,0.6))

                WsplitBox = Wsplit.createBox(ROBOT_FOOT_HEIGHT,ROBOT_VOLUME_MIN_HEIGHT)

                p2 = ProjectPolytopesDownInsideBox(pobjects, wsurfaces[i], WsplitBox)
                for k in range(0,len(p2)):
                        W2split = WalkableSurface.fromVertices(ap,bp,p2[k],iObject)
                        Wsurfaces_decomposed.append(W2split)

###############################################################################
# Visualize traversable surfaces
###############################################################################

for i in range(0,len(Wsurfaces_decomposed)):
        plot.polytopeFromPolygonVertices3D(Wsurfaces_decomposed[i].getVertexRepresentation(),0.2,(0,1,0,0.6))

###############################################################################
# Find topology of traversable surfaces
###############################################################################
N = len(Wsurfaces_decomposed)
WD = np.zeros((N,N))
WM = np.zeros((N,N))
for i in range(0,N):
        for j in range(i+1,N):
                WD[i,j]=WD[j,i]=distanceWalkableSurfaceWalkableSurface(Wsurfaces_decomposed[i], Wsurfaces_decomposed[j])
                WM[i,j]=WM[j,i]=(0 if WD[i,j]>0.1 else 1)

        WM[i,i]=1

print np.around(WD,2)
print WM
print np.around(WM,2)

###############################################################################
# Plan path on top of traversable surfaces
###############################################################################
## START AND GOAL CONTACTS
xstart = np.array((-0.5,2,0))
xgoal = np.array((0,-2,0))

G = nx.Graph()
for i in range(0,N):
        for j in range(0,N):
                if not WD[i,j]>0.01:
                        G.add_edge(i,j)

#nx.path.bidirectional_dijkstra(G,1,2)
plot.graphLayout(G)

plot.show()
