import numpy as np
import networkx as nx
from networkx import graphviz_layout

from src.polytope import Polytope
from src.walkable import WalkableSurface, WalkableSurfacesFromPolytopes, ProjectPolytopesDownInsideBox
from src.walkable import getStartGoalWalkableSurfaces
from src.urdfparser import URDFtoPolytopes
from src.linalg import intersection
from src.linalg import distanceWalkableSurfaceWalkableSurface
from src.linalg import distancePointPolytope
from src.linalg import distanceWalkableSurfacePolytope
from src.robotspecifications import ROBOT_FOOT_RADIUS, ROBOT_MAX_SLOPE, ROBOT_SPHERE_RADIUS
from src.robotspecifications import ROBOT_FOOT_HEIGHT
from src.robotspecifications import ROBOT_VOLUME_MIN_HEIGHT, ROBOT_VOLUME_MAX_HEIGHT
from src.robotspecifications import ROBOT_MAX_HEAD_DISPLACEMENT, ROBOT_HEAD_SIZE
from src.pathoptimizer import optimizePath
from plotter import Plotter

## START AND GOAL CONTACTS OF THE ROBOT
xstart = np.array((-0.5,2,0))
xgoal = np.array((0,-2,0))
env_fname = "wall.urdf"

plot = Plotter()

pobjects = URDFtoPolytopes(env_fname)
wsurfaces = WalkableSurfacesFromPolytopes(pobjects)

###############################################################################
# create footBox, box over Sip with height of one foot
###############################################################################
footBox=[]
for i in range(0,len(wsurfaces)):
        footBox.append(wsurfaces[i].createBox(0, ROBOT_FOOT_HEIGHT))

###############################################################################
# Visualize Complete Scene
###############################################################################
for i in range(0,len(pobjects)):
        plot.polytopeFromVertices(pobjects[i].getVertexRepresentation(),(0,0,0,0.05))

###############################################################################
# Project Objects in footBox down, create clipped surfaces
###############################################################################

Wsurfaces_candidates= []
for i in range(0,len(wsurfaces)):

        ap = wsurfaces[i].ap
        bp = wsurfaces[i].bp
        iObject = wsurfaces[i].iObject
        p = ProjectPolytopesDownInsideBox(pobjects, wsurfaces[i], footBox[i])

        ###############################################################################
        # Create smaller walkable surfaces from the clipped surfaces
        ###############################################################################

        for j in range(0,len(p)):
                Wsplit = WalkableSurface.fromVertices(ap,bp,p[j],iObject)
                Wsurfaces_candidates.append(Wsplit)

###############################################################################
# Get all paths on top of walkable surfaces, and discard all other surfaces
###############################################################################
Wsurfaces_decomposed=[]
N_candidates = len(Wsurfaces_candidates)

[xstartI, xstartProj, xgoalI, xgoalProj] = getStartGoalWalkableSurfaces(Wsurfaces_candidates, xstart, xgoal)
#print "START:# ",xstart[0],"->",xstartProj[0]
#print "      # ",xstart[1],"->",xstartProj[1]
#print "      # ",xstart[2],"->",xstartProj[2]
#print "GOAL: # ",xgoal[0],"->",xgoalProj[0]
#print "      # ",xgoal[1],"->",xgoalProj[1]
#print "      # ",xgoal[2],"->",xgoalProj[2]

psize = 80
plot.point(xstartProj, size=psize,color=(1,0,0,0.9))
plot.point(xgoalProj,  size=psize,color=(1,0,0,0.9))
plot.point(xstart,     size=psize,color=(1,0,0,0.9))
plot.point(xgoal,      size=psize,color=(1,0,0,0.9))

WD = np.zeros((N_candidates,N_candidates))
WM = np.zeros((N_candidates,N_candidates))
for i in range(0,N_candidates):
        for j in range(i+1,N_candidates):
                WD[i,j]=WD[j,i]=distanceWalkableSurfaceWalkableSurface(Wsurfaces_candidates[i], Wsurfaces_candidates[j])
                WM[i,j]=WM[j,i]=(0 if WD[i,j]>0.1 else 1)

        WM[i,i]=1
print np.around(WD,2)
print np.around(WM,2)

G = nx.Graph()
for i in range(0,N_candidates):
        for j in range(0,N_candidates):
                if not WD[i,j]>0.01:
                        G.add_edge(i,j)

paths=list(nx.all_simple_paths(G, source=xstartI, target=xgoalI))

indices = []
print "----------------------------------------------------------------"
print "found ",len(paths)," paths between walkable surfaces:"
for i in range(0,len(paths)):
        print paths[i]
        for k in range(0,len(paths[i])):
                indices.append(paths[i][k])
print "----------------------------------------------------------------"

indices = sorted(set(indices))
print indices
for I in range(0,len(indices)):
        Wsurfaces_decomposed.append(Wsurfaces_candidates[indices[I]])

###############################################################################
# Visualize decomposed walkable surfaces for the foot
###############################################################################
Wsurface_box_compl = []
for i in range(0,len(Wsurfaces_decomposed)):
        plot.walkableSurface( \
                        Wsurfaces_decomposed[i].getVertexRepresentation(),\
                        fcolor=(0,1,0,0.2))
        Wsurface_box_compl.append(Wsurfaces_decomposed[i].createBox(0, \
                2*ROBOT_VOLUME_MAX_HEIGHT))


###############################################################################
# For each walkable surface part, inspect the head topology
###############################################################################
head_surfaces=[]
for i in range(0,len(Wsurfaces_decomposed)):
                WsplitBox = Wsurfaces_decomposed[i].createBox( \
                                ROBOT_VOLUME_MIN_HEIGHT-ROBOT_HEAD_SIZE,  \
                                ROBOT_VOLUME_MIN_HEIGHT, \
                                DeltaSide=ROBOT_MAX_HEAD_DISPLACEMENT)

                #plot.polytopeFromVertices( \
                #                WsplitBox.getVertexRepresentation(), \
                #                fcolor=(0,1,0,0.6))

                p2 = ProjectPolytopesDownInsideBox(pobjects, wsurfaces[i], WsplitBox)

                for k in range(0,len(p2)):
                        HeadSplit = WalkableSurface.fromVertices(ap,bp,p2[k],iObject)

                        dwp = distanceWalkableSurfacePolytope(HeadSplit,\
                                        Wsurface_box_compl[i].A,\
                                        Wsurface_box_compl[i].b)

                        if dwp <= 0.001:
                                head_surfaces.append(HeadSplit)

###############################################################################
# Visualize head surfaces
###############################################################################

for i in range(0,len(head_surfaces)):
        plot.walkableSurface(head_surfaces[i].getVertexRepresentation(),0.1,(0,0,1,0.6))

###############################################################################
# Find topology of traversable surfaces
###############################################################################
#N = len(Wsurfaces_decomposed)
#WD = np.zeros((N,N))
#WM = np.zeros((N,N))
#for i in range(0,N):
#        for j in range(i+1,N):
#                WD[i,j]=WD[j,i]=distanceWalkableSurfaceWalkableSurface(Wsurfaces_decomposed[i], Wsurfaces_decomposed[j])
#                WM[i,j]=WM[j,i]=(0 if WD[i,j]>0.1 else 1)
#
#        WM[i,i]=1
#
#print np.around(WD,2)
#print np.around(WM,2)
#
#optimizePath(xstartProj, xgoalProj, path, Wsurfaces_decomposed)

plot.showEnvironment()
