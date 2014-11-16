import numpy as np
import networkx as nx
from networkx import graphviz_layout
import pickle

from src.polytope import Polytope
from src.walkable import WalkableSurface, WalkableSurfacesFromPolytopes
from src.walkable import ProjectPolytopesDownInsideBox
from src.walkable import getStartGoalWalkableSurfaces
from src.urdfparser import URDFtoPolytopes
from src.linalg import intersection
from src.linalg import distanceWalkableSurfaceWalkableSurface
from src.linalg import distancePointPolytope
from src.linalg import distancePolytopePolytope
from src.linalg import distanceWalkableSurfacePolytope
from src.linalg import projectPointOntoHyperplane
from scipy.spatial import ConvexHull
from src.pathoptimizer import optimizePath
from plotter import Plotter

import sys

###############################################################################
# CONFIGURE / PARAMETERS 
###############################################################################
from src.robotspecifications import ROBOT_FOOT_RADIUS 
from src.robotspecifications import ROBOT_MAX_SLOPE
from src.robotspecifications import ROBOT_SPHERE_RADIUS
from src.robotspecifications import ROBOT_FOOT_HEIGHT
from src.robotspecifications import ROBOT_VOLUME_MIN_HEIGHT
from src.robotspecifications import ROBOT_VOLUME_MAX_HEIGHT
from src.robotspecifications import ROBOT_MAX_HEAD_DISPLACEMENT
from src.robotspecifications import ROBOT_HEAD_SIZE
from src.robotspecifications import ROBOT_MAX_UPPER_BODY_DISTANCE_FROM_FOOT
## minimal distance between walkable surfaces, such that we consider them
## connected
MIN_DISTANCE_WALKABLE_SURFACES=0.01

## START AND GOAL CONTACTS OF THE ROBOT
xstart = np.array((-0.5,2,0))
xgoal = np.array((0,-2,0))
env_fname = "wall.urdf"

## COLORS
colorScene=(0.3,0.3,0.3,0.2)
colorBodyBox=(1,0,0,0.2)
colorWalkableSurface=(1,0,1,0.8)
###############################################################################

plot = Plotter()

pobjects = URDFtoPolytopes(env_fname)
print "----------------------------------------------------------------"
print "Loaded environment",env_fname
print " -- Objects:",len(pobjects)

wsurfaces = WalkableSurfacesFromPolytopes(pobjects)

###############################################################################
# create footBox, box over Sip with height of one foot
###############################################################################

footBoxCandidate=[]
for i in range(0,len(wsurfaces)):
        footBoxCandidate.append(wsurfaces[i].createBox(0.01, ROBOT_FOOT_HEIGHT))

###############################################################################
# Visualize Complete Scene
###############################################################################

#for i in range(0,len(pobjects)):
#        V = pobjects[i].getVertexRepresentation()
#
#        middle = True
#        for j in range(0,len(V)):
#                if abs(V[j][1])>0.3:
#                        middle=False
#
#        if middle:
#                plot.polytopeFromVertices(\
#                                pobjects[i].getVertexRepresentation(),\
#                                fcolor=colorScene)

###############################################################################
# Project Objects in footBox down, create clipped surfaces
###############################################################################

Wsurfaces_candidates= []
for i in range(0,len(wsurfaces)):

        ap = wsurfaces[i].ap
        bp = wsurfaces[i].bp
        iObject = wsurfaces[i].iObject
        p = ProjectPolytopesDownInsideBox(pobjects, wsurfaces[i],\
                        footBoxCandidate[i])

        #######################################################################
        # Create smaller walkable surfaces from the clipped surfaces
        #######################################################################

        for j in range(0,len(p)):
                Wsplit = WalkableSurface.fromVertices(ap,bp,p[j],iObject)
                Wsurfaces_candidates.append(Wsplit)

print "splitted",len(wsurfaces),"walkable surfaces into",\
                len(Wsurfaces_candidates),"(reasoning about foot placement)"
print "----------------------------------------------------------------"


###############################################################################
# Get all paths on top of walkable surfaces, and discard all other surfaces
###############################################################################

Wsurfaces_decomposed=[]
N_candidates = len(Wsurfaces_candidates)

[xstartI, xstartProj, xgoalI, xgoalProj] = \
                getStartGoalWalkableSurfaces(Wsurfaces_candidates, xstart, xgoal)

psize = 80
plot.point(xstartProj, size=psize,color=(1,0,0,0.9))
plot.point(xgoalProj,  size=psize,color=(1,0,0,0.9))
plot.point(xstart,     size=psize,color=(1,0,0,0.9))
plot.point(xgoal,      size=psize,color=(1,0,0,0.9))

WD = np.zeros((N_candidates,N_candidates))
WM = np.zeros((N_candidates,N_candidates))
for i in range(0,N_candidates):
        for j in range(i+1,N_candidates):
                WD[i,j]=WD[j,i]=distanceWalkableSurfaceWalkableSurface(\
                                Wsurfaces_candidates[i], \
                                Wsurfaces_candidates[j])
                WM[i,j]=WM[j,i]=(0 if WD[i,j]>MIN_DISTANCE_WALKABLE_SURFACES else 1)

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
print "found",len(paths),"paths between walkable surfaces:"

for i in range(0,len(paths)):
        print paths[i]
        for k in range(0,len(paths[i])):
                indices.append(paths[i][k])

indices = sorted(set(indices))

print "----------------------------------------------------------------"
print "relevant walkable surfaces:", indices
print "----------------------------------------------------------------"

for I in range(0,len(indices)):
        Wsurfaces_decomposed.append(Wsurfaces_candidates[indices[I]])

###############################################################################
# compute paths inside of subgraph to get indices right
###############################################################################
N_w = len(Wsurfaces_decomposed)
WD = np.zeros((N_w,N_w))
WM = np.zeros((N_w,N_w))

[xstartI, xstartProj, xgoalI, xgoalProj] = \
                getStartGoalWalkableSurfaces(Wsurfaces_decomposed, xstart, xgoal)

for i in range(0,N_w):
        for j in range(i+1,N_w):
                WD[i,j]=WD[j,i]=distanceWalkableSurfaceWalkableSurface(\
                                Wsurfaces_decomposed[i], \
                                Wsurfaces_decomposed[j])
                WM[i,j]=WM[j,i]=(0 if WD[i,j]>MIN_DISTANCE_WALKABLE_SURFACES else 1)

        WM[i,i]=1
print np.around(WD,2)
print np.around(WM,2)

G = nx.Graph()
for i in range(0,N_w):
        for j in range(0,N_w):
                if not WD[i,j]>0.01:
                        G.add_edge(i,j)

paths=list(nx.all_simple_paths(G, source=xstartI, target=xgoalI))

print "----------------------------------------------------------------"
print "found",len(paths),"paths between walkable surfaces:"

for i in range(0,len(paths)):
        print paths[i]

## compute neighbors of each walkable surface
path_neighbors = []
for i in range(0,N_w):
        path_neighbors.append([])

for i in range(0,len(paths)):
        for k in range(0,len(paths[i])):
                curNode = paths[i][k]
                if k==0:
                        nextNode = paths[i][k+1]
                        path_neighbors[curNode].append(nextNode)
                else:
                        if k<len(paths[i])-1:
                                lastNode = paths[i][k-1]
                                nextNode = paths[i][k+1]
                                path_neighbors[curNode].append(nextNode)
                                path_neighbors[curNode].append(lastNode)
                                #path_neighbors[curNode] = [lastNode,nextNode]
                        else:
                                lastNode = paths[i][k-1]
                                path_neighbors[curNode].append(lastNode)

print "----------------------------------------------------------------"
print "Neighbor computation"

for i in range(0,len(path_neighbors)):
        path_neighbors[i] = sorted(set(path_neighbors[i]))
        curNode = i
        nNodes = path_neighbors[i]
        print "WS",curNode," has neighbors",nNodes



###############################################################################
# Visualize decomposed walkable surfaces for the foot
###############################################################################

for i in range(0,N_w):
        plot.walkableSurface( \
                        Wsurfaces_decomposed[i].getVertexRepresentation(),\
                        fcolor=colorWalkableSurface)

###############################################################################
# create complete box over walkable surface, in which swept volume has to lie
###############################################################################

Wsurface_box = []
Wsurface_box_extension = []
for i in range(0,N_w):
        box = Wsurfaces_decomposed[i].createBox(0.01, \
                        ROBOT_VOLUME_MAX_HEIGHT)
        box_ext = Wsurfaces_decomposed[i].createBox(0.01, \
                        ROBOT_VOLUME_MAX_HEIGHT,\
                        DeltaSide=ROBOT_MAX_UPPER_BODY_DISTANCE_FROM_FOOT)
        Wsurface_box.append(box)
        Wsurface_box_extension.append(box_ext)

Wsurface_objects = []
for i in range(0,N_w):
        objects_in_box_i = []
        for j in range(0,len(pobjects)):
                O = pobjects[j]
                B = Wsurface_box_extension[i]
                d = distancePolytopePolytope(O,B)
                if d < 0.0001:
                        objects_in_box_i.append(O)
        Wsurface_objects.append(objects_in_box_i)
print "----------------------------------------------------------------"
print "Objects on top of walkable surfaces:"
print "----------------------------------------------------------------"
for i in range(0,N_w):
        print "WS",i,"has",len(Wsurface_objects[i]),"objects associated"
print "----------------------------------------------------------------"

###############################################################################
# Compute stack of boxes on top of each walkable surface
# each stack will act as a constraint on the swept volume path optimization
###############################################################################

Wsurface_box_vstack = []
delta = ROBOT_VOLUME_MAX_HEIGHT/16

#for i in range(2,N_w-1):
for i in range(0,N_w):
        bottomHeight = ROBOT_FOOT_HEIGHT
        W = Wsurfaces_decomposed[i]
        objs = Wsurface_objects[i]
        for j in range(0,len(path_neighbors[i])):
                k = path_neighbors[i][j]
                print "removing",k,"from",i
                box = Wsurface_box[k]
                objs.append(box)
        #objs.append(box)
        #objs.append(box)

        #for j in range(0,len(objs)):
        #        plot.polytopeFromVertices(\
        #                objs[j].getVertexRepresentation(),\
        #                fcolor=colorScene)

        Wi_box_vstack=[]
        foot_box = W.createBox(0,bottomHeight)
        foot_stack = []
        foot_stack.append(foot_box)
        Wi_box_vstack.append(foot_stack)

        while bottomHeight<ROBOT_VOLUME_MAX_HEIGHT:
                ap = W.ap
                bp = W.bp+bottomHeight

                iObject=W.iObject

                delta_box = W.createBox(\
                                bottomHeight,\
                                bottomHeight+delta,\
                                DeltaSide=ROBOT_MAX_UPPER_BODY_DISTANCE_FROM_FOOT)

                pprime = ProjectPolytopesDownInsideBox(\
                                objs,\
                                W, \
                                delta_box)

                hstack = []
                for k in range(0,len(pprime)):

                        Ksplit = WalkableSurface.fromVertices(\
                                        ap,bp,pprime[k],iObject)

                        ##new polytope can be outside of walkable surface box
                        ## in which case there is no connection, and we have to
                        ## discard it
                        dwp = distanceWalkableSurfacePolytope(Ksplit,\
                                        Wsurface_box[i].A,\
                                        Wsurface_box[i].b-0.01)

                        if dwp <= 0.001:
                                ### TODO: check if object is changing topology of the box
                                box = Ksplit.createBox(0,delta)
                                hstack.append(box)

                if not len(hstack)>0:
                        ## layer is not existent, i.e. there is an object which
                        ## blocks the layer completely, for example a roof
                        ## => break and goto next walkable surface
                        break

                Wi_box_vstack.append(hstack)
                bottomHeight+=delta

        Wsurface_box_vstack.append(Wi_box_vstack)
###############################################################################
# Store xstart, xgoal, paths, wsurfaces and the vstack
###############################################################################

pickle.dump( xstartProj, open( "data/xstart.dat", "wb" ) )
pickle.dump( xgoalProj, open( "data/xgoal.dat", "wb" ) )
pickle.dump( paths, open( "data/paths.dat", "wb" ) )
pickle.dump( Wsurfaces_decomposed, open( "data/wsurfaces.dat", "wb" ) )
pickle.dump( Wsurface_box_vstack, open( "data/wsurfaces_vstack.dat", "wb" ) )

xstartProj = pickle.load( open( "data/xstart.dat", "rb" ) )
xgoalProj = pickle.load( open( "data/xgoal.dat", "rb" ) )
paths = pickle.load( open( "data/paths.dat", "rb" ) )
Wsurfaces_decomposed = pickle.load( open( "data/wsurfaces.dat", "rb" ) )
Wsurface_box_vstack = pickle.load( open( "data/wsurfaces_vstack.dat", "rb" ) )

###############################################################################
### print summary of vstacks on each surface
###############################################################################
for i in range(0,len(Wsurface_box_vstack)):
        ### TODO: remove, just for visualizing
        if i==3:
                continue
        vstack = Wsurface_box_vstack[i]
        print "WS",i,"has",len(vstack),"layers"
        for j in range(0,len(vstack)):
                hstack = vstack[j]
                print "  layer",j,"is decomposed into",len(hstack),"boxes"
                for k in range(0,len(hstack)):
                        plot.polytopeFromVertices(\
                                hstack[k].getVertexRepresentation(),\
                                fcolor=colorBodyBox)

plot.showEnvironment()
