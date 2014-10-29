import numpy as np
import networkx as nx
from networkx import graphviz_layout

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
## minimal distance between walkable surfaces, such that we consider them
## connected
MIN_DISTANCE_WALKABLE_SURFACES=0.01

## START AND GOAL CONTACTS OF THE ROBOT
xstart = np.array((-0.5,2,0))
xgoal = np.array((0,-2,0))
env_fname = "wall.urdf"

## COLORS
colorScene=(0.3,0.3,0.3,0.1)
colorBodyBox=(1,0,0,0.1)
colorFootBox=(1,1,0,0.1)
colorHeadBox=(0.2,0.2,0.2,0.1)
colorWalkableSurface=(1,0,1,0.8)
###############################################################################

plot = Plotter()

pobjects = URDFtoPolytopes(env_fname)
wsurfaces = WalkableSurfacesFromPolytopes(pobjects)

###############################################################################
# create footBox, box over Sip with height of one foot
###############################################################################

footBoxCandidate=[]
for i in range(0,len(wsurfaces)):
        footBoxCandidate.append(wsurfaces[i].createBox(0, ROBOT_FOOT_HEIGHT))

###############################################################################
# Visualize Complete Scene
###############################################################################

for i in range(0,len(pobjects)):
        plot.polytopeFromVertices(\
                        pobjects[i].getVertexRepresentation(),\
                        fcolor=colorScene)

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
# Visualize decomposed walkable surfaces for the foot
###############################################################################

Wsurface_box_compl = []
for i in range(0,len(Wsurfaces_decomposed)):
        plot.walkableSurface( \
                        Wsurfaces_decomposed[i].getVertexRepresentation(),\
                        fcolor=colorWalkableSurface)
        Wsurface_box_compl.append(Wsurfaces_decomposed[i].createBox(0, \
                2*ROBOT_VOLUME_MAX_HEIGHT))

###############################################################################
# Visualize foot boxes
###############################################################################
foot_boxes = []
for i in range(0,len(Wsurfaces_decomposed)):
        box = Wsurfaces_decomposed[i].createBox(0, ROBOT_FOOT_HEIGHT)
        foot_boxes.append( box )
        plot.polytopeFromVertices(box.getVertexRepresentation(),\
                        fcolor=colorFootBox)

###############################################################################
# For each walkable surface part, inspect the head topology
###############################################################################

head_boxes=[]
for i in range(0,len(Wsurfaces_decomposed)):
        W = Wsurfaces_decomposed[i]
        ap = W.ap
        bp = W.bp+ROBOT_VOLUME_MIN_HEIGHT-ROBOT_HEAD_SIZE

        head_boxes_i=[]
        WsplitBox = Wsurfaces_decomposed[i].createBox( \
                        ROBOT_VOLUME_MIN_HEIGHT-ROBOT_HEAD_SIZE,  \
                        ROBOT_VOLUME_MIN_HEIGHT, \
                        DeltaSide=ROBOT_MAX_HEAD_DISPLACEMENT)

        p2 = ProjectPolytopesDownInsideBox(\
                        pobjects,\
                        Wsurfaces_decomposed[i], \
                        WsplitBox)

        for k in range(0,len(p2)):
                HeadSplit = WalkableSurface.fromVertices(\
                                ap,bp,p2[k],iObject)

                dwp = distanceWalkableSurfacePolytope(HeadSplit,\
                                Wsurface_box_compl[i].A,\
                                Wsurface_box_compl[i].b)

                if dwp <= 0.001:
                        HeadSplitBox = HeadSplit.createBox(0,ROBOT_HEAD_SIZE)
                        head_boxes_i.append(HeadSplitBox)

        head_boxes.append(head_boxes_i)

###############################################################################
# Visualize head surfaces
###############################################################################

print "----------------------------------------------------------------"
for i in range(0,len(head_boxes)):
        print "walkable surface",indices[i],"has",\
                        len(head_boxes[i]),"head homotopy classes"
print "----------------------------------------------------------------"

for i in range(0,len(head_boxes)):
        for j in range(0,len(head_boxes[i])):
                plot.walkableSurface( \
                                head_boxes[i][j].getVertexRepresentation(),\
                                0.1,fcolor=colorHeadBox)

###############################################################################
# Analyse body box to obtain the body homotopy classes
###############################################################################

## get body box by projecting head surface down onto walkable surface, then
## taking the convex hull, and building up a new box on top of this convex hull

        #############################################################################
        ### Creating body boxes
        #############################################################################
body_boxes = []
for i in range(0,len(head_boxes)):
        W = Wsurfaces_decomposed[i]
        ap = W.ap
        bp = W.bp
        Vw = W.getVertexRepresentation()

        for j in range(0,len(head_boxes[i])):
                V = head_boxes[i][j].getVertexRepresentation()
                Vp = []
                for k in range(0,len(V)):
                        vproj = projectPointOntoHyperplane(V[k],ap,bp)
                        Vp.append(vproj)

                Vunion_all=np.vstack((Vw,Vp))

                hull = ConvexHull(Vunion_all[:,0:2])
                Vunion=Vunion_all[hull.vertices,:]
                w_union = WalkableSurface.fromVertices(ap,bp,Vunion,W.iObject)
                #plot.walkableSurface(w_union.getVertexRepresentation(),\
                                #thickness=0.2,fcolor=(1,1,0,0.8))
                body_boxes.append(\
                                w_union.createBox(
                                        ROBOT_FOOT_HEIGHT,\
                                        ROBOT_VOLUME_MIN_HEIGHT-\
                                        ROBOT_HEAD_SIZE))  \
        #############################################################################
        # Visualize body box
        #############################################################################

for i in range(0,len(body_boxes)):
        plot.polytopeFromVertices( \
                        body_boxes[i].getVertexRepresentation(), \
                        fcolor=colorBodyBox)

        #############################################################################
        # Obtain list of object inside body box (used for planning)
        #############################################################################
print "----------------------------------------"
for i in range(0,len(body_boxes)):
        B = body_boxes[i]
        B_obj = []
        for j in range(0,len(pobjects)):
                O = pobjects[j]
                d = distancePolytopePolytope(O, B)
                if d < 0.001:
                        Ointer = O.intersectWithPolytope(B)
                        B_obj.append(Ointer)
                        plot.polytopeFromVertices(\
                                        Ointer.getVertexRepresentation(),\
                                        fcolor=(1,0,0,0.6))
        print "Body box ",i,"has",len(B_obj),"objects associated"

        ### compute distance between objects inside B
        D_obj = np.zeros((len(B_obj),len(B_obj)))
        for j in range(0,len(B_obj)):
                for k in range(0,len(B_obj)):
                        Bj = B_obj[j]
                        Bk = B_obj[k]
                        D_obj[j,k]=D_obj[k,j]=distancePolytopePolytope(Bj,Bk)
                D_obj[j,j]=1000

        print np.around(D_obj,3)

        homotopyChangingObjects = 0
        for j in range(0,len(B_obj)):
                md = np.nanmin(D_obj[j,:])
                print D_obj[j,:]
                if md > ROBOT_SPHERE_RADIUS:
                        ## object changes topology!
                        B_obj[i].changesTopology(True)
                        homotopyChangingObjects += 1
                        print " object",j,"changes topology of body box!"

        print " => ",homotopyChangingObjects,"topological holes"
        print "----------------------------------------"


###############################################################################
# Optimization of paths in each homotopy class
###############################################################################
for k in range(0,len(paths)):
        optimizePath(xstartProj, xgoalProj, paths[k], Wsurfaces_decomposed)

plot.showEnvironment()
