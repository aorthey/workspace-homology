import numpy as np
import cvxpy as cvx
import pickle
import sys
sys.path.append("..")
from plotter import Plotter
from numpy import inf,array,zeros
from cvxpy import *
from math import tan,pi
from src.util import *

from src.linalg import *
from src.walkable import WalkableSurface, WalkableSurfacesFromPolytopes
from src.robotspecifications import *

start= pickle.load( open( "data/xstart.dat", "rb" ) )
goal= pickle.load( open( "data/xgoal.dat", "rb" ) )
path_candidates = pickle.load( open( "data/paths.dat", "rb" ) )
Wsurfaces_decomposed = pickle.load( open( "data/wsurfaces.dat", "rb" ) )
Wsurface_box_vstack = pickle.load( open( "data/wsurfaces_vstack.dat", "rb" ) )

HeadName = "data/xspacemanifold-same-axes/headersamples.dat"
[Npts, VSTACK_DELTA, heights] = pickle.load( open( HeadName, "rb" ) )

Aname = "data/polytopes/A.dat"
ARKname = "data/polytopes/Ark.dat"
bname = "data/polytopes/b.dat"

Aflat = pickle.load( open( Aname, "rb" ) )
Aleftrightconv = pickle.load( open( ARKname, "rb" ) )
bflat = pickle.load( open( bname, "rb" ) )

XspaceDimension = Aflat[0].shape[1]
XspaceMinima = len(Aflat)
print XspaceDimension
###############################################################################
# Check which path has the best chances of finding a feasible trajectory
###############################################################################
path = path_candidates[1]
N_walkablesurfaces = len(path)
print "path over ",path,"walkable surfaces"

###############################################################################
# parameters
###############################################################################

A = np.eye(3);
A[2,2]=0;

cknee = np.array((0,0,tan(ROBOT_APERTURE_KNEE_FOOT)))
chipfoot = np.array((0,0,tan(ROBOT_APERTURE_HIP_FOOT)))
chip= np.array((0,0,tan(ROBOT_APERTURE_HIP_KNEE)))
cwaist= np.array((0,0,tan(ROBOT_APERTURE_WAIST_HIP)))
cneck = np.array((0,0,tan(ROBOT_APERTURE_NECK_WAIST)))
chead = np.array((0,0,tan(ROBOT_APERTURE_HEAD_NECK)))

###############################################################################
# optimize over the connection between surfaces
###############################################################################
plot = Plotter()

## compute connectors
connector = []
upperBodyConnector = []

for i in range(0,N_walkablesurfaces-1):
        pcurrent = path[i]
        pnext = path[i+1]
        print "forming connection between ",pcurrent,"and",pnext
        W=Wsurfaces_decomposed[pcurrent]
        Wnext=Wsurfaces_decomposed[pnext]

        Wstack = Wsurface_box_vstack[pcurrent]
        WstackNext = Wsurface_box_vstack[pnext]

        #### create foot connector
        boxW = W.createBox(0,ROBOT_FOOT_HEIGHT,DeltaSide=0.02)
        boxWnext = Wnext.createBox(0,ROBOT_FOOT_HEIGHT)
        binter = boxW.intersectWithPolytope(Wnext)
        ap = Wnext.ap
        bp = Wnext.bp
        connectorWtoWnext = WalkableSurface(ap,bp,binter.A, binter.b, Wnext.iObject)
        plot.walkableSurface( \
                        connectorWtoWnext.getVertexRepresentation(),\
                        fcolor=[1,0,0,0.4],\
                        thickness=0.01)
        connector.append(connectorWtoWnext)

        ### create upper body connectors
        colorScene=(0.8,0.1,0.1,0.5)
        upperBodyConnectorStack = []
        for j in range(0,XspaceDimension):
                if j < len(Wstack):
                        Winter = Wstack[j][0].intersectWithPolytope(WstackNext[j][0])
                        upperBodyConnectorStack.append(Winter)
                        plot.polytopeFromVertices(\
                                        Winter.getVertexRepresentation(),\
                                        fcolor=colorScene)
                else:
                        ## create helper box for higher dimensions
                        Whelper = connectorWtoWnext.createTinyHelperBox(heights[j],heights[j]+VSTACK_DELTA)
                        upperBodyConnectorStack.append(Whelper)
                        plot.polytopeFromVertices(\
                                        Whelper.getVertexRepresentation(),\
                                        fcolor=colorScene)

        upperBodyConnector.append(upperBodyConnectorStack)
###############################################################################
# building variables
###############################################################################
## the x positions on the connection between two walkable surfaces
x_goal = Variable(3,1)
x_start = Variable(3,1)

x_connection = []
x_connection_upper_body = []
for i in range(0,N_walkablesurfaces-1):
        x_connection.append(Variable(3,1))
        yycon = []
        for j in range(0,Npts):
                yycon.append(Variable(3,1))
        x_connection_upper_body.append(yycon)

## compute number of points per walkable surface, depending on the distance to
## the next connector segment
M_w = []

N_c = len(connector)
dI = distancePointWalkableSurface(start, connector[0])
M_w.append(distanceInMetersToNumberOfPoints(dI))

for i in range(0,N_c-1):
        C = connector[i]
        Cnext = connector[i+1]
        dcc = distanceWalkableSurfaceWalkableSurface(C,Cnext)
        M_w.append(distanceInMetersToNumberOfPoints(dcc))

dG = distancePointWalkableSurface(goal, connector[N_c-1])
M_w.append(distanceInMetersToNumberOfPoints(dG))

print M_w

x_WS = []
for i in range(0,N_walkablesurfaces):
        print "WS ",i,"has",M_w[i],"points to optimize"
        x_WS.append(Variable(3,M_w[i]))


###############################################################################
# building constraints
###############################################################################

constraints = []

## start/goal regions constraints
constraints.append(norm(x_start - start) <= PATH_RADIUS_START_REGION)
constraints.append(norm(x_goal - goal) <= PATH_RADIUS_GOAL_REGION)

for i in range(0,N_walkablesurfaces-1):
        ### constraints for points on the connection of the WS
        C = connector[i]
        y = x_connection[i]
        constraints.append( np.matrix(C.A)*y <= C.b )
        constraints.append( np.matrix(C.ap)*y == C.bp)

        ### constraints for upper body constraints at the intersection pts
        #yy = x_connection_upper_body[i]
        #Cstack = upperBodyConnector[i]
        #for j in range(0,XspaceDimension):
        #        constraints.append( np.matrix(Cstack[j].A)*yy[j] == Cstack[j].b)

###############################################################################
# building objective
###############################################################################
v2 = np.array((0,0,1))
v1 = np.array((1,0,0))
v3 = np.array((0,1,0))
for i in range(0,XspaceMinima):
        objective = Minimize(norm(x_connection[0]))
        Ae = Aflat[i]
        be = bflat[i]
        mincon = []
        Ai = Aleftrightconv[i]

        rho = Variable(XspaceDimension)
        mincon.append( np.matrix(Ae)*rho <= be)
        mincon.append( np.matrix(Ae)*np.matrix(Ai)*rho <= be)

        for j in range(0,N_walkablesurfaces-1):
                Cstack = upperBodyConnector[j]
                for k in range(0,Npts):
                        vv = x_connection[j] - rho[k]*v1 + heights[k]*v2
                        mincon.append( np.matrix(Cstack[k].A)*vv <= Cstack[k].b)

        for j in range(0,len(constraints)):
                mincon.append(constraints[j])


        ###############################################################################
        # solve
        ###############################################################################
        #print "optimizing path from",start,"to",goal
        prob = Problem(objective, mincon)
        #prob.solve(solver=cvx.SCS, verbose=True)
        prob.solve(solver=cvx.SCS)
        #prob.solve(solver=cvx.CVXOPT)
        print prob.value
        if prob.value < inf:
                print "minima",i,"admits a solution"
                break

###############################################################################
# plot
###############################################################################
if prob.value < inf:
        plot.point(x_goal.value,color=(0,0,0,0.9))
        plot.point(x_start.value)
        for i in range(0,N_walkablesurfaces-1):
                if x_connection[i].value is not None:
                        plot.point(x_connection[i].value)
                        for k in range(0,Npts):
                                offset = 0.2
                                pt = x_connection[i].value+(-rho.value[k]*v1+heights[k]*v2+offset*v3).T
                                plot.point(pt,size=100,color=(0,0,0,1))

                        rhoR = np.matrix(Ai)*rho.value

                        for k in range(0,Npts):
                                offset = 0.2
                                pt = x_connection[i].value+(-rhoR[k]*v1+heights[k]*v2+offset*v3).T
                                plot.point(pt,size=100,color=(0,0,0,1))

plot.set_view(90,0)
plot.showEnvironment()
