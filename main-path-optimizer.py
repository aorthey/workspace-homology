import numpy as np
import cvxpy as cvx
import pickle
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

###############################################################################
# Check which path has the best chances of finding a feasible trajectory
###############################################################################
path = path_candidates[1]
N_w = len(path)
print "path over ",path,"walkable surfaces"

###############################################################################
# parameters
###############################################################################
startRegionRadius = 0.05
goalRegionRadius = 0.05
distanceWaypoints = 0.5

robotMaxHeight = 1.539
robotMinKneeHeight = 0.2
robotMaxKneeHeight = 0.3
robotKneeAperture = pi/8
robotHipAperture = pi/200
robotNeckAperture = pi/16
robotHeadAperture = pi/16
robotKneeHeadAperture = pi/32
robotRadiusSweptVolume = 0.1

A = np.eye(3);
A[2,2]=0;
cknee = np.array((0,0,tan(robotKneeAperture)))
chip = np.array((0,0,tan(robotHipAperture)))
cneck = np.array((0,0,tan(robotNeckAperture)))
chead = np.array((0,0,tan(robotHeadAperture)))
ckneehead = np.array((0,0,tan(robotKneeHeadAperture)))

###############################################################################
# optimize over the connection between surfaces
###############################################################################
plot = Plotter()

## compute connectors
connector = []

for i in range(0,N_w-1):
        pi = path[i]
        pnext = path[i+1]
        print "forming connection between ",pi,"and",pnext
        W=Wsurfaces_decomposed[pi]
        Wnext=Wsurfaces_decomposed[pnext]

        ##create box, which intersects the next box by 2cm
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

###############################################################################
# building variables
###############################################################################
## the x positions on the connection between two walkable surfaces
x_goal = Variable(3,1)
x_start = Variable(3,1)

x_connection = []
for i in range(0,N_w-1):
        x_connection.append(Variable(3,1))


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
for i in range(0,N_w):
        print "WS ",i,"has",M_w[i],"points to optimize"
        x_WS.append(Variable(3,M_w[i]))

## for each foot contact point add the upper body points, one for each layer
###############################################################################
# building constraints
###############################################################################

constraints = []

## start/goal regions constraints
constraints.append(norm(x_goal - goal) <= goalRegionRadius)
constraints.append(norm(x_start - start) <= startRegionRadius)

for i in range(0,N_w-1):
        C = connector[i]
        y = x_connection[i]
        constraints.append( np.matrix(C.A)*y <= C.b )
        constraints.append( np.matrix(C.ap)*y == C.bp)

###############################################################################
# building objective
###############################################################################
objective = Minimize(norm(x_goal))

###############################################################################
# solve
###############################################################################
print "optimizing path from",start,"to",goal
prob = Problem(objective, constraints)
#prob.solve(solver=cvx.SCS, verbose=True)
prob.solve(solver=cvx.SCS)
#prob.solve(solver=cvx.CVXOPT)
print prob.value

###############################################################################
# plot
###############################################################################
if prob.value < inf:
        plot.point(x_goal.value)
        plot.point(x_start.value)
        for i in range(0,N_w-1):
                if x_connection[i].value is not None:
                        plot.point(x_connection[i].value)

plot.showEnvironment()
