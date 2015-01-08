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

DEBUG=1

svLeftColor = (0.5,0,0.5,1)
svRightColor = (0.5,0,0,1)
colorScene=(0.6,0.6,0.6,0.2)
svPointSize=50
### start/goal direction
start_normal = np.array((1,1,0))
goal_normal = np.array((1,0,0))

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
print "found",len(path_candidates),"paths"
for i in range(0,len(path_candidates)):
        print i,":",path_candidates[i]

print "choosing path",1,"over ",path,"walkable surfaces"

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
        upperBodyConnectorStack = []
        for j in range(0,XspaceDimension):
                if j < len(Wstack):
                        Winter = Wstack[j][0].intersectWithPolytope(WstackNext[j][0])
                        upperBodyConnectorStack.append(Winter)
                else:
                        ## create helper box for higher dimensions
                        Whelper = connectorWtoWnext.createTinyHelperBox(heights[j],heights[j]+VSTACK_DELTA)
                        upperBodyConnectorStack.append(Whelper)

        upperBodyConnector.append(upperBodyConnectorStack)
###############################################################################
# building variables
###############################################################################
## the x positions on the connection between two walkable surfaces
x_goal = Variable(3,1)
x_start = Variable(3,1)

x_connection = []
#x_connection_upper_body = []
for i in range(0,N_walkablesurfaces-1):
        x_connection.append(Variable(3,1))
        #yycon = []
        #for j in range(0,Npts):
        #        yycon.append(Variable(3,1))
        #x_connection_upper_body.append(yycon)

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

x_WS = []
for i in range(0,N_walkablesurfaces):
        print "WS ",i,"has",M_w[i],"points to optimize"

        x_WS_tmp = []
        for j in range(0,M_w[i]):
                x_WS_tmp.append(Variable(3,1))
        x_WS.append(x_WS_tmp)
        #x_WS.append(Variable(3,M_w[i]))


###############################################################################
# building constraints
###############################################################################

constraints = []

## start/goal regions constraints
constraints.append(norm(x_start - start) <= PATH_RADIUS_START_REGION)
constraints.append(norm(x_goal - goal) <= PATH_RADIUS_GOAL_REGION)

pstart = path[0]
pend = path[len(path)-1]
W=Wsurfaces_decomposed[pstart]
constraints.append( np.matrix(W.A)*x_start <= W.b )
constraints.append( np.matrix(W.ap)*x_start == W.bp)
W=Wsurfaces_decomposed[pend]
constraints.append( np.matrix(W.A)*x_goal <= W.b )
constraints.append( np.matrix(W.ap)*x_goal == W.bp)


for i in range(0,N_walkablesurfaces-1):
        ### constraints for points on the connection of the WS
        C = connector[i]
        y = x_connection[i]
        constraints.append( np.matrix(C.A)*y <= C.b )
        constraints.append( np.matrix(C.ap)*y == C.bp)

for i in range(0,N_walkablesurfaces):
        #### x_WS should contain only points on the surface
        W=Wsurfaces_decomposed[path[i]]
        plot.walkableSurface( \
                        W.getVertexRepresentation(),\
                        fcolor=colorScene,\
                        thickness=0.01)
        for j in range(0,M_w[i]):
                constraints.append( np.matrix(W.A)*x_WS[i][j] <= W.b)
                constraints.append( np.matrix(W.ap)*x_WS[i][j] == W.bp)

### constraint: x_WS should lie in the same functional space

### constraint: x_WS connections
constraints.append( x_WS[0][0] == x_start)
constraints.append( x_WS[0][M_w[0]-1] == x_connection[0])

for i in range(1,N_walkablesurfaces-1):
        Lws = len(x_WS[i])-1
        constraints.append( x_WS[i][0] == x_connection[i-1])
        constraints.append( x_WS[i][Lws] == x_connection[i])

Lws = len(x_WS[N_walkablesurfaces-1])-1
constraints.append( x_WS[N_walkablesurfaces-1][0] == x_connection[N_walkablesurfaces-2])
constraints.append( x_WS[N_walkablesurfaces-1][Lws] == x_goal)

### constraint: distance between x_WS
for i in range(0,N_walkablesurfaces):
        for j in range(1,len(x_WS[i])):
                constraints.append(norm(x_WS[i][j] - x_WS[i][j-1]) <= PATH_DIST_WAYPOINTS_MAX)

###############################################################################
# building objective
###############################################################################
v2 = np.array((0,0,1))
v1 = np.array((1,0,0))
v3 = np.array((0,1,0))
minimaIter = 0

if DEBUG:
        startMinima = 20
else:
        startMinima = 0

for i in range(startMinima,XspaceMinima):
        Ae = Aflat[i]
        be = bflat[i]
        mincon = []
        Ark = Aleftrightconv[i]
        rho = Variable(XspaceDimension)

        ## constraint: points only from manifold flat inside X
        objective = Minimize(norm(rho)+norm(np.matrix(Ark)*rho))
        mincon.append( np.matrix(Ae)*rho <= be)

        ## how many dimensions are there until the linear subspace starts?
        maxNonzeroDim = np.max(np.nonzero(Ark)[0])

        ## constraint: all intersection points have to be inside of an environment box

        for j in range(0,N_walkablesurfaces-1):
                Ebox = upperBodyConnector[j]
                for k in range(0,maxNonzeroDim):
                        vv = x_connection[j] + rho[k]*v1 + heights[k]*v2
                        mincon.append( np.matrix(Ebox[k].A)*vv <= Ebox[k].b)
                        rhoR = np.matrix(Ark)*rho
                        vvR = x_connection[j] + rhoR[k]*v1 + heights[k]*v2
                        mincon.append( np.matrix(Ebox[k].A)*vvR <= Ebox[k].b)

        for j in range(0,len(constraints)):
                mincon.append(constraints[j])

        ###############################################################################
        # solve
        ###############################################################################
        prob = Problem(objective, mincon)
        #prob.solve(solver=cvx.SCS, verbose=True)
        prob.solve(solver=cvx.SCS)
        #prob.solve(solver=cvx.CVXOPT)
        print prob.value,"(minima",i,"/",XspaceMinima,")"
        if prob.value < inf:
                ## BICONVEX condition: check second convex problem on feasibility
                ### constraint: all points have to be inside of an environment box
                ## inbetween Rho
                ibRho = []
                mincon = []
                objective=[]
                vp = []
                for j in range(0,N_walkablesurfaces):
                        #W=Wsurfaces_decomposed[path[j]]
                        W = Wsurface_box_vstack[path[j]]
                        ibRho_tmp=[]
                        vp_tmp=[]
                        for p in range(0,len(x_WS[j])):
                                ibRho_tmp.append(Variable(XspaceDimension))
                                ## vp: normal to tangent at x_WS[j][p]
                                vp_tmp.append(v1)

                        for p in range(0,len(x_WS[j])):
                                mincon.append( np.matrix(Ae)*ibRho_tmp[p] <= be)
                                for k in range(0,maxNonzeroDim):
                                        vv = x_WS[j][p].value + ibRho_tmp[p][k]*vp_tmp[p] + heights[k]*v2
                                        mincon.append( np.matrix(W[k][0].A)*vv <= W[k][0].b)
                                        rhoR = np.matrix(Ark)*ibRho_tmp[p]
                                        vvR = x_WS[j][p].value + rhoR[k]*vp_tmp[p] + heights[k]*v2
                                        mincon.append( np.matrix(W[k][0].A)*vvR <= W[k][0].b)

                        ibRho.append(ibRho_tmp)
                        vp.append(vp_tmp)

                objective = Minimize(norm(ibRho[0][0]))
                prob = Problem(objective, mincon)
                prob.solve(solver=cvx.SCS)
                print prob.value,"(minima",i,"/",XspaceMinima,")"

                if prob.value<inf:
                        minimaIter = i
                        print "minima",i,"admits a solution"
                        break

###############################################################################
# plot
###############################################################################
if prob.value < inf:
        print "plotting workspace swept volume approximation"
        Ark = Aleftrightconv[minimaIter]
        rhoR = np.matrix(Ark)*rho.value
        plot.point(x_goal.value,color=(0,0,0,0.9))
        plot.point(x_start.value)
        maxNonzeroDim = np.max(np.nonzero(Ark)[0])
        firstIS = None
        lastIS = None
        for i in range(0,N_walkablesurfaces-1):
                if x_connection[i].value is not None:

                        #### plot intersection foot pos
                        plot.point(x_connection[i].value)

                        #### plot intersection SV boundary left
                        offset = 0.15 ## for visibility reasons
                        for k in range(0,maxNonzeroDim):
                                pt = x_connection[i].value+(rho.value[k]*v1+heights[k]*v2+offset*v3).T
                                plot.point(pt,size=svPointSize,color=svLeftColor)

                        #### plot intersection SV boundary right
                        rhoR = np.matrix(Ark)*rho.value
                        for k in range(0,maxNonzeroDim):
                                pt = x_connection[i].value+(rhoR[k]*v1+heights[k]*v2+offset*v3).T 
                                plot.point(pt,size=svPointSize,color=svRightColor)

                        #### plot intersection E-boxes
                        for k in range(0,maxNonzeroDim):
                                W = upperBodyConnector[i][k]
                                plot.polytopeFromVertices(\
                                                W.getVertexRepresentation(),\
                                                fcolor=colorScene)

        ### plot goal/start swept volume boundary
        for k in range(0,maxNonzeroDim):
                pt = x_start.value+(rho.value[k]*start_normal+heights[k]*v2).T
                plot.point(pt,size=svPointSize,color=svLeftColor)
                pt = x_goal.value+(rho.value[k]*goal_normal+heights[k]*v2).T
                plot.point(pt,size=svPointSize,color=svLeftColor)
        for k in range(0,maxNonzeroDim):
                rhoR = np.matrix(Ark)*rho.value
                pt = x_start.value+(rhoR[k]*start_normal+heights[k]*v2).T
                plot.point(pt,size=svPointSize,color=svRightColor)
                pt = x_goal.value+(rhoR[k]*goal_normal+heights[k]*v2).T
                plot.point(pt,size=svPointSize,color=svRightColor)

        ### plot paths on each WS
        for i in range(0,N_walkablesurfaces):
                for j in range(0,len(x_WS[i])):
                        plot.point(x_WS[i][j].value,size=svPointSize,color=svLeftColor)
                        for k in range(0,maxNonzeroDim):
                                pt = x_WS[i][j].value.T+(ibRho[i][j][k].value*vp[i][j]+heights[k]*v2).T
                                plot.point(np.array(pt).flatten(),size=svPointSize,color=svLeftColor)
                        ibRhoR = np.matrix(Ark)*ibRho[i][j].value
                        for k in range(0,maxNonzeroDim):
                                pt = x_WS[i][j].value+(ibRhoR[k]*vp[i][j]+heights[k]*v2).T
                                plot.point(np.array(pt).flatten(),size=svPointSize,color=svRightColor)

        plot.set_view(90,0)
        plot.showEnvironment()
else:
        print "problem not feasible"
