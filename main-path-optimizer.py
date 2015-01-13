from timeit import default_timer as timer
import numpy as np
import cvxpy as cvx
import pickle
import sys
sys.path.append("..")
from plotter import Plotter,rotFromRPY
from numpy import inf,array,zeros
from cvxpy import *
from math import tan,pi
from src.util import *
from src.linalg import *
from src.walkable import WalkableSurface, WalkableSurfacesFromPolytopes
from xspace.htoq import *
from src.robotspecifications import * 
start= pickle.load( open( "data/xstart.dat", "rb" ) )
goal= pickle.load( open( "data/xgoal.dat", "rb" ) )
path_candidates = pickle.load( open( "data/paths.dat", "rb" ) )
Wsurfaces_decomposed = pickle.load( open( "data/wsurfaces.dat", "rb" ) )
Wsurface_box_vstack = pickle.load( open( "data/wsurfaces_vstack.dat", "rb" ) )

DEBUG=1
if DEBUG:
        goal[0]=0.4

v2 = np.array((0,0,1))
v1 = np.array((1,0,0))
v3 = np.array((0,1,0))

svLeftColor = (0.5,0,0.5,1)
svRightColor = (0.5,0,0,1)
colorScene=(0.6,0.6,0.6,0.2)
svPointSize=50

### start/goal direction
start_normal = np.array((1,1,0))
goal_normal = np.array((1,0,0))

startNormal = np.array((1,1,0))
goalNormal = np.array((1,0,0))
startNormalNormal = np.dot(rotFromRPY(0,0,-pi/2),startNormal)
goalNormalNormal = np.dot(rotFromRPY(0,0,-pi/2),goalNormal)


HeadName = "data/xspacemanifold-same-axes/headersamples.dat"
[Npts, VSTACK_DELTA, heights] = pickle.load( open( HeadName, "rb" ) )

Aname = "data/polytopes/A.dat"
ARKname = "data/polytopes/Ark.dat"
bname = "data/polytopes/b.dat"
HName = "data/polytopes/H.dat"

Harray = pickle.load( open( HName, "rb" ) )
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
for i in range(0,N_walkablesurfaces-1):
        x_connection.append(Variable(3,1))

###############################################################################
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


### we have to know the intersection normals, and the normal to those normals
connectorNormalNormal = []
connectorNormal= []
for i in range(0,N_walkablesurfaces-1):
        ## TODO: use normal of intersection, and rotate it via
        ## rotfromrpy(0,0,pi) to the right 
        connectorNormalNormal.append(np.array((1,0,0)))
        connectorNormal.append(np.array((0,1,0)))

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

###############################################################################
for i in range(0,N_walkablesurfaces-1):
        ### constraints for points on the connection of the WS
        C = connector[i]
        y = x_connection[i]
        constraints.append( np.matrix(C.A)*y <= C.b )
        constraints.append( np.matrix(C.ap)*y == C.bp)

###############################################################################
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

###############################################################################
### constraint: x_WS should lie in the same functional space

###############################################################################
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

###############################################################################
### constraint: distance between x_WS
for i in range(0,N_walkablesurfaces):
        for j in range(1,len(x_WS[i])):
                constraints.append(norm(x_WS[i][j] - x_WS[i][j-1]) <= PATH_DIST_WAYPOINTS_MAX)

###############################################################################
### constraint: x_WS before the intersection should have the same orientation as
### the intersection

for i in range(0,N_walkablesurfaces-1):
        v = connectorNormal[i]
        Lws = len(x_WS[i])-1
        ## trois points: x before the intersection, x at the intersection and x
        ## after the intersection. They all should lie on a line perpendicular to
        ## the normalnormal of the intersection
        xbefore = x_WS[i][Lws-1]
        xconnect = x_WS[i+1][0]
        xafter = x_WS[i+1][1]

        gammaA = Variable(1)
        gammaB = Variable(1)
        constraints.append( gammaA*v + xconnect == xafter )
        constraints.append( gammaB*v + xconnect == xbefore )

###############################################################################
### constraint: the first point after start should have same orientation as start
v = startNormalNormal
xnext = x_WS[0][1]
gammaStart = Variable(1)
constraints.append( gammaStart*v + x_WS[0][0] == xnext )

### constraint: the first point before goal should have same orientation as goal
v = goalNormalNormal
N = len(x_WS)
M = len(x_WS[N-1])
xbefore = x_WS[N-1][M-2]
xgoal = x_WS[N-1][M-1]

gammaGoal = Variable(1)
constraints.append( gammaGoal*v + xgoal == xbefore )


###############################################################################
# building objective
###############################################################################
minimaIter = 0

startMinima = 594


allValuesFirst = []
allValuesSecond = []

start = timer()

bestMinima = 0
bestMinimaValue = inf
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
        print "minima",i,"/",XspaceMinima," => ",prob.value
        allValuesFirst.append(prob.value)
        if prob.value < inf and prob.value > -inf:
                ## BICONVEX condition: check second convex problem on feasibility
                ### constraint: all points have to be inside of an environment box
                ## inbetween Rho
                ibRho = []
                mincon = []
                objective=[]
                vp = []
                for j in range(0,N_walkablesurfaces):
                        W = Wsurface_box_vstack[path[j]]
                        ibRho_tmp=[]
                        vp_tmp=[]
                        for p in range(0,len(x_WS[j])):
                                ibRho_tmp.append(Variable(XspaceDimension))
                                ## vp: normal to tangent at x_WS[j][p]
                                if p+1<len(x_WS[j]):
                                        xc = x_WS[j][p].value
                                        xn = x_WS[j][p+1].value
                                else:
                                        xc = x_WS[j][p-1].value
                                        xn = x_WS[j][p].value
                                v = xn-xc
                                v = v/np.linalg.norm(v)
                                vr = np.dot(rotFromRPY(0,0,pi/2),v)
                                vp_tmp.append(vr)

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
                print "minima",i,"/",XspaceMinima," (2nd cvx problem) => ",prob.value

                allValuesSecond.append(prob.value)
                if prob.value<inf:
                        minimaIter = i
                        print "minima",i,"admits a solution"
                        if bestMinimaValue > prob.value:
                                bestMinimaValue = prob.value
                                bestMinima = i
                        #break
        else:
                allValuesSecond.append(prob.value)
        if i%100==0:
                end = timer()
                ts= np.around(end - start,2)
                print "================================================================"
                print "Time elapsed after checking",i,"minima:"
                print ts,"s"
                print "================================================================"


end = timer()
ts= np.around(end - start,2)

print "================================================================"
print "Time elapsed for checking all",XspaceMinima,"minima"
print "================="
print ts,"s"
print "================================================================"
###############################################################################
# statistics
###############################################################################
inf = float('inf')

validMinima = np.sum(np.array(allValuesFirst) < inf)
validMinimatwo = np.sum(np.array(allValuesSecond) < inf)

pp = float(validMinima)/float(XspaceMinima)
pptwo = float(validMinimatwo)/float(XspaceMinima)

print validMinima,"of",XspaceMinima,"are valid (",pp*100,"%)"
print validMinima,"of",XspaceMinima,"second minima are valid (",pptwo*100,"%)"
print "best minima:",bestMinima,"with value",bestMinimaValue

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

        ## plot intersection environment boxes
        for i in range(0,N_walkablesurfaces-1):
                if x_connection[i].value is not None:
                        for k in range(0,maxNonzeroDim):
                                W = upperBodyConnector[i][k]
                                plot.polytopeFromVertices(\
                                                W.getVertexRepresentation(),\
                                                fcolor=colorScene)

        ### plot paths on each WS
        ###  build one path for each dimension:
        svPathPoints = 0
        for i in range(0,N_walkablesurfaces):
                svPathPoints = svPathPoints + len(x_WS[i])

        svPathsLeft = np.zeros((maxNonzeroDim, svPathPoints, 3))
        svPathsRight = np.zeros((maxNonzeroDim, svPathPoints, 3))
        svPathsMiddle = np.zeros((maxNonzeroDim, svPathPoints, 3))
        thetaV = np.zeros((svPathPoints, 5))
        for k in range(0,maxNonzeroDim):
                ctr = 0
                for i in range(0,N_walkablesurfaces):
                        for j in range(0,len(x_WS[i])):
                                svPathsMiddle[k][ctr] = x_WS[i][j].value.T
                                pt = x_WS[i][j].value+(ibRho[i][j][k].value*vp[i][j].T+heights[k]*v2).T
                                svPathsLeft[k][ctr] = np.array(pt).flatten()
                                ibRhoR = np.matrix(Ark)*ibRho[i][j].value
                                pt = x_WS[i][j].value+(np.array(ibRhoR[k]).flatten()[0]*vp[i][j].T+heights[k]*v2).T
                                svPathsRight[k][ctr] = np.array(pt).flatten()
                                ctr = ctr+1 
        ctr=0
        for i in range(0,N_walkablesurfaces):
                for j in range(0,len(x_WS[i])):
                        [k,h1,h2,h3] = Harray[minimaIter]
                        thetaV[ctr] = htoq(k,h1,h2,h3)[1]
                        ctr = ctr+1 

        svLeftFname= "data/path/xpathL.dat"
        svRightFname = "data/path/xpathR.dat"
        svMiddleFname = "data/path/xpathM.dat"
        svQValuesFname = "data/path/xpathQ.dat"


        pickle.dump( svPathsLeft, open( svLeftFname, "wb" ) )
        pickle.dump( svPathsRight, open( svRightFname, "wb" ) )
        pickle.dump( svPathsMiddle, open( svMiddleFname, "wb" ) )
        pickle.dump( thetaV, open( svQValuesFname, "wb" ) )
        plot.lines(svPathsLeft,'-or')
        plot.lines(svPathsRight,'-om')
        plot.lines(svPathsMiddle,'-ok')
        plot.set_view(90,0)
        plot.showEnvironment()

else:
        print "problem not feasible"
