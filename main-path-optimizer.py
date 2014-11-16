import numpy as np
import cvxpy as cvx
import pickle
from numpy import inf,array,zeros
from cvxpy import *
from math import tan,pi

xstart= pickle.load( open( "data/xstart.dat", "rb" ) )
xgoal= pickle.load( open( "data/xgoal.dat", "rb" ) )
path_candidates = pickle.load( open( "data/paths.dat", "rb" ) )
Wsurfaces_decomposed = pickle.load( open( "data/wsurfaces.dat", "rb" ) )
Wsurface_box_vstack = pickle.load( open( "data/wsurfaces_vstack.dat", "rb" ) )

###############################################################################
# Check which path has the best chances of finding a feasible trajectory
###############################################################################
path = path_candidates[0]
N_w = len(path)
print "path over ",path,"walkable surfaces"

###############################################################################
# parameters
###############################################################################

robotMaxHeight = 1.539
robotMinKneeHeight = 0.2
robotMaxKneeHeight = 0.3
robotKneeAperture = pi/8
robotHipAperture = pi/200
robotNeckAperture = pi/16
robotHeadAperture = pi/16
robotKneeHeadAperture = pi/32
distanceWaypoints = 0.5
robotRadiusSweptVolume = 0.1

A = np.eye(3);
A(2,2)=0;
cknee = np.array((0,0,tan(robotKneeAperture)))
chip = np.array((0,0,tan(robotHipAperture)))
cneck = np.array((0,0,tan(robotNeckAperture)))
chead = np.array((0,0,tan(robotHeadAperture)))
ckneehead = np.array((0,0,tan(robotKneeHeadAperture)))

###############################################################################
# optimize over the connection between surfaces
###############################################################################

###############################################################################
# building variables
###############################################################################
## the x positions on the conenction between two walkable surfaces

x_connection = []

for i in range(0,N_w-1):
        x_connection.append(Variable(3,1))

for i in range(0,N_w):
        ## estimate how many points we need per surface

###############################################################################
# building constraints
###############################################################################

constraints = []
constraints.append( ap*x == bp )
constraints.append( ap*y == bp )

###############################################################################
# building objective
###############################################################################
objective = Minimize(sum_squares(x-y))

###############################################################################
# solve
###############################################################################
print "optimizing path from",xstart,"to",xgoal
prob = Problem(objective, constraints)
#prob.solve(solver=cvx.SCS, verbose=True)
prob.solve(solver=cvx.SCS)
#prob.solve(solver=cvx.CVXOPT)
print prob.value
