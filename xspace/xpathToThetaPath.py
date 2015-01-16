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
from src.robotspecifications import *

dk=ROBOT_DIST_FOOT_SOLE
d0=ROBOT_DIST_KNEE_FOOT
d1=ROBOT_DIST_HIP_KNEE
d2=ROBOT_DIST_WAIST_HIP 
d3=ROBOT_DIST_NECK_WAIST
d4=ROBOT_DIST_HEAD_NECK

svLeftFname= "../data/path/xpathL.dat"
svRightFname = "../data/path/xpathR.dat"
svMiddleFname = "../data/path/xpathM.dat"
qname = "../data/path/xpathQ.dat"
svPathsLeft = pickle.load( open( svLeftFname, "rb" ) )
svPathsRight = pickle.load( open( svRightFname, "rb" ) )
svPathsMiddle = pickle.load( open( svMiddleFname, "rb" ) )
qv = pickle.load( open( qname, "rb" ) )

[XspaceDim, M, D] = svPathsLeft.shape
print "dimensions on X:",XspaceDim
print "points on workspace path:",M
pathFootR = np.zeros((M,4))
pathFootL = np.zeros((M,4))
pathTheta = np.zeros((M,5))
for i in range(0,M):
        svL = svPathsLeft[:,i,:]

        #######################################################################
        ## compute q0 from x points
        #######################################################################
        print qv[i]

        #######################################################################
        ##### compute foot position and orientation at point M
        #######################################################################
        xFm = svPathsMiddle[0,i,:]
        xFl = svPathsRight[0,i,:]
        vF = xFl - xFm

        ## clockwise 90 rot of V
        yy = np.array((0,1,0))
        xx = np.array((1,0,0))

        vFn = vF/np.linalg.norm(vF)
        drot = np.dot(vFn.T,yy)
        t = np.dot(vFn.T,xx)
        if drot > 0:
                ##counterclockwise rot
                thetaF = acos(t)
        else:
                thetaF = -acos(t)

        pathFootR[i] = [xFm[0],xFm[1],xFm[2],thetaF]

        vv = np.dot(rotFromRPY(0,0,thetaF),xx)
        vback = np.dot(rotFromRPY(0,0,pi/2),vv)
        vback = vback/np.linalg.norm(vback)
        xFleft = xFm + 0.2*vback
        pathFootR[i] = [xFm[0],xFm[1],xFm[2],thetaF]
        pathFootL[i] = [xFleft[0],xFleft[1],xFleft[2],thetaF]
        pathTheta[i] = qv[i]
        print xFm,thetaF,"(rot of (1,0,0) around z-axis, + -> counterclockwise)"
        print xFleft,thetaF,"(rot of (1,0,0) around z-axis, + -> counterclockwise)"

pathThetaName = "../data/path/pathTheta.dat"
pathFootNameR = "../data/path/pathFootR.dat"
pathFootNameL = "../data/path/pathFootL.dat"
pickle.dump( pathFootR, open( pathFootNameR, "wb" ) )
pickle.dump( pathFootL, open( pathFootNameL, "wb" ) )
pickle.dump( pathTheta, open( pathThetaName, "wb" ) )
