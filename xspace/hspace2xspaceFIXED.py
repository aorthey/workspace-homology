from pylab import *
from timeit import default_timer as timer
import sys
sys.path.append("..")

import random as rnd
from math import cos,sin,tan,pi,asin,acos,atan2,atan
import pickle
import numpy as np
from src.robotspecifications import *
from xspace.htoq import *

from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection 
import matplotlib.pyplot as plt 
import os 

Npts = 40
heights = np.zeros((Npts,1))
heights[0]=0

##DEBUG:
DEBUG=1

for i in range(1,len(heights)):
        heights[i] = VSTACK_DELTA+heights[i-1]

def hspace2xspace(k,h1,h2,h3):

        dk=ROBOT_DIST_FOOT_SOLE
        d0=ROBOT_DIST_KNEE_FOOT
        d1=ROBOT_DIST_HIP_KNEE
        d2=ROBOT_DIST_WAIST_HIP 
        d3=ROBOT_DIST_NECK_WAIST
        d4=ROBOT_DIST_HEAD_NECK

        dtankle = ROBOT_THICKNESS_FOOT_SOLE*0.5
        dt0 = ROBOT_THICKNESS_KNEE_FOOT*0.5
        dt1 = ROBOT_THICKNESS_HIP_KNEE*0.5
        dt2 = ROBOT_THICKNESS_WAIST_HIP*0.5
        dt3 = ROBOT_THICKNESS_NECK_WAIST*0.5
        dt4 = ROBOT_THICKNESS_HEAD_NECK*0.5

        [q,theta] = htoq(k,h1,h2,h3)

        ###############################################################################
        ### check limits
        ###############################################################################

        ## real limits from the hardware (relative to the previous joint)
        thetaL = np.array((-1.31,-0.03,-2.18,-0.09,-0.52))
        thetaU = np.array((0.73,2.62,0.73,1.05,0.79))
        if not((theta<=thetaU).all() and (theta>=thetaL).all()):
                return [None,None,None]

        ## artificial limits imposed by tangens (-pi/2,pi/2)
        tlimit = pi/2
        qaL = np.array((-tlimit,-pi/2.0,-tlimit,-tlimit,-tlimit))
        qaU = np.array((tlimit,pi/2.0,tlimit,tlimit,tlimit))

        if not((q<=qaU).all() and (q>=qaL).all()):
                return [None,None,None]

        ###############################################################################
        ## X \in \R^Npts
        ###############################################################################

        ## given h1, compute the 
        ## a : distance from knee to the main axis through foot, hip and waist.
        ## b : distance from neck to the main axis through foot, hip and waist.
        ## http://mathworld.wolfram.com/Circle-CircleIntersection.html

        ## compute q -> x
        xL = np.zeros((Npts,1))
        xM = np.zeros((Npts,1))
        xR = np.zeros((Npts,1))

        knee_height = d0*cos(q[0])+dk
        hip_height = knee_height+d1*cos(q[1])
        waist_height = hip_height+d2*cos(q[2])
        neck_height = waist_height+d3*cos(q[3])
        head_height = neck_height+d4*cos(q[4])

        xctr=1

        xL[0]=-dtankle
        xR[0]=dtankle
        xM[0]=0

        t0 = tan((q[0]))
        t1 = tan((q[1]))
        t2 = tan((q[2]))
        t3 = tan((q[3]))
        t4 = tan((q[4]))

        ###############################################################################
        ### foot-to-knee path
        ###############################################################################
        while heights[xctr] <= dk:
                xL[xctr] = xL[0]
                xR[xctr] = xR[0]
                xM[xctr] = xM[0]
                xctr=xctr+1

        while heights[xctr] <= knee_height:
                x = (heights[xctr]-dk)*t0
                xL[xctr] = x - dt0
                xR[xctr] = x + dt0
                xM[xctr]=x
                xctr=xctr+1

        ################################################################################
        #### knee-to-hip path
        ################################################################################
        offset = (knee_height-dk)*t0
        kneepos = offset
        while heights[xctr] < hip_height:
                x = (heights[xctr]-knee_height)*t1+offset
                xL[xctr] = x - dt1
                xR[xctr] = x + dt1
                xM[xctr]=x
                xctr=xctr+1

        ################################################################################
        #### hip-to-waist path
        ################################################################################

        offset = (knee_height-dk)*t0+(hip_height-knee_height)*t1
        hippos = offset

        while heights[xctr] < waist_height:
                x = (heights[xctr]-hip_height)*t2+offset
                xL[xctr] = x - dt2
                xR[xctr] = x + dt2
                xM[xctr]=x
                xctr=xctr+1

        ################################################################################
        #### waist-to-neck path
        ################################################################################
        offset = (knee_height-dk)*t0\
                        +(hip_height-knee_height)*t1\
                        +(waist_height-hip_height)*t2

        waistpos = offset

        while heights[xctr] < neck_height:
                x = (heights[xctr]-waist_height)*t3+offset
                xL[xctr] = x - dt3
                xR[xctr] = x + dt3
                xM[xctr]=x
                xctr=xctr+1
        ################################################################################
        #### neck-to-head path
        ################################################################################
        offset = (knee_height-dk)*t0\
                        +(hip_height-knee_height)*t1\
                        +(waist_height-hip_height)*t2\
                        +(neck_height-waist_height)*t3
        neckpos = offset


        while xctr<len(heights) and heights[xctr] < head_height:
                x = (heights[xctr]-neck_height)*t4+offset
                xL[xctr] = x - dt4
                xR[xctr] = x + dt4
                xM[xctr]=x
                xctr=xctr+1

        headpos = (knee_height-dk)*t0\
                        +(hip_height-knee_height)*t1\
                        +(waist_height-hip_height)*t2\
                        +(neck_height-waist_height)*t3\
                        +(head_height-neck_height)*t4


        if DEBUG:
                v = hip_height-h1
                d = sqrt(v*v)
                if d > 0.001:
                        print "hip height not equal to h1, error in algorithm"
                        print hip_height,h1
                        sys.exit(0)
                v = head_height-h3
                d = sqrt(v*v)
                if d > 0.001:
                        print "head height not equal to h3, error in algorithm"
                        print head_height,h3
                        sys.exit(0)
        return [xL,xM,xR]

def xspaceDisplay(xL,xM,xR):
        xspacePlot(xL,xM,xR)
        plt.pause(0.1)

def xspacePlot(xL,xM,xR):
        fig=figure(1)
        fig.clf()
        ax = fig.gca()

        ax.scatter(xL,heights,marker='o',c='r')
        plot(xL,heights,'-r')
        ax.scatter(xR,heights,marker='o',c='r')
        plot(xR,heights,'-r')
        ax.scatter(xM,heights,marker='o',c='r')
        plot(xM,heights,'-r')
        lenlines=0.6

        plt.gca().set_aspect('equal', adjustable='box')
        for i in range(0,len(heights)):
                plot([-lenlines,lenlines],[heights[i],heights[i]],'-b')
        
        env=np.zeros((29,2))
        env[ 0 ]=[ -0.0899856966569 , 0.319985696657 ]
        env[ 1 ]=[ -0.0899856950575 , 0.319985696295 ]
        env[ 2 ]=[ -0.0299852151502 , 0.439985215036 ]
        env[ 3 ]=[ 1.40941414287e-05 , 0.449985873245 ]
        env[ 4 ]=[ 1.4103279401e-05 , 0.449985906137 ]
        env[ 5 ]=[ 1.43041187898e-05 , 0.409985696549 ]
        env[ 6 ]=[ -0.119981741564 , 0.379981741604 ]
        env[ 7 ]=[ -0.159985206983 , 0.189985207802 ]
        env[ 8 ]=[ -0.209985515515 , 0.189985515328 ]
        env[ 9 ]=[ -0.239985931476 , 0.189985931674 ]
        env[ 10 ]=[ -0.24998562944 , 0.209985629052 ]
        env[ 11 ]=[ -0.249985906096 , 0.199985903555 ]
        env[ 12 ]=[ -0.249985905854 , 0.19998590603 ]
        env[ 13 ]=[ -0.23998521525 , 0.229985214738 ]
        env[ 14 ]=[ -0.239982886739 , 0.249982886714 ]
        env[ 15 ]=[ -0.23997761043 , 0.289977610651 ]
        env[ 16 ]=[ -0.309191053821 , 0.289991116353 ]
        env[ 17 ]=[ -0.309162274999 , 0.539162274991 ]
        env[ 18 ]=[ -0.309162274994 , 0.539162274992 ]
        env[ 19 ]=[ -0.279957721233 , 0.529957721157 ]
        env[ 20 ]=[ -0.129985629093 , 0.329985629091 ]
        env[ 21 ]=[ -0.129992770764 , 0.209992770792 ]
        env[ 22 ]=[ -0.139975759854 , 0.179975760194 ]
        env[ 23 ]=[ -0.139975759858 , 0.179975760135 ]
        env[ 24 ]=[ -0.149975998286 , 0.159975997964 ]
        env[ 25 ]=[ 0.0400163308533 , 0.0599836888132 ]
        env[ 26 ]=[ 0.280016315611 , 0.299983688732 ]
        env[ 27 ]=[ 0.0200163113524 , 0.0399836886161 ]
        env[ 27 ]=[ 0.230016311961 , 0.249983688729 ]
        env[ 28 ]=[ -0.309184600543 , 0.0399845991986 ]
        env[ 28 ]=[ 0.230016311344 , 0.249983677602 ]

        for i in range(0,len(env)):
                plot([env[i][0],env[i][1]],[heights[i],heights[i]],'-r')
        
def xspaceToImage(xL,xM,xR,did):
        xspacePlot(xL,xM,xR)
        fname = "../data/xspaceWalk/xspaceWalk"+str(did)+".png"
        savefig(fname, bbox_inches='tight')
