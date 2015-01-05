from pylab import *
from timeit import default_timer as timer
import sys
sys.path.append("..")
import random as rnd
from math import cos,sin,tan,pi,asin,acos,atan2,atan
import pickle
import numpy as np
from src.robotspecifications import *

from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection 
import matplotlib.pyplot as plt 

import os 

Npts = 40

heights = np.zeros((Npts,1))
heights[0]=0
for i in range(1,len(heights)):
        heights[i] = VSTACK_DELTA+heights[i-1]

def hspace2xspace(k,h1,h2,h3):
        asign = 1
        bsign = 1
        if k==0:
                asign=1
                bsign=1
        if k==1:
                asign=1
                bsign=-1
        if k==2:
                asign=-1
                bsign=-1
        if k==3:
                asign=-1
                bsign=1
        if k>3 or k<0:
                print "[WARNING] k not in {0,1,2,3}"
                return [None,None,None]

        dankle=ROBOT_DIST_FOOT_SOLE
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


        ## real limits from the hardware (but they apply relative to the
        ## previous joint)
        #qL = np.array((-1.31,-0.03,-2.18,-0.09,-0.52))
        #qU = np.array((0.73,2.62,0.73,1.05,0.79))

        ## artificial limits imposed by tangens (-pi/2,pi/2)
        tlimit = pi/3
        qaL = np.array((-tlimit,-tlimit,-tlimit,-tlimit,-tlimit))
        qaU = np.array((tlimit,tlimit,tlimit,tlimit,tlimit))

        Npts = 40
        ###############################################################################
        ## X \in \R^Npts
        ###############################################################################

        q = np.array((0.0,0.0,0.0,0.0,0.0))

        ## given h1, compute the 
        ## a : distance from knee to the main axis through foot, hip and waist.
        ## b : distance from neck to the main axis through foot, hip and waist.

        ## http://mathworld.wolfram.com/Circle-CircleIntersection.html
        l = h3 - h1 - d2
        d5 = sqrt(h2*h2+l*l)

        l0 = h1-dankle
        sqrtA = 4*l0*l0*d0*d0-pow(l0*l0-d1*d1+d0*d0,2)
        sqrtB = 4*d5*d5*d3*d3-pow(d5*d5-d4*d4+d3*d3,2)

        if sqrtA < 0:
                return [None,None,None]
        if sqrtB < 0:
                return [None,None,None]
        a=asign*0.5*(1/l0)*sqrt(sqrtA) 
        b=bsign*0.5*(1/d5)*sqrt(sqrtB) 

        if math.isnan(b):
                #print "triangle inequality d3+d4 > d5 is not fulfilled"
                #print d3,d4,d5
                return [None,None,None]
        if math.isnan(a):
                #print "triangle inequality d0+d1 > dist(foot,hip) is not fulfilled"
                #print d0,d1
                return [None,None,None]

        q = np.array((0.0,0.0,0.0,0.0,0.0))
        if abs(-a/d0) > 1 or abs(a/d1) > 1:
                print "fatal error: knee must be below foot, not allowed"
                return [None,None,None]

        q[0] = asin(-a/d0)
        q[1] = asin(a/d1)
        q[2] = 0.0

        q[3] = asin(h2/d5)+asin(b/d3)
        v = np.array((h2-d3*sin(q[3]),l-d3*cos(q[3])))
        zaxis = np.array((0,1))
        vn = v/np.linalg.norm(v)

        q[4]=acos(np.dot(zaxis,vn))*vn[0]/abs(vn[0])

        if not((q<=qaU).all() and (q>=qaL).all()):
                #print "not in range of tan configuration",q
                return [None,None,None]

        ## compute q -> x
        xL = np.zeros((Npts,1))
        xM = np.zeros((Npts,1))
        xR = np.zeros((Npts,1))

        knee_height = d0*cos(q[0])+dankle
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
        while heights[xctr] <= dankle:
                xL[xctr] = xL[0]
                xR[xctr] = xR[0]
                xM[xctr]=xM[0]
                xctr=xctr+1

        while heights[xctr] <= knee_height:
                x = (heights[xctr]-dankle)*t0
                xL[xctr] = x - dt0
                xR[xctr] = x + dt0
                xM[xctr]=x
                xctr=xctr+1

        ################################################################################
        #### knee-to-hip path
        ################################################################################
        offset = (knee_height-dankle)*t0
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

        offset = (knee_height-dankle)*t0+(hip_height-knee_height)*t1
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
        offset = (knee_height-dankle)*t0\
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
        offset = (knee_height-dankle)*t0\
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

        headpos = (knee_height-dankle)*t0\
                        +(hip_height-knee_height)*t1\
                        +(waist_height-hip_height)*t2\
                        +(neck_height-waist_height)*t3\
                        +(head_height-neck_height)*t4

        ##DEBUG:
        DEBUG=1

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
        plt.pause(0.1)

def xspaceToImage(xL,xM,xR,did):
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
        #plt.pause(0.1)
        fname = "../data/xspaceWalk/xspaceWalk"+str(did)+".png"
        savefig(fname, bbox_inches='tight')
