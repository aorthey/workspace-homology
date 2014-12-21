import random as rnd
from math import cos,sin,tan,pi
import numpy as np
from src.robotspecifications import *

from pylab import *
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

d0=ROBOT_DIST_KNEE_FOOT
d1=ROBOT_DIST_HIP_KNEE
d2=ROBOT_DIST_WAIST_HIP 
d3=ROBOT_DIST_NECK_WAIST
d4=ROBOT_DIST_HEAD_NECK

q = np.array((0.0,0.0,0.0,0.0,0.0))

#q = np.array((-pi/16,+pi/16,-pi/4,pi/16,0.8))


## real limits from the hardware
qL = np.array((-1.31,-0.03,-2.18,-0.09,-0.52))
qU = np.array((0.73,2.62,0.73,1.05,0.79))

## artificial limits imposed by tangens
tlimit = pi/3
qaL = np.array((-tlimit,-tlimit,-tlimit,-tlimit,-tlimit))
qaU = np.array((tlimit,tlimit,tlimit,tlimit,tlimit))

##generate uniform random config qL <= q <= qU
Nsamples=100000
qarray=[]
for ctr in range(0,Nsamples):
        for i in range(0,len(qL)):
                q[i]=rnd.uniform(qL[i], qU[i])

        if not((q<=qU).all() and (q>=qL).all()):
                #print "not feasible configuration"
                continue
        if not((q<=qaU).all() and (q>=qaL).all()):
                #print "not in range of tan configuration"
                continue
        ## compute q -> x

        Npts = 40
        x = np.zeros((Npts,1))
        heights = np.zeros((Npts,1))
        heights[0]=0
        heights[1]=ROBOT_FOOT_HEIGHT

        for i in range(1,len(heights)):
                heights[i] = VSTACK_DELTA+heights[i-1]

        knee_height = d0*cos(q[0])
        hip_height = knee_height+d1*cos(q[1])
        waist_height = hip_height+d2*cos(q[2])
        neck_height = waist_height+d3*cos(q[3])
        head_height = neck_height+d4*cos(q[4])

        xctr=1
        xd=0

        x[0]=0
        t0 = tan((q[0]))
        t1 = tan((q[1]))
        t2 = tan((q[2]))
        t3 = tan((q[3]))
        t4 = tan((q[4]))
        ###############################################################################
        ### foot-to-knee path
        ###############################################################################
        while heights[xctr] <= knee_height:
                x[xctr] = heights[xctr]*t0
                xctr=xctr+1

        ################################################################################
        #### knee-to-hip path
        ################################################################################
        offset = knee_height*t0
        kneepos = offset
        while heights[xctr] < hip_height:
                x[xctr] = (heights[xctr]-knee_height)*t1+offset
                xctr=xctr+1

        ################################################################################
        #### hip-to-waist path
        ################################################################################

        offset = knee_height*t0+(hip_height-knee_height)*t1
        hippos = offset

        while heights[xctr] < waist_height:
                x[xctr] = (heights[xctr]-hip_height)*t2+offset
                xctr=xctr+1

        ################################################################################
        #### waist-to-neck path
        ################################################################################
        offset = knee_height*t0\
                        +(hip_height-knee_height)*t1\
                        +(waist_height-hip_height)*t2

        waistpos = offset

        while heights[xctr] < neck_height:
                x[xctr] = (heights[xctr]-waist_height)*t3+offset
                xctr=xctr+1
        ################################################################################
        #### neck-to-head path
        ################################################################################
        offset = knee_height*t0\
                        +(hip_height-knee_height)*t1\
                        +(waist_height-hip_height)*t2\
                        +(neck_height-waist_height)*t3
        neckpos = offset


        while xctr<len(heights) and heights[xctr] < head_height:
                x[xctr] = (heights[xctr]-neck_height)*t4+offset
                xctr=xctr+1

        headpos = knee_height*t0\
                        +(hip_height-knee_height)*t1\
                        +(waist_height-hip_height)*t2\
                        +(neck_height-waist_height)*t3\
                        +(head_height-neck_height)*t4

        if abs(headpos) >= 0.05:
                continue
        if abs(hippos) >= 0.005:
                continue
        if abs(waistpos) >= 0.005:
                continue

        qarray.append(q)
        ### display x
        fig=figure(1)
        fig.clf()
        ax = fig.gca()

        y=heights
        ax.scatter(x,y,marker='o',c='r')
        plot(x,y,'-r')
        lenlines=0.6
        plot([-lenlines,lenlines],[head_height,head_height],'-k',linewidth=3.0)
        plot([-lenlines,lenlines],[neck_height,neck_height],'-k',linewidth=2.0)
        plot([-lenlines,lenlines],[waist_height,waist_height],'-k',linewidth=2.0)
        plot([-lenlines,lenlines],[hip_height,hip_height],'-k',linewidth=2.0)
        plot([-lenlines,lenlines],[knee_height,knee_height],'-k',linewidth=2.0)

        for i in range(0,len(heights)):
                plot([-lenlines,lenlines],[heights[i],heights[i]],'-b')
        #print q
        plt.pause(0.0001)

print len(qarray),"samples from Q (from",Nsamples,"random shots)"
