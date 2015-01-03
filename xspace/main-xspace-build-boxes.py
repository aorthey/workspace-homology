from timeit import default_timer as timer
import sys
sys.path.append("..")

import random as rnd
from math import cos,sin,tan,pi,asin,acos,atan2,atan
import pickle
import numpy as np
from src.robotspecifications import *
from xspace.hspace2xspaceHRP import *

from pylab import *
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection 
import matplotlib.pyplot as plt 
import os 

def memory_usage_psutil():
    # return the memory usage in MB
    import psutil
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0] / float(2 ** 20)
    return mem

dankle=ROBOT_DIST_FOOT_SOLE
d0=ROBOT_DIST_KNEE_FOOT
d1=ROBOT_DIST_HIP_KNEE
d2=ROBOT_DIST_WAIST_HIP 
d3=ROBOT_DIST_NECK_WAIST
d4=ROBOT_DIST_HEAD_NECK


start = timer()
## real limits from the hardware
#qL = np.array((-1.31,-0.03,-2.18,-0.09,-0.52))
#qU = np.array((0.73,2.62,0.73,1.05,0.79))
## artificial limits imposed by tangens (-pi/2,pi/2)
tlimit = pi/3
qaL = np.array((-tlimit,-tlimit,-tlimit,-tlimit,-tlimit))
qaU = np.array((tlimit,tlimit,tlimit,tlimit,tlimit))

###############################################################################
## X \in \R^Npts, the space of points, 
## each constraint to one environment (E) box
###############################################################################

## h3 is the height of the head, which we fix here to investigate the factors
## which define the surface which induces the same linear subspace in X
## Let h1 be the height of the hip

minH1 = 100000
maxH1 = 0
minH2 = 100000
maxH2 = 0
minH3 = 100000
maxH3 = 0

h1low = 0.305
h1high = 0.620
h2low = -0.425
h2high = 0.435
h3low = 0.950
h3high = 1.520

h1step = 0.005
h2step = 0.005
h3step = 0.005

h3=h3low
NCtr = 0
XLarray = []
XRarray = []
XMarray = []
Harray = []

while h3 < h3high:
        h3 = h3+h3step

        ##iterate over homotopy classes (HC)
        for k in range(0,4):

                h2 = h2low

                while h2<=h2high:

                        h2=h2+h2step
                        h1 = h1low

                        while h1 <= h1high:

                                NCtr = NCtr+1
                                h1 = h1+h1step

                                [xL,xM,xR] = hspace2xspace(k,h1,h2,h3)

                                if h2 >= maxH2:
                                        maxH2 = h2
                                if h2 <= minH2:
                                        minH2 = h2
                                if h1 >= maxH1:
                                        maxH1 = h1
                                if h1 <= minH1:
                                        minH1 = h1

                                if xL is not None:
                                        XLarray.append(xL)
                                        XRarray.append(xL)
                                        XMarray.append(xL)
                                        Harray.append([k,h1,h2,h3])
                                        #xspaceDisplay(xL,xM,xR)
                                #### display x

        #print "for h1 in",hmin,hmax
        Nsamples=len(XLarray)
        mem = memory_usage_psutil()
        print "h3=",h3," and points=",Nsamples," (memory usage:",mem,")"
        if mem > 14000:
                print "memory limit reached:",mem
                sys.exit(1)

NfeasibleCtr = len(XLarray)

XLname = "../data/xspacemanifold-same-axes/xsamplesL.dat"
XRname = "../data/xspacemanifold-same-axes/xsamplesR.dat"
XMname = "../data/xspacemanifold-same-axes/xsamplesM.dat"
Hname = "../data/xspacemanifold-same-axes/hsamples.dat"
HeadName = "../data/xspacemanifold-same-axes/headersamples.dat"

pickle.dump( XLarray, open( XLname, "wb" ) )
pickle.dump( XRarray, open( XRname, "wb" ) )
pickle.dump( XMarray, open( XMname, "wb" ) )
pickle.dump( Harray, open( Hname, "wb" ) )

pickle.dump( [Npts, VSTACK_DELTA, heights], open( HeadName, "wb" ) )

print "======================================================================="
print "Samples:"
print "  N   = "+str(NCtr)
print "  N_f = "+str(NfeasibleCtr)
print "======================================================================="
print "interval of H values"
print "h1=\["+str(minH1)+","+str(maxH1)+"\]"
print "h2=\["+str(minH2)+","+str(maxH2)+"\]"
print "h3=\["+str(minH3)+","+str(maxH3)+"\]"
print "======================================================================="
end = timer()

ts= np.around(end - start,2)

print "================================================================"
print "Time elapsed for building flats"
print "================="
print ts,"s"
print "================================================================"
