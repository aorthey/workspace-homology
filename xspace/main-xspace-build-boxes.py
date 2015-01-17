from timeit import default_timer as timer
import sys
sys.path.append("..")

import random as rnd
from math import cos,sin,tan,pi,asin,acos,atan2,atan
import pickle
import numpy as np
from src.robotspecifications import *
from xspace.hspace2xspaceFIXED import *

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

h1low = 0.2
h1high = 0.8
h2low = -0.2
h2high = 0.8
h3low = 0.7
h3high = 1.7

#h1step = 0.005
#h2step = 0.005
#h3step = 0.005

h1step = 0.05
h2step = 0.05
h3step = 0.05
VIDEO_DEBUG = 0

h3=h3low
NCtr = 0
NfeasibleCtr = 0
XLarray = []
XRarray = []
XMarray = []
Harray = []
THETAarray=[]
Qarray=[]

imgCtr=0
while h3 <= h3high:
        h3 = h3+h3step

        ## check the size of array (if only one, then continue)
        M = 0

        ##iterate over homotopy classes (HC)
        for k in range(0,4):

                h2 = h2low

                while h2<=h2high:

                        h2=h2+h2step
                        h1 = h1low

                        while h1 <= h1high:

                                NCtr = NCtr+1
                                h1 = h1+h1step

                                [xL,xM,xR,q,theta] = hspace2xspace(k,h1,h2,h3)


                                if xL is not None:
                                        if h1 >= maxH1:
                                                maxH1 = h1
                                        if h1 <= minH1:
                                                minH1 = h1
                                        if h2 >= maxH2:
                                                maxH2 = h2
                                        if h2 <= minH2:
                                                minH2 = h2
                                        if h3 >= maxH3:
                                                maxH3 = h3
                                        if h3 <= minH3:
                                                minH3 = h3

                                        NfeasibleCtr = NfeasibleCtr+1
                                        M = M+1

                                        XLarray.append(xL)
                                        XRarray.append(xR)
                                        XMarray.append(xM)
                                        THETAarray.append(theta)
                                        Qarray.append(q)
                                        Harray.append([k,h1,h2,h3])

                                        imgCtr=imgCtr+1
                                        if VIDEO_DEBUG:
                                                xspaceToImage(xL,xM,xR,imgCtr)
                                #### display x

        #print "for h1 in",hmin,hmax
        Nsamples=len(XLarray)
        mem = memory_usage_psutil()
        print "h3=",h3," and points=",M," (total:",Nsamples,", memory usage:",mem,")"
        if mem > 14000:
                print "memory limit reached:",mem
                sys.exit(1)

NfeasibleCtrReduced = len(XLarray)

XLname = "../data/xspace/xsamplesL.dat"
XMname = "../data/xspace/xsamplesM.dat"
XRname = "../data/xspace/xsamplesR.dat"
THETAname = "../data/xspace/xsamplesTheta.dat"
Qname = "../data/xspace/xsamplesQ.dat"
Hname = "../data/xspace/hsamples.dat"
HeadName = "../data/xspace/headersamples.dat"

pickle.dump( XLarray, open( XLname, "wb" ) )
pickle.dump( XMarray, open( XMname, "wb" ) )
pickle.dump( XRarray, open( XRname, "wb" ) )
pickle.dump( Harray, open( Hname, "wb" ) )
pickle.dump( Qarray, open( Qname, "wb" ) )
pickle.dump( THETAarray, open( THETAname, "wb" ) )

pickle.dump( [Npts, VSTACK_DELTA, heights], open( HeadName, "wb" ) )

print "======================================================================="
print "Samples:"
print "  N   = "+str(NCtr)
print "  N_f = "+str(NfeasibleCtr)
print "  N_f_r = "+str(NfeasibleCtrReduced)
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

if VIDEO_DEBUG and NfeasibleCtrReduced > 0:
        folder = "xspaceWalk"
        rmrfstr = "rm -rf ../data/"+folder+"/out.mp4"
        os.system(rmrfstr)
        ffmpegstr = "ffmpeg -y -framerate 8 -start_number 0 -i ../data/"+folder+"/xspaceWalk%d.png -pix_fmt yuv420p ../data/"+folder+"/out.mp4"
        vlcstr = "vlc ../data/"+folder+"/out.mp4"
        os.system(ffmpegstr)
        os.system(vlcstr)
