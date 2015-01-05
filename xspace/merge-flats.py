import os
from scipy.spatial import ConvexHull
import random as rnd
import sys
sys.path.append("..")
from src.polytope import *
from src.linalg import *

from pylab import *
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection 
import matplotlib.pyplot as plt 
from timeit import default_timer as timer
from xspace.hspace2xspaceHRP import *
import numpy as np
import pickle

def memory_usage_psutil():
    # return the memory usage in MB
    import psutil
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0] / float(2 ** 20)
    return mem

Nr = 10000
XLname = "../data/xspacemanifold-same-axes/xsamplesL-reduced-"+str(Nr)+".dat"
XRname = "../data/xspacemanifold-same-axes/xsamplesR-reduced-"+str(Nr)+".dat"
XMname = "../data/xspacemanifold-same-axes/xsamplesM-reduced-"+str(Nr)+".dat"
Hname = "../data/xspacemanifold-same-axes/hsamples-reduced-"+str(Nr)+".dat"

XLname = "../data/xspacemanifold-same-axes/xsamplesL.dat"
XRname = "../data/xspacemanifold-same-axes/xsamplesR.dat"
XMname = "../data/xspacemanifold-same-axes/xsamplesM.dat"
Hname = "../data/xspacemanifold-same-axes/hsamples.dat"

#Xname = "../data/xspacemanifold-same-axes/xsamples.dat"
#Hname = "../data/xspacemanifold-same-axes/hsamples.dat"
HeadName = "../data/xspacemanifold-same-axes/headersamples.dat"
start = timer()
XLarray = pickle.load( open( XLname, "rb" ) )
XRarray = pickle.load( open( XRname, "rb" ) )
XMarray = pickle.load( open( XMname, "rb" ) )
Harray = pickle.load( open( Hname, "rb" ) )
[Npts, VSTACK_DELTA, heights] = pickle.load( open( HeadName, "rb" ) )
end = timer()

ts= np.around(end - start,2)
print "================================================================"
print "Time elapsed for loading Xspace configurations"
print "================="
print ts,"s"
print "================================================================"

print len(XLarray)
N_f = len(XLarray)

#print heights
Llimit = -0.5
Rlimit = 0.5

Aarray = []
ARKarray = []
barray = []
polytopes = []

omem = memory_usage_psutil()

h3cur = 0
XLarraySameH3 = []

for i in range(0,N_f):
        x = XLarray[i]
        ## build polytope from flat description
        [k,h1,h2,h3] = Harray[i]

        if h3 > h3cur:
                ## reset if h3 changes, and perform linear regression
                print "=== h3",h3,"(",len(XLarraySameH3),"samples)"
                h3cur = h3
                N = len(XLarraySameH3)

                XLarraySameH3 = sorted(XLarraySameH3,key=lambda tmp: tmp[1])
                ## check all h2 in the same group 

                ###### Split XLarraySameH3 into M arrays, whereby each subarray
                ###### has the same h2 value (could use some optimization)
                h2cur = -5
                XLarraySameH3SameH2 = []
                for j in range(0,N):
                        [xj,h2j] = XLarraySameH3[j]
                        if j==N-1:
                                XLarraySameH3SameH2.append([xj])
                        if h2j > h2cur or j==N-1:
                                if len(XLarraySameH3SameH2)>Npts:
                                        print "h2:",h2cur,"has",len(XLarraySameH3SameH2),"samples"
                                        ## plot the first two PCA's
                                        [X,Y,Z,S]=PCAprojectionOnList(XLarraySameH3SameH2)
                                        print S
                                        #projectOnSmallestSubspace(XLarraySameH3)
                                



                                h2cur = h2j
                                XLarraySameH3SameH2 = []
                                XLarraySameH3SameH2.append([xj])
                        else:
                                XLarraySameH3SameH2.append([xj])

                ## reset for next round
                #if N>0:
                        #[X,Y,Z,S]=PCAprojectionOnList(XLarraySameH3)
                        #sys.exit(0)
                XLarraySameH3 = []
                XLarraySameH3.append([x,h2])
        else:
                XLarraySameH3.append([x,h2])

        A = np.zeros((2*Npts,Npts))
        b = np.zeros((2*Npts))
        #print k,h1,h2,h3
        for j in range(0,Npts):
                e = np.zeros((Npts))
                e[j] = 1
                if heights[j] <= h3:
                        A[j,:] = e
                        b[j] = x[j]

                        A[j+Npts,:] = -e
                        b[j+Npts] = -x[j]
                else:
                        A[j,:] = e
                        b[j] = Rlimit
                        A[j+Npts,:] = -e
                        b[j+Npts] = -Llimit

        #Aarray.append(A)
        #barray.append(b)

        V = (np.dot(A,x.flatten()) <= b.flatten())

        if not V.all():
                print "[ERROR] polytope does not include its center member"
                print A,x.flatten(),b.flatten(),V
                for j in range(0,len(x)):
                        xj = x[j]
                        if not V[j]:
                                print j,V[j],xj,A[j,:],b[j]
                        if not V[j+Npts]:
                                print j+Npts,V[j+Npts],xj,A[j+Npts,:],b[j+Npts],heights[j],h3
                sys.exit(0)

        ### A defines hyperrectangle
        Aarray.append(A)
        barray.append(b.flatten())
        P = Polytope(A,b.flatten(),x.flatten())
        polytopes.append(P)

        xr = XRarray[i]
        ##find flat conversion matrix, A*x = xr

        [Ark,d] = findProjectionMatrix(x,xr)
        if d > 0.001:
                print x,xr
                print "error to big between xl and xr:",d
                sys.exit(0)

        ARKarray.append(Ark)

        ###DEBUG:
        #xspaceDisplay(x,x,x)

end = timer()
ts= np.around(end - start,2)

Aname = "../data/polytopes/A.dat"
ARKname = "../data/polytopes/Ark.dat"
bname = "../data/polytopes/b.dat"
pickle.dump( Aarray, open( Aname, "wb" ) )
pickle.dump( barray, open( bname, "wb" ) )
pickle.dump( ARKarray, open( ARKname, "wb" ) )
print "================================================================"
print "Time elapsed for computing polytopes from flats"
print "================="
print ts,"s"
print "================================================================"
print N_f,"samples have a memory footprint of",np.around(memory_usage_psutil()-omem,2),"MB"

