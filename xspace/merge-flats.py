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

folder="../data/xspace"
XLname =   folder+"/xsamplesL.dat"
XRname =   folder+"/xsamplesR.dat"
XMname =   folder+"/xsamplesM.dat"
Hname =    folder+"/hsamples.dat"
HeadName = folder+"/headersamples.dat"

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
HVarray = []

for i in range(0,N_f):
        x = XLarray[i]
        xl = XLarray[i]
        ## build polytope from flat description
        [k,h1,h2,h3] = Harray[i]

        A = np.zeros((2*Npts,Npts))
        b = np.zeros((2*Npts))
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
        [Ark,d] = findProjectionMatrix(xl,xr)

        xrProj = np.array(np.matrix(Ark)*xl)
        vtt = xrProj-xr
        d = np.dot(vtt.T,vtt)

        if d > 0.001:
                print x,xr
                print "error to big between xl and xr:",d
                sys.exit(0)

        ARKarray.append(Ark)
        HVarray.append(Harray[i])
        if i%100==0:
                print i,"/",N_f,"converted to polytopes"

end = timer()
ts= np.around(end - start,2)

Hname = "../data/polytopes/H.dat"
Aname = "../data/polytopes/A.dat"
ARKname = "../data/polytopes/Ark.dat"
bname = "../data/polytopes/b.dat"
pickle.dump( HVarray, open( Hname, "wb" ) )
pickle.dump( Aarray, open( Aname, "wb" ) )
pickle.dump( barray, open( bname, "wb" ) )
pickle.dump( ARKarray, open( ARKname, "wb" ) )
print "================================================================"
print "Time elapsed for computing polytopes from flats"
print "================="
print ts,"s"
print "================================================================"
print N_f,"samples have a memory footprint of",np.around(memory_usage_psutil()-omem,2),"MB"

