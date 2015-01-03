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
import numpy as np
import pickle

def memory_usage_psutil():
    # return the memory usage in MB
    import psutil
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0] / float(2 ** 20)
    return mem

#Xname = "../data/xspacemanifold-same-axes/xsamples.dat"
#Hname = "../data/xspacemanifold-same-axes/hsamples.dat"
Nr = 1000
Xname = "../data/xspacemanifold-same-axes/xsamples-reduced-"+str(Nr)+".dat"
Hname = "../data/xspacemanifold-same-axes/hsamples-reduced-"+str(Nr)+".dat"
HeadName = "../data/xspacemanifold-same-axes/headersamples.dat"
start = timer()
Xarray = pickle.load( open( Xname, "rb" ) )
Harray = pickle.load( open( Hname, "rb" ) )
[Npts, VSTACK_DELTA, heights] = pickle.load( open( HeadName, "rb" ) )
end = timer()

ts= np.around(end - start,2)
print "================================================================"
print "Time elapsed for loading Xspace configurations"
print "================="
print ts,"s"
print "================================================================"

print len(Xarray)
N_f = len(Xarray)

#print heights
Llimit = -0.5
Rlimit = 0.5
#Aarray = []
#barray = []
polytopes = []

omem = memory_usage_psutil()

h3cur = 0
XarraySameH3 = []
for i in range(0,N_f):
        x = Xarray[i]
        ## build polytope from flat description
        [k,h1,h2,h3] = Harray[i]

        if h3 > h3cur:
                ## reset if h3 changes, and perform linear regression
                h3cur = h3
                ## check the size of array (if only one, then continue)
                if len(XarraySameH3)>0:
                        ctr=0
                        while heights[ctr] < h3cur:
                                ctr=ctr+1
                        ctr=ctr-1

                        ## take the points below H3
                        ## compute the dimension of the reduced space Xr
                        ## take boxes in X as points in Xr
                        ## do k-means or regression to compute the 
                        Nh = len(XarraySameH3)
                        print "h3=",h3,"|",len(XarraySameH3)," -> ",Nh,"samples"
                        print "h3=",h3,"|",len(XarraySameH3[0])," -> ",ctr,"samples"

                        d = np.zeros((Nh,Nh))
                        for j in range(0,Nh):
                                for k in range(j+1,Nh):
                                        xj = XarraySameH3[j]
                                        xk = XarraySameH3[k]
                                        v = xj-xk
                                        d[j,k] = np.dot(v.T,v)
                                        if d[j,k] <= 0.0001:
                                                del XarraySameH3[k]

                        print "reduced from ",Nh,"to",len(XarraySameH3),"samples"
                        #print np.around(d,2)


                XarraySameH3 = []
                XarraySameH3.append(x)
        else:
                XarraySameH3.append(x)

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
                print A,x.flatten(),b.flatten()

        P = Polytope(A,b.flatten(),x.flatten())
        polytopes.append(P)

        ##we start with h3 fixed, meaning we can check which points lie on a
        ## line!

        ###DEBUG:
        if 0:
                fig=figure(1)
                fig.clf()
                ax = fig.gca()

                y=heights
                #ax.scatter(x,y,marker='o',c='r')
                #plot(x,y,'-r')

                #lenlines=0.6
                #plt.gca().set_aspect('equal', adjustable='box')
                #for j in range(0,len(heights)):
                #        plot([-lenlines,lenlines],[heights[j],heights[j]],'-b')

                #plt.pause(0.1)
                ## generate random points inside the polytope
                for j in range(0,1):
                        fig.clf()
                        ax = fig.gca()
                        x_r = np.zeros((Npts))

                        for k in range(0,Npts):
                                r=b[k]
                                l=-b[k+Npts]
                                x_r[k]=rnd.uniform(l,r)
                                #print l,"<=",np.dot(A[k],x_r),",",np.dot(A[k+Npts],x_r),"<=",r," (",np.dot(A[k],x_r)<=b[k],",",np.dot(A[k+Npts],x_r)<=b[k+Npts],")"

                        ax.scatter(x_r,y,marker='o',c='r')
                        #plot(x_r,y,'-r')
                        V = (np.dot(A,x_r) <= b.flatten())
                        V = (np.dot(A,x.flatten()) <= b.flatten())
                        print V.all()

                        lenlines=0.6
                        plt.gca().set_aspect('equal', adjustable='box')
                        for i in range(0,len(heights)):
                                plot([-lenlines,lenlines],[heights[i],heights[i]],'-b')

                        plt.pause(0.001)

end = timer()
ts= np.around(end - start,2)
print "================================================================"
print "Time elapsed for computing polytopes from flats"
print "================="
print ts,"s"
print "================================================================"
print N_f,"samples have a memory footprint of",np.around(memory_usage_psutil()-omem,2),"MB"
