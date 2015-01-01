from pylab import *
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection 
import matplotlib.pyplot as plt 
from timeit import default_timer as timer
import numpy as np
import pickle
#Xname = "../data/xspacemanifold-same-axes/xsamples.dat"
#Hname = "../data/xspacemanifold-same-axes/hsamples.dat"
Xname = "../data/xspacemanifold-same-axes/xsamples-reduced.dat"
Hname = "../data/xspacemanifold-same-axes/hsamples-reduced.dat"
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

print heights
Llimit = -0.5
Rlimit = 0.5
for i in range(0,1):
        x = Xarray[i]
        ## build polytope from flat description
        [k,h1,h2,h3] = Harray[i]
        A = []
        b = []
        print k,h1,h2,h3
        for j in range(0,Npts):
                e = np.zeros((Npts,1))
                e[j] = 1
                if heights[j] < h3:
                        A.append(e)
                        b.append(x[j])
                        A.append(-e)
                        b.append(x[j])
                else:
                        A.append(e)
                        b.append(Rlimit)
                        A.append(-e)
                        b.append(Llimit)

        fig=figure(1)
        fig.clf()
        ax = fig.gca()

        y=heights
        ax.scatter(x,y,marker='o',c='r')
        plot(x,y,'-r')

        lenlines=0.6
        plt.gca().set_aspect('equal', adjustable='box')
        for i in range(0,len(heights)):
                plot([-lenlines,lenlines],[heights[i],heights[i]],'-b')

        plt.pause(0.1)
