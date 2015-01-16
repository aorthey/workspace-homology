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
from math import fabs


#folder="xspace"
#folder="xspace/h1_0_01_h2_0_01_h3_0_01"
folder="xspace/h1_0_01_h2_0_005_h3_0_005"
XLname   = "../data/"+folder+"/xsamplesL.dat"
XRname   = "../data/"+folder+"/xsamplesR.dat"
XMname   = "../data/"+folder+"/xsamplesM.dat"
Hname    = "../data/"+folder+"/hsamples.dat"
HeadName = "../data/"+folder+"/headersamples.dat"

print XLname
XLarray = pickle.load( open( XLname, "rb" ) )
#XRarray = pickle.load( open( XRname, "rb" ) )
#XMarray = pickle.load( open( XMname, "rb" ) )
Harray = pickle.load( open( Hname, "rb" ) )
[Npts, VSTACK_DELTA, heights] = pickle.load( open( HeadName, "rb" ) )

M = len(XLarray)
N = len(XLarray[0])

Mk=0

h3value = 1.4
h1value = 0.4
h2value = 0.0

[k,h1,h2,h3] = Harray[0]
Nk = np.sum(heights<h3value)
Nk=N
def conditionOnH(k,h1,h2,h3):
        return k==0 and fabs(h3-h3value)<0.0001

for i in range(0,M):
        x = XLarray[i]
        [k,h1,h2,h3] = Harray[i]
        if conditionOnH(k,h1,h2,h3):
                Mk+=1


xN = np.zeros((Mk,Nk))
print "number of pts:",Mk,"/",M,"with dimension",Nk,"/",N

ctr=0
for i in range(0,M):
        ##choose one dimension:
        x = XLarray[i]
        [k,h1,h2,h3] = Harray[i]

        if conditionOnH(k,h1,h2,h3):
                xN[ctr,:]=x[0:Nk].T
                ctr=ctr+1

[X,Y,Z,S] = PCAprojection(xN)
print S

if len(xN) < Nk:
        print "less points then Xspace dimension"
        sys.exit(0)

fig=figure(1)
fig.clf()
ax = fig.gca(projection='3d')
#ax = fig.gca()
#f, (ax1,ax2,ax3) = plt.subplots(3, sharex=True)
ax.plot(X,Y,Z, '-or', markersize=2)
plt.show()
