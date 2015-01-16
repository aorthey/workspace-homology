import os
from scipy.spatial import ConvexHull
import random as rnd
import sys
sys.path.append("..")
from src.polytope import *
from src.linalg import *
import networkx as nx   
import random

from pylab import *
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection 
import matplotlib.pyplot as plt 
from timeit import default_timer as timer
from xspace.hspace2xspaceHRP import *
import numpy as np
import pickle
from math import fabs

HeadName = "../data/xspace/headersamples.dat"
Aname    = "../data/polytopes/A.dat"
ARKname  = "../data/polytopes/Ark.dat"
bname    = "../data/polytopes/b.dat"
HName    = "../data/polytopes/H.dat"
[N, VSTACK_DELTA, heights] = pickle.load( open( HeadName, "rb" ) )

Harray = pickle.load( open( HName, "rb" ) )
Aflat = pickle.load( open( Aname, "rb" ) )
Aleftrightconv = pickle.load( open( ARKname, "rb" ) )
bflat = pickle.load( open( bname, "rb" ) )

M = len(Aflat)

def connected(P1,P2):
        d = distancePolytopePolytope(P1,P2)
        if d < 0.5:
                return True
        else:
                return False

def buildDistanceMatrixFromPolytopes(Aflat, bflat, Dfname):
        print "infering graph from",M,"samples"
        D = np.ones((M,M))
        for i in range(0,M):
                D[i,i]=0
                for j in range(i+1,M):
                        P1 = Polytope(Aflat[i],bflat[i])
                        P2 = Polytope(Aflat[j],bflat[j])
                        d = distancePolytopePolytope(P1,P2)
                        D[i,j]=d
                        D[j,i]=d

        pickle.dump( D, open( Dfname, "wb" ) )

def getDistanceMatrixFromFile(Dfname):
        D = pickle.load( open( Dfname, "rb" ) )
        return D

Dfname = "../data/graph/distMatrix.dat"
#buildDistanceMatrixFromPolytopes(Aflat,bflat,Dfname)
D = getDistanceMatrixFromFile(Dfname)
M = D.shape[0]

G = nx.Graph()
for i in range(0,M):
        for j in range(i+1,M):
                if D[i,j]<0.2:
                        G.add_edge(i,j)

#pos=nx.random_layout(G)
#pos=nx.spectral_layout(G)
pos=nx.spring_layout(G,dim=2)
nx.draw(G,pos)
nx.draw_networkx_edges(G,pos,font_size=30)
import pylab as plt
plt.show()


