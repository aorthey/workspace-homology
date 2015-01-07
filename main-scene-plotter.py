from src.urdfparser import URDFtoPolytopes
from plotter import Plotter
import pickle
import numpy as np
from math import sqrt
import sys
sys.path.append("..")

## COLORS
colorScene=(0.1,0.1,0.1,0.5)
colorBodyBox=(1,0,0,0.5)
colorWalkableSurface=(1,0,0,0.7)

xstartProj = pickle.load( open( "data/xstart.dat", "rb" ) )
xgoalProj = pickle.load( open( "data/xgoal.dat", "rb" ) )
paths = pickle.load( open( "data/paths.dat", "rb" ) )
Wsurfaces_decomposed = pickle.load( open( "data/wsurfaces.dat", "rb" ) )
Wsurface_box_vstack = pickle.load( open( "data/wsurfaces_vstack.dat", "rb" ) )

plot=Plotter()
env_fname = "wall.urdf"
pobjects = URDFtoPolytopes(env_fname)

for i in range(0,len(pobjects)):
        V = pobjects[i].getVertexRepresentation()

        middle = True
        for j in range(0,len(V)):
                if abs(V[j][1])>0.3:
                        middle=False

        if middle:
                plot.polytopeFromVertices(\
                                pobjects[i].getVertexRepresentation(),\
                                fcolor=colorScene)
N_w = len(Wsurfaces_decomposed)
for i in range(0,N_w):
        V = Wsurfaces_decomposed[i].getVertexRepresentation() - 0.025*np.array((0,0,1))
        #V = Wsurfaces_decomposed[i].getVertexRepresentation()
        plot.walkableSurface( \
                        V,\
                        fcolor=colorWalkableSurface, thickness=0.05)

###############################################################################
### print summary of vstacks on each surface
###############################################################################
for i in range(0,len(Wsurface_box_vstack)):
        vstack = Wsurface_box_vstack[i]
        print "WS",i,"has",len(vstack),"layers"
        for j in range(0,len(vstack)):
                hstack = vstack[j]
                #print "  layer",j,"is decomposed into",len(hstack),"boxes"
                for k in range(0,len(hstack)):
                        V = hstack[k].getVertexRepresentation()
                        middle = True
                        minDistance = 0.0
                        minY = 10000
                        maxY = -10000
                        for p in range(0,len(V)):
                                if abs(V[p][1])>0.3:
                                        middle=False
                                if V[p][0]>maxY:
                                        maxY=V[p][0]
                                if V[p][0]<minY:
                                        minY=V[p][0]
                                

                        if middle:
                                print "env[",j,"]=[",minY,",",maxY,"]"
                                plot.polytopeFromVertices(\
                                        V, fcolor=colorBodyBox)
plot.showEnvironment()
