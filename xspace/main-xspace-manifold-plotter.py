import os.path
import pickle
import scipy.io
import numpy as np

from pylab import *
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

L = 900

folder = "xspacemanifold-same-axes"
while L < 2000:
        L = L+1

        Lstr = str(L).zfill(4)
        fnameX = "data/"+folder+"/cspaceX_"+Lstr+".dat"
        fnameY = "data/"+folder+"/cspaceY_"+Lstr+".dat"
        fnameZ = "data/"+folder+"/cspaceZ_"+Lstr+".dat"
        fnameS = "data/"+folder+"/cspaceEigenvalues_"+Lstr+".dat"

        if os.path.exists(fnameX):
                XX = pickle.load( open( fnameX, "rb" ) )
                YY = pickle.load( open( fnameY, "rb" ) )
                ZZ = pickle.load( open( fnameZ, "rb" ) )
                S = pickle.load( open( fnameS, "rb" ) )

                fig=figure(1)
                fig.clf()

                N=len(XX)/4
                print "showing L=",float(L)/100,"(",len(XX),"samples)"

                ax = fig.gca(projection='3d')
                ax.scatter(XX,YY,ZZ,marker='o',c='r',s=5)

                variability = sum(S[0:3])/sum(S)
                #print "first three coordinates account for",variability,"percent of variability in the data"

                ax.set_xlim3d(-3, 3)
                ax.set_ylim3d(-3, 3)
                ax.set_zlim3d(-1, 1)

                ax.text(-4,-3,2, "height="+str(np.around(float(L)/1000,3)), None)
                ax.text(-4,-5,2, "variability="+str(np.around(variability,2)), None)
                #scipy.io.savemat('scatter.mat', dict(X=XX,Y=YY,Z=ZZ))
                fnameFig = "data/"+folder+"/xspace"+Lstr+".png"
                savefig( fnameFig, dpi=gcf().dpi)

ffmpegstr = "ffmpeg -y -framerate 25 -start_number 0953 -i data/"+folder+"/xspace%04d.png -pix_fmt yuv420p data/"+folder+"/out.mp4"
vlcstr = "vlc data/"+folder+"/out.mp4"
os.system(ffmpegstr)
os.system(vlcstr)
