import numpy as np
from pylab import *
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
def getMainPCAaxes(Xarray,h3,Npts):
        if len(Xarray)>0:
                from sklearn import preprocessing
                ## Xarray to MxN numpy array
                xx = np.zeros((Npts,len(Xarray)))
                for i in range(0,len(Xarray)):
                        for j in range(0,Npts):
                                xx[j,i] = Xarray[i][j]

                for i in range(0,Npts):
                        xd = 0
                        for j in range(0,len(Xarray)):
                                xd = xd+xx[i,j]
                        xd = xd/len(Xarray)
                        for j in range(0,len(Xarray)):
                                xx[i,j]=xx[i,j]-xd

                xx = preprocessing.scale(xx)

                [U,S,V]=np.linalg.svd(xx)
                #print np.around(S,1)
                uu = np.around(U,2)

                X1 = uu[:,0]
                X2 = uu[:,1]
                X3 = uu[:,2]

                return [X1,X2,X3]

def plotDataMainAxes(Xarray,PCaxes,h3,Npts,fname):
        if len(Xarray)>0:
                from sklearn import preprocessing
                ## Xarray to MxN numpy array
                xx = np.zeros((Npts,len(Xarray)))
                for i in range(0,len(Xarray)):
                        for j in range(0,Npts):
                                xx[j,i] = Xarray[i][j]

                for i in range(0,Npts):
                        xd = 0
                        for j in range(0,len(Xarray)):
                                xd = xd+xx[i,j]
                        xd = xd/len(Xarray)
                        for j in range(0,len(Xarray)):
                                xx[i,j]=xx[i,j]-xd

                xx = preprocessing.scale(xx)

                [X1,X2,X3] = PCaxes
                Xproj = np.zeros((3,len(Xarray)))
                for i in range(0,len(Xarray)):
                        x = np.dot(X1.T,xx[:,i])
                        y = np.dot(X2.T,xx[:,i])
                        z = np.dot(X3.T,xx[:,i])
                        Xproj[0,i] = x
                        Xproj[1,i] = y
                        Xproj[2,i] = z

                X = Xproj[0,:]
                Y = Xproj[1,:]
                Z = Xproj[2,:]

                fig=figure(1)
                fig.clf()
                ax = fig.gca(projection='3d')
                ax.scatter(X,Y,Z,marker='o',c='r',s=5)

                #ax.set_xlim3d(-6, 6)
                #ax.set_ylim3d(-6, 6)
                #ax.set_zlim3d(-1, 1)

                ax.text(-4,-3,2, "height="+str(np.around(float(h3),3)), None)
                #ax.text(-4,-5,2, "variability="+str(np.around(variability,2)), None)
                plt.show()
                savefig( fname, dpi=gcf().dpi)


