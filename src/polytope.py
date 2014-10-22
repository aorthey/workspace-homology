import numpy as np
from numpy import dot
from itertools import combinations
from sympy import *

class Polytope:

        def __init__(self, A, b, xyz=[]):
                self.A=A
                self.b=b
                self.xyz=xyz
                self.V=[]
                self.mean = np.zeros((3,1))

        def __str__(self):
                out = ""
                out += "------------------------------\n"
                out += "Polytope\n"
                out += "------------------------------\n"
                out += "A shape: "+str(self.A.shape)+"\n"
                out += "b shape: "+str(self.b.shape)+"\n"
                out += "xyz pos: "+str(self.xyz)+"\n"
                out += "------------------------------\n"
                return out

        def numberOfHalfspaces(self):
                return self.A.shape[0]

        def getMean(self):
                if not self.V:
                        self.getVertexRepresentation()
                return self.mean

        ## returns A,b
        def getHyperplaneRepresentation(self):
                return self.A,self.b

        ##returns V
        def getVertexRepresentation(self):
                M = self.A.shape[0]
                N = self.A.shape[1]

                vertices = []
                for rowlist in combinations(range(M), N):
                        Ap = self.A[np.ix_(rowlist,range(0,N))]
                        bp = self.b[np.ix_(rowlist)]
                        if np.linalg.det(Ap) != 0:
                                xp = np.linalg.solve(Ap,bp)
                                P = np.less_equal(dot(self.A,xp),self.b)
                                if P.all():
                                        vertices.append(xp)
                if len(vertices)==0:
                        #print "number of vertices for object"
                        #print self
                        #print np.around(self.A,3)
                        #print np.around(self.b,3)
                        #print "is NULL"
                        return []

                V = np.zeros((len(vertices),3))
                theta = np.zeros((len(vertices),1))
                mean = np.zeros((3,1))
                for i in range(0,len(vertices)):
                        mean[0] = mean[0]+vertices[i][0]
                        mean[1] = mean[1]+vertices[i][1]
                        mean[2] = mean[2]+vertices[i][2]

                mean[0]=mean[0]/len(vertices)
                mean[1]=mean[1]/len(vertices)
                mean[2]=mean[2]/len(vertices)

                self.mean = mean

                for i in range(0,len(vertices)):
                        V[i,0]=vertices[i][0]
                        V[i,1]=vertices[i][1]
                        V[i,2]=vertices[i][2]
                        theta[i] = atan2(V[i,1]-mean[1],V[i,0]-mean[0])

                ## sort vertices clockwise order:
                Iv = np.argsort(theta.T)
                self.V = V[Iv][0]
                return self.V
