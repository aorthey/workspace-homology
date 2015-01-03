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
                self.isTopologyChanging_=False

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

        def intersectWithPolytope(self, rhs):
                Ai = self.A
                bi = self.b
                Aj = rhs.A
                bj = rhs.b
                Aij = np.vstack((Ai,Aj))
                bij = np.vstack((bi,bj))
                return Polytope(Aij,bij)

        def unionWithPolytope(self, rhs):
                Ai = self.A
                bi = self.b
                Aj = rhs.A
                bj = rhs.b
                Aij = np.vstack((Ai,Aj))
                bij = np.vstack((bi,bj))
                return Polytope(Aij,bij)


        def isTopologyChanging(self):
                return self.isTopologyChanging_

        def changesTopology(self, T):
                self.isTopologyChanging_=T

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
                                #keep care of numerical instabilities 
                                # by adding an offset to b
                                P = np.less_equal(dot(self.A,xp),self.b+0.0001)
                                if P.all():
                                        vertices.append(xp)

                if len(vertices)==0:
                        #print "[WARNING] number of vertices for object is NULL"
                        return []

                V = np.zeros((len(vertices),3))
                theta = np.zeros((len(vertices),1))

                from src.linalg import getMeanFromVerticesList
                self.mean = getMeanFromVerticesList(vertices)

                for i in range(0,len(vertices)):
                        V[i,0]=vertices[i][0]
                        V[i,1]=vertices[i][1]
                        V[i,2]=vertices[i][2]
                        theta[i] = atan2(V[i,1]-self.mean[1],V[i,0]-self.mean[0])

                ## sort vertices clockwise order:
                Iv = np.argsort(theta.T)
                self.V = V[Iv][0]
                return self.V

        def getVertexRepresentationUnsorted(self):
                M = self.A.shape[0]
                N = self.A.shape[1]

                vertices = []
                for rowlist in combinations(range(M), N):
                        Ap = self.A[np.ix_(rowlist,range(0,N))]
                        bp = self.b[np.ix_(rowlist)]
                        if np.linalg.det(Ap) != 0:
                                xp = np.linalg.solve(Ap,bp)
                                #keep care of numerical instabilities 
                                # by adding an offset to b
                                P = np.less_equal(dot(self.A,xp),self.b+0.0001)
                                if P.all():
                                        vertices.append(xp)

                if len(vertices)==0:
                        #print "[WARNING] number of vertices for object is NULL"
                        return []

                return vertices
