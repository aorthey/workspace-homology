import numpy as np
from numpy import dot
from cvxpy import *
from math import acos,cos,sin,atan2
from src.polytope import Polytope

def intersection(p1, p2):
        N1=p1.numberOfHalfspaces()
        N2=p2.numberOfHalfspaces()
        ##A intersection object box A_iob
        A_iob = np.zeros((N1+N2,3))
        b_iob = np.zeros((N1+N2,1))

        for j in range(0,N1):
                A_iob[j,:]=p1.A[j]
                b_iob[j]=p1.b[j]

        for j in range(0,N2):
                A_iob[j+N1,:]=p2.A[j]
                b_iob[j+N1] = p2.b[j]

        return Polytope(A_iob, b_iob)

def sortVertices2D(vertices):
        mean = np.zeros((2,1))
        V = np.zeros((len(vertices),2))
        theta = np.zeros((len(vertices),1))
        for i in range(0,len(vertices)):
                mean[0] = mean[0]+vertices[i][0]
                mean[1] = mean[1]+vertices[i][1]
        mean[0]=mean[0]/len(vertices)
        mean[1]=mean[1]/len(vertices)

        for i in range(0,len(vertices)):
                V[i,0]=vertices[i][0]
                V[i,1]=vertices[i][1]
                theta[i] = atan2(V[i,1]-mean[1],V[i,0]-mean[0])
        ## sort vertices clockwise order:
        Iv = np.argsort(theta.T)
        return V[Iv][0]

def getRotationMatrixAligningHyperplaneAndXYPlane(ap, bp):
        z=np.zeros((3,1))
        z[2]=1
        y=np.zeros((3,1))
        y[1]=1
        x=np.zeros((3,1))
        x[0]=1
        #########################
        axy = ap - (dot(ap.T,z))*z
        axynorm = np.linalg.norm(axy)
        if axynorm > 0:
                axy = axy/axynorm
                dya = dot(y.T,axy)
                if dya > 0.01:
                        txy = acos(dya)
                else:
                        txy = 0
        else:
                txy = 0
        #########################
        azy = ap - (dot(ap.T,x))*x
        azynorm = np.linalg.norm(azy)
        if azynorm > 0:
                azy = azy/azynorm
                dza = dot(z.T,azy)
                if dza > 0.01:
                        tzy = acos(dza)
                else:
                        tzy = 0
        else:
                tzy = 0
        #########################
        RX = np.zeros((3,3))
        RX[0,0]=1
        RX[1,1]=cos(txy)
        RX[1,2]=-sin(txy)
        RX[2,1]=sin(txy)
        RX[2,2]=cos(txy)

        RZ = np.zeros((3,3))
        RZ[2,2]=1
        RZ[0,0]=cos(tzy)
        RZ[0,1]=-sin(tzy)
        RZ[1,0]=sin(tzy)
        RZ[1,1]=cos(tzy)
        R = dot(RX,RZ)
        return R

def getMeanFromVerticesList(V):
        M = len(V[0])
        N = len(V)
        mean = np.zeros((M,1))
        for i in range(0,N):
                for j in range(0,M):
                        mean[j] += V[i][j]
        mean /= N
        return mean

def getMeanFromVerticesNumpy(V):
        mean = np.zeros((V.shape[1],1))
        for i in range(0,V.shape[0]):
                for j in range(0,V.shape[1]):
                        mean[j] += V[i,j]

        mean /= V.shape[0]
        return mean

def projectPointOntoHyperplane(v, a, b):
        #a=a[0]
        assert len(a)==3
        return v - (dot(v,a) - b)*a

def distancePointHyperplane(v, a, b):
        assert len(a)==3
        print v,a,b
        vprime = v - (dot(v,a) - b)*a
        return np.linalg.norm(vprime-v)

def distancePointWalkableSurface(v, W):
        xob = Variable(3)
        objective = Minimize(sum_squares(xob  - v ))

        AsurfaceX = W.ap
        bsurfaceX = W.bp
        ApolyX =    W.A
        bpolyX =    W.b

        constraints = []

        constraints.append(np.matrix(ApolyX)*xob <= bpolyX)
        constraints.append(np.matrix(AsurfaceX)*xob == bsurfaceX)

        prob = Problem(objective, constraints)
        return sqrt(abs(prob.solve())).value

def projectPointOntoPolytope(v, Ai, bi):
        xob = Variable(3)
        objective = Minimize(sum_squares(xob  - v))
        constraints = [np.matrix(Ai)*xob <= bi]
        prob = Problem(objective, constraints)
        prob.solve()
        return xob.value

def projectPointOntoWalkableSurface(v, W):
        xob = Variable(3)
        objective = Minimize(sum_squares(xob  - v ))

        AsurfaceX = W.ap
        bsurfaceX = W.bp
        ApolyX =    W.A
        bpolyX =    W.b

        constraints = []

        constraints.append(np.matrix(ApolyX)*xob <= bpolyX)
        constraints.append(np.matrix(AsurfaceX)*xob == bsurfaceX)

        prob = Problem(objective, constraints)
        prob.solve()
        return xob.value

def distancePointPolytope(v, A, b):
        xob = Variable(3)
        objective = Minimize(sum_squares(xob  - v ))
        constraints = [np.matrix(A)*xob <= b]
        prob = Problem(objective, constraints)
        return sqrt(abs(prob.solve())).value

def distancePolytopePolytope(Ai, bi, Aj, bj):
        xob = Variable(3)
        yob = Variable(3)
        objective = Minimize(sum_squares(xob  - yob ))
        constraints = [np.matrix(Ai)*xob <= bi,np.matrix(Aj)*yob <= bj]
        prob = Problem(objective, constraints)
        return sqrt(abs(prob.solve())).value

def distanceWalkableSurfacePolytope(W, A, b):
        xob = Variable(3)
        yob = Variable(3)
        objective = Minimize(sum_squares(xob  - yob ))

        AsurfaceX = W.ap
        bsurfaceX = W.bp
        ApolyX =    W.A
        bpolyX =    W.b

        constraints = []

        constraints.append(np.matrix(ApolyX)*xob <= bpolyX)
        constraints.append(np.matrix(AsurfaceX)*xob == bsurfaceX)

        constraints.append(np.matrix(A)*yob <= b)

        prob = Problem(objective, constraints)
        return sqrt(abs(prob.solve())).value

def distanceWalkableSurfaceWalkableSurface(Wi, Wj):

        xob = Variable(3)
        yob = Variable(3)
        objective = Minimize(sum_squares(xob  - yob ))

        AsurfaceX = Wi.ap
        bsurfaceX = Wi.bp
        ApolyX =    Wi.A
        bpolyX =    Wi.b

        AsurfaceY = Wj.ap
        bsurfaceY = Wj.bp
        ApolyY =    Wj.A
        bpolyY =    Wj.b

        constraints = []

        constraints.append(np.matrix(ApolyX)*xob <= bpolyX)
        constraints.append(np.matrix(AsurfaceX)*xob == bsurfaceX)

        constraints.append(np.matrix(ApolyY)*yob <= bpolyY)
        constraints.append(np.matrix(AsurfaceY)*yob == bsurfaceY)

        prob = Problem(objective, constraints)
        return sqrt(abs(prob.solve())).value

def distanceWalkableSurfaceHyperplane(W, ai, bi):
        xob = Variable(3)
        yob = Variable(3)
        objective = Minimize(sum_squares(xob  - yob ))

        AsurfaceX = W.ap
        bsurfaceX = W.bp
        ApolyX =    W.A
        bpolyX =    W.b

        constraints = []

        constraints.append(np.matrix(ApolyX)*xob <= bpolyX)
        constraints.append(np.matrix(AsurfaceX)*xob == bsurfaceX)

        #constraints.append(ai[0]*yob[0]+ai[1]*yob[1]+ai[2]*yob[2]== bi)
        constraints.append(np.matrix(ai)*yob == bi)

        prob = Problem(objective, constraints)
        d = sqrt(abs(prob.solve())).value
        xx = None
        if xob is not None:
                xx = np.zeros((3,1))
                x=xob.value
                xx[0]=x[0]
                xx[1]=x[1]
                xx[2]=x[2]
        return [d,xx]

