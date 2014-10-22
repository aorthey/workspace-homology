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

        print "------------------------------------"

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
        print ap[0],bp
        z=np.zeros((3,1))
        z[2]=1
        y=np.zeros((3,1))
        y[1]=1
        x=np.zeros((3,1))
        x[0]=1
        axy = ap - (dot(ap.T,z))*z
        axy = axy/np.linalg.norm(axy)
        azy = ap - (dot(ap.T,x))*x
        azy = azy/np.linalg.norm(azy)
        print axy,ap
        #########################
        dya = dot(y.T,axy)
        if dya > 0.01:
                txy = acos(dya)
        else:
                txy = 0
        dza = dot(z.T,azy)
        if dza > 0.01:
                tzy = acos(dza)
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

def projectPointOntoHyperplane(v, a, b):
        a=a[0]
        return v - (dot(v,a) - b)*a

def distancePointHyperplane(v, a, b):
        a=a[0]
        vprime = v - (dot(v,a) - b)*a
        return np.linalg.norm(vprime-v)

def projectPointOntoPolytope(v, Ai, bi):
        xob = Variable(3)
        objective = Minimize(sum_squares(xob  - v))
        constraints = [np.matrix(Ai)*xob <= bi]
        prob = Problem(objective, constraints)
        prob.solve()
        return xob.value

def distancePolytopePolytope(Ai, bi, Aj, bj):
        xob = Variable(3)
        yob = Variable(3)
        objective = Minimize(sum_squares(xob  - yob ))
        constraints = [np.matrix(Ai)*xob <= bi,np.matrix(Aj)*yob <= bj]
        prob = Problem(objective, constraints)
        return sqrt(abs(prob.solve())).value

def distanceWalkableSurfacePolytope(Wi, Ai, bi):
        xob = Variable(3)
        yob = Variable(3)
        objective = Minimize(sum_squares(xob  - yob ))

        AsurfaceX = Wi[0]
        bsurfaceX = Wi[1]
        ApolyX =    Wi[2]
        bpolyX =    Wi[3]

        constraints = []

        constraints.append(np.matrix(ApolyX)*xob <= bpolyX)
        constraints.append(np.matrix(AsurfaceX)*xob == bsurfaceX)

        constraints.append(np.matrix(Ai)*yob <= bi)

        prob = Problem(objective, constraints)
        return sqrt(abs(prob.solve())).value

def distanceWalkableSurfaceWalkableSurface(Wi, Wj):

        xob = Variable(3)
        yob = Variable(3)
        objective = Minimize(sum_squares(xob  - yob ))

        AsurfaceX = Wi[0]
        bsurfaceX = Wi[1]
        ApolyX =    Wi[2]
        bpolyX =    Wi[3]

        AsurfaceY = Wj[0]
        bsurfaceY = Wj[1]
        ApolyY =    Wj[2]
        bpolyY =    Wj[3]

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

