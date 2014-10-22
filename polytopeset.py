from BeautifulSoup import BeautifulSoup
from scipy.spatial import ConvexHull
import numpy as np
from numpy import dot
## vertex enumeration
from sympy import *
from itertools import combinations
###
import os
import re
import math
from copy import copy 
from cvxpy import *
from pyhull.halfspace import Halfspace, HalfspaceIntersection
from plotter import Plotter

DEBUG = 1

class PolytopeSet:
        """Set of polytopes, specified by A,b and optional centroid xyz"""
        ROBOT_SPHERE_RADIUS = 0.3 #the max sphere embedded in the irreducible volume
        ROBOT_MAX_SLOPE = 5 # in degrees
        ROBOT_FOOT_RADIUS = 0.15 # in m

        def __init__(self):
                self.A=[]
                self.b=[]
                self.xyz=[]
                self.D=[]
                self.M=[]
                self.WD=[]
                self.W=[] ##walkable surface

        ## Vertex enumeration problem:
        ## brute forcing algorithm: check the solutions to all 3x3 submatrices S of A: Sx=b. 
        ## If x is feasible, then it is a vertex of Ax <= b

        def getVertices(self, A, b):
                A = Matrix(A)
                b = Matrix(b+0.01) 
                M = A.rows
                N = A.cols

                vertices = []
                for rowlist in combinations(range(M), N):
                        Ap = A.extract(rowlist, range(N))
                        bp = b.extract(rowlist, [0])
                        if Ap.det() != 0:
                                xp = np.linalg.solve(Ap,bp)
                                P = np.less_equal(A*xp,b)
                                if P.all():
                                        vertices.append(xp)

                V = np.zeros((len(vertices),3))
                theta = np.zeros((len(vertices),1))

                mean = np.zeros((2,1))
                for i in range(0,len(vertices)):
                        mean[0] = mean[0]+vertices[i][0]
                        mean[1] = mean[1]+vertices[i][1]
                mean[0]=mean[0]/len(vertices)
                mean[1]=mean[1]/len(vertices)

                for i in range(0,len(vertices)):
                        V[i,0]=vertices[i][0]
                        V[i,1]=vertices[i][1]
                        V[i,2]=vertices[i][2]
                        theta[i] = atan2(V[i,1]-mean[1],V[i,0]-mean[0])

                ## sort vertices clockwise order:
                Iv = np.argsort(theta.T)
                return V[Iv][0]

        def sortVertices(self,vertices):
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

        def getRotationMatrixAligningHyperplaneAndXYPlane(self, ap, bp):
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

        def projectPointOntoHyperplane(self, v, a, b):
                a=a[0]
                return v - (dot(v,a) - b)*a

        def distancePointHyperplane(self, v, a, b):
                a=a[0]
                vprime = v - (dot(v,a) - b)*a
                return np.linalg.norm(vprime-v)

        def projectPointOntoPolytope(self, v, Ai, bi):
                xob = Variable(3)
                objective = Minimize(sum_squares(xob  - v))
                constraints = [np.matrix(Ai)*xob <= bi]
                prob = Problem(objective, constraints)
                prob.solve()
                return xob.value

        def distancePolytopePolytope(self, Ai, bi, Aj, bj):
                xob = Variable(3)
                yob = Variable(3)
                objective = Minimize(sum_squares(xob  - yob ))
                constraints = [np.matrix(Ai)*xob <= bi,np.matrix(Aj)*yob <= bj]
                prob = Problem(objective, constraints)
                return sqrt(abs(prob.solve())).value

        def distanceWalkableSurfacePolytope(self, Wi, Ai, bi):
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

        def distanceWalkableSurfaceWalkableSurface(self, Wi, Wj):

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

        def distanceWalkableSurfaceHyperplane(self, Wi, ai, bi):
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

        def computeDistanceMatrix(self):
                N = self.N
                self.M=np.zeros((N,N))
                self.D=np.zeros((N,N))
                for i in range(0,N):
                        for j in range(i+1,N):
                                self.D[i,j] = self.distancePolytopePolytope(self.A[i],self.b[i],self.A[j],self.b[j])
                                self.D[j,i] = self.D[i,j]
                                if self.D[i,j] < self.ROBOT_SPHERE_RADIUS:
                                        self.M[i,j] = 1
                                        self.M[j,i] = 1
                                else:
                                        self.M[i,j] = 0
                                        self.M[j,i] = 0
                        print i,"/",N

                self.D=np.around(self.D,3)
                print self.D

        def distanceWalkableSurfaceMatrix(self):
                if not self.W:
                        print "No walkable surfaces! Did you call getWalkableSurfaces first?"
                        exit

                N = len(self.W)
                self.WD = np.zeros((N,N))
                self.WM = np.zeros((N,N))
                for i in range(0,N):
                        self.WM[i,i]=1
                        for j in range(i+1,N):

                                self.WD[i,j]=self.WD[j,i]=self.distanceWalkableSurfaceWalkableSurface(self.W[i],self.W[j])
                                if self.WD[i,j] < self.ROBOT_FOOT_RADIUS:
                                        self.WM[i,j]=self.WM[j,i]=1

                self.WD=np.around(self.WD,3)
                print self.WD
                print self.WM

        def createWalkableSimplicialComplex(self):
                C2candidates=[]
                N = len(self.W)
                for i in range(0,N):
                        for j in range(i+1,N):
                                if self.WM[i,j]==1:
                                        for k in range(j+1,N):
                                                if self.WM[i,k]+self.WM[j,k]==2:
                                                        C2candidates.append([i,j,k])
                C0=[]
                C1=[]
                for i in range(0,N):
                        C0.append([i])
                        for j in range(i+1,N):
                                if self.WM[i,j]==1:
                                        C1.append([i,j])
                C2=[]
                for p in range(0,len(C2candidates)):
                        [i,j,k]=C2candidates[p]
                        xob = Variable(3)
                        yob = Variable(3)
                        zob = Variable(3)
                        objective = Minimize(sum_squares(xob-yob)+sum_squares(xob-yob)+sum_squares(yob-zob))
                        constraints = [A[i]*xob <= b[i],A[j]*yob <= b[j],A[k]*zob <= b[k]]
                        prob = Problem(objective, constraints)
                        dist = sqrt(abs(prob.solve())).value
                        if dist < ROBOT_SPHERE_RADIUS:
                                C2.append([i,j,k])
                        print "2-cells ",p,"/",len(C2candidates)
                print C0
                print C1
                print C2
                self.WC0 = C0
                self.WC1 = C1
                self.WC2 = C2

                np.save("WC0.simcomplex",C0)
                np.save("WC1.simcomplex",C1)
                np.save("WC2.simcomplex",C2)

        def computeProjectableObjectCandidates(self, surfaceElement):
                if surfaceElement >= len(self.W):
                        print "exceeds number of walkable surfaces!"
                        exit
                W = self.W[surfaceElement]
                Wobj = np.zeros((self.N))
                for i in range(0,self.N):
                        Wobj[i] = self.distanceWalkableSurfacePolytope(W, self.A[i], self.b[i])

                H = []
                for i in range(0,len(W[2])):
                        h = Halfspace(W[2][i], W[3][i])
                        H.append(h)
                hi = HalfspaceIntersection(H, self.projectPointOntoHyperplane(W[4], W[0],W[1]) )
                print hi.vertices
                print np.around(Wobj,3)

        def fromWalkableSurfaceComputeBoxElement(self, surfaceElement):
                RobotFootHeight = 0.1
                ##introduce some offset to remove the objects which are adjacent
                if surfaceElement >= len(self.W):
                        print "exceeds number of walkable surfaces!"
                        exit
                W = self.W[surfaceElement]
                ap = W[0]
                objectBelongingToWalkableSurface=W[5]
                apclean = np.zeros((3,1))
                apclean[0]=ap[0][0]
                apclean[1]=ap[0][1]
                apclean[2]=ap[0][2]
                bp = W[1]
                Rxy = self.getRotationMatrixAligningHyperplaneAndXYPlane(apclean,bp)
                print Rxy
                ######################################################
                ######################################################

                A_box =[]
                b_box =[]
                ##surface hyperplane, but opposite direction
                A_box.append(-ap)
                b_box.append(-bp)
                ##distance from surface hyperplane, pointing outside
                A_box.append(ap)
                b_box.append(bp+RobotFootHeight)
                print "ROBOTFOOTHEIGHT", bp+RobotFootHeight, " <<<"

                for j in range(0,len(W[2])):
                        aj = W[2][j]
                        bj = W[3][j]
                        if np.dot(ap,aj) >0.99: 
                                ##hard alignment, either
                                ##parallel or equal => discard
                                continue
                        [value, x0] = self.distanceWalkableSurfaceHyperplane(W,aj,bj)
                        if value < 0.0001:
                                #project hyperplane
                                ajp = aj - (np.dot(ap,aj) - bp)*ap

                                bjp = dot(x0.T,np.array(ajp).T)
                                A_box.append(ajp)
                                b_box.append(bjp)

                A_clean = np.zeros((len(A_box),3))
                b_clean = np.zeros((len(b_box),1))
                for j in range(0,len(A_box)):
                        A_clean[j,0] = A_box[j][0][0]
                        A_clean[j,1] = A_box[j][0][1]
                        A_clean[j,2] = A_box[j][0][2]
                        print b_box[j]
                        b_clean[j] = b_box[j]

                A_box = A_clean
                b_box = b_clean

                #############################################################
                ## compute distance between box and objects in the scene
                #############################################################
                N = len(self.A)
                WD = []
                print "-----------------------------------------------"
                print "Distance between Box over walkable surface and"
                print "object in the environment"
                print "-----------------------------------------------"

                plot = Plotter()

                proj_objects=[]
                for i in range(0,N):
                        if i==objectBelongingToWalkableSurface:
                                continue
                        A_obj = self.A[i]
                        b_obj = self.b[i]
                        ##clean b_obj
                        b_clean = np.zeros((len(b_obj),1))
                        for j in range(0,len(b_obj)):
                                b_clean[j] = b_obj[j]
                        b_obj = b_clean

                        d=self.distancePolytopePolytope(A_obj,b_obj,A_box,b_box)
                        if d < 0.0001:
                                WD.append(d)
                                N_obj=len(A_obj)
                                N_box=len(A_box)
                                ##A intersection object box A_iob
                                A_iob = np.zeros((N_obj+N_box,3))
                                b_iob = np.zeros((N_obj+N_box,1))

                                for j in range(0,N_obj):
                                        A_iob[j,:]=A_obj[j]
                                        b_iob[j]=b_obj[j]

                                for j in range(0,N_box):
                                        A_iob[j+N_obj,:]=A_box[j]
                                        b_iob[j+N_obj] = b_box[j]

                                print "------------------------------------"
                                print "Object ",i," distance ",d

                                v_obj = self.getVertices(A_obj,b_obj)
                                plot.polytopeFromVertices(v_obj)
                                v_iob = self.getVertices(A_iob,b_iob-0.001)
                                plot.polytopeFromVertices(v_iob)
                                v_iob_prime = np.zeros((len(v_iob),3))

                                for j in range(0,len(v_iob)):
                                        v_prime = self.projectPointOntoHyperplane(v_iob[j], ap, bp)
                                        v_iob_prime[j] = np.dot(Rxy,v_prime)
                                        print v_iob_prime[j][0],v_iob_prime[j][1],v_iob_prime[j][2]
                                proj_objects.append(v_iob_prime)

                print "-----------------------------------------------"
                print "Number of objects which have to be projected: ",len(WD)
                print "-----------------------------------------------"
                print np.around(WD,3)
                #######################################################
                ## write to file for convex decomposition
                #######################################################

                v_box = self.getVertices(A_box,b_box)
                v_on_surface = np.zeros((len(v_box),1))
                segmentCtr=0
                verticesCtr=0
                verticesToWrite=[]
                segmentsToWrite=[]

                ## get box vertices
                for j in range(0,len(v_box)):
                        d = self.distancePointHyperplane(v_box[j],ap,bp)
                        v_on_surface[j] = False
                        if d <= 0.02:
                                v_on_surface[j] = True

                v_box = self.getVertices(A_box,b_box)
                v_box_prime = []
                for j in range(0,len(v_box)):
                        if v_on_surface[j]:
                                v_box_prime.append(v_box[j])

                firstVertex = verticesCtr
                polygonBoxV = []
                for j in range(0,len(v_box_prime)):
                        ## use only x,y component, since we will do polygonal
                        ## decomposition
                        x = np.around(v_box_prime[j][0],2)
                        y = np.around(v_box_prime[j][1],2)
                        verticesToWrite.append([verticesCtr,x,y])

                        if j==len(v_box_prime)-1:
                                segmentsToWrite.append([segmentCtr,verticesCtr,firstVertex])
                        else:
                                segmentsToWrite.append([segmentCtr,verticesCtr,verticesCtr+1])
                        segmentCtr=segmentCtr+1
                        verticesCtr=verticesCtr+1
                        polygonBoxV.append((x,y))

                # get vertices of objects
                objectSegments=[]
                objectSegmentsNumber=0
                polygonObjArray = []
                meanProjObjects=np.zeros((len(proj_objects),2))
                for j in range(0,len(proj_objects)):
                        vp = proj_objects[j]
                        nonDoubleCtr=0
                        lastNonDouble=0
                        nonDouble = np.zeros((len(vp),1))
                        for k in range(0,len(vp)):
                                xk = np.around(vp[k][0],2)
                                yk = np.around(vp[k][1],2)
                                meanProjObjects[j][0]=meanProjObjects[j][0]+xk/len(vp)
                                meanProjObjects[j][1]=meanProjObjects[j][1]+yk/len(vp)
                                doubleV=False
                                for l in range(0,k)[::-1]:
                                        xl = np.around(vp[l][0],2)
                                        yl = np.around(vp[l][1],2)
                                        dlk = sqrt((xk-xl)**2+(yk-yl)**2).value
                                        if dlk <= 0.012:
                                                doubleV=True
                                nonDouble[k]=False
                                if not doubleV:
                                        nonDoubleCtr=nonDoubleCtr+1
                                        nonDouble[k]=True
                                        lastNonDouble=k

                        firstVertex=verticesCtr
                        polygonObjV=[]
                        for k in range(0,len(vp)):
                                if nonDouble[k]:
                                        xk = np.around(vp[k][0],2)
                                        yk = np.around(vp[k][1],2)
                                        polygonObjV.append((xk,yk))
                                        verticesToWrite.append([verticesCtr,xk,yk])
                                        if k==lastNonDouble:
                                                segmentsToWrite.append([segmentCtr,verticesCtr,firstVertex])
                                        else:
                                                segmentsToWrite.append([segmentCtr,verticesCtr,verticesCtr+1])
                                        verticesCtr=verticesCtr+1
                                        segmentCtr=segmentCtr+1

                        polygonObjArray.append(polygonObjV)

                #######################################################
                ## Create Polygons
                #######################################################
                import Polygon, Polygon.IO

                pbox = Polygon.Polygon( polygonBoxV )
                qbox = pbox
                p = []

                for j in range(0,len(polygonObjArray)):
                        pobj = Polygon.Polygon( polygonObjArray[j] )
                        qbox = qbox - pobj
                        p.append(pobj)
                print qbox
                qdecomp = qbox.triStrip()
                print qdecomp
                for j in range(0,len(qdecomp)):
                        qdecomp[j]=self.sortVertices(qdecomp[j])
                        plot.polytopeFromPolygonVertices(qdecomp[j])
                Polygon.IO.writeSVG("poly.img", qdecomp)

                #######################################################
                ## Write vertices to file
                #######################################################
                fh = open('walkable-projection.poly', 'w')

                fh.write(str(len(verticesToWrite))+" 2 0 0\n")

                for j in range(0,len(verticesToWrite)):
                        n = verticesToWrite[j][0]
                        x = verticesToWrite[j][1]
                        y = verticesToWrite[j][2]
                        fh.write(str(n)+" "+str(x)+" "+str(y)+"\n")

                fh.write(str(len(segmentsToWrite))+" 0\n")
                for j in range(0,len(segmentsToWrite)):
                        n = segmentsToWrite[j][0]
                        d = segmentsToWrite[j][1]
                        dd = segmentsToWrite[j][2]
                        fh.write(str(n)+" "+str(d)+" "+str(dd)+"\n")


                fh.write(str(len(meanProjObjects))+"\n")
                for j in range(0,len(meanProjObjects)):
                        fh.write(str(j)+" "+str(meanProjObjects[j][0])+" "+str(meanProjObjects[j][1])+"\n")

                fh.close()
                #######################################################
                #######################################################

                plot.polytopeFromVertices(v_box)
                plot.show()


        def getWalkableSurfaces(self):
                self.W = []
                ## iterate over all objects and extract information if it is a
                ## walkable surface

                ##gravity vector
                vg = np.array((0,0,1))
                coneD = float(np.sqrt((2-2*math.cos(self.ROBOT_MAX_SLOPE*math.pi/180.0))))
                ctrW = 0
                print "-----------------------------------------------"
                print "Walkable surfaces"
                print "-----------------------------------------------"
                for i in range(0,self.N):
                        ##iterate over all polytopes
                        for j in range(0,len(self.A[i])):
                                ## iterate over all surface patches and check the
                                ## two conditions on walkability
                                A = copy(self.A[i])
                                b = copy(self.b[i])
                                K = len(A)
                                a = np.matrix(A[j])
                                if np.linalg.norm(a-vg) <= coneD:
                                        ## second condition: check if we can put a foot 
                                        ## inside the surface
                                        R = Variable(1)
                                        x = Variable(3)
                                        constraints = []
                                        for k in range(0,j)+range(j+1,K):
                                                aprime = A[k] - np.dot(A[k],A[j])*A[j]
                                                anorm = np.linalg.norm(aprime)
                                                if anorm>0.001:
                                                        ## not parallel hyperplanes
                                                        aprime = aprime/anorm
                                                        v = np.dot(A[k],aprime)
                                                        constraints.append(A[k][0]*x[0]+ A[k][1]*x[1]+A[k][2]*x[2] + R*v <= b[k])

                                        constraints.append( A[j][0]*x[0]+ A[j][1]*x[1]+A[j][2]*x[2]== b[j])
                                        constraints.append(R>=0)

                                        objective = Maximize(R)
                                        self.prob = Problem(objective, constraints)
                                        solver_output = self.prob.solve(solver=ECOS)
                                        radius = self.prob.value
                                        if radius >= self.ROBOT_FOOT_RADIUS:
                                                print ctrW,": radius on surface: ",radius
                                                ##surface is walkable
                                                self.W.append([np.array(a),b[j],self.A[i],self.b[i],self.xyz[i],i])
                                                ctrW = ctrW + 1
                print "-----------------------------------------------"


        def computeSimplex(self):
                C2candidates=[]
                C3candidates=[]
                for i in range(0,N):
                        for j in range(i+1,N):
                                if M[i,j]==1:
                                        for k in range(j+1,N):
                                                if M[i,k]+M[j,k]==2:
                                                        foundFour=0
                                                        for l in range(k+1,N):
                                                                if M[i,l]+M[j,l]+M[k,l]==3:
                                                                        C3candidates.append([i,j,k,l])
                                                                        foundFour=1
                                                        if foundFour==0:
                                                                #connection between 3 cells
                                                                C2candidates.append([i,j,k])

                ###############################################################################
                ## Create simplicial complex from distance matrices
                ###############################################################################

                C0=[]
                C1=[]
                for i in range(0,N):
                        C0.append([i])
                        for j in range(i+1,N):
                                if M[i,j]==1:
                                        C1.append([i,j])
                C2=[]
                for p in range(0,len(C2candidates)):
                        [i,j,k]=C2candidates[p]
                        xob = Variable(3)
                        yob = Variable(3)
                        zob = Variable(3)
                        objective = Minimize(sum_squares(xob-yob)+sum_squares(xob-yob)+sum_squares(yob-zob))
                        constraints = [A[i]*xob <= b[i],A[j]*yob <= b[j],A[k]*zob <= b[k]]
                        prob = Problem(objective, constraints)
                        dist = sqrt(abs(prob.solve())).value
                        if dist < ROBOT_SPHERE_RADIUS:
                                C2.append([i,j,k])
                        print "2-cells ",p,"/",len(C2candidates)

                ###############################################################################
                ## delete subgraphs of size 4, which are fully connected, 
                ## and add two 2-cells instead
                ###############################################################################

                print C0
                print C1
                print C2

                np.save("C0.simcomplex",C0)
                np.save("C1.simcomplex",C1)
                np.save("C2.simcomplex",C2)

        def fromURDF(self,urdf_fname):
                soup = BeautifulSoup(open(urdf_fname))
                links = soup.robot.findAll("collision")
                K=[]
                for i in range(0,len(links)):
                        L=links[i]
                        st=L.geometry.box["size"]
                        size = re.split(' ',st)
                        sx = float(size[0])
                        sy = float(size[1])
                        sz = float(size[2])

                        pos=L.origin["xyz"]
                        pos = re.split(' ',pos)
                        ori=L.origin["rpy"]
                        ori = re.split(' ',ori)

                        x = float(pos[0])
                        y = float(pos[1])
                        z = float(pos[2])
                        ro = float(ori[0])
                        po = float(ori[1])
                        yo = float(ori[2])

                        ## prune small boxes (DEBUG MODE)
                        if DEBUG:
                                if sx+sy > 0.5:
                                        K.append([sx,sy,sz,x,y,z,ro,po,yo])
                        else:
                                K.append([sx,sy,sz,x,y,z,ro,po,yo])

                self.N = len(K)
                self.A=[]
                self.b=[]
                self.xyz=[]
                for i in range(0,self.N):
                        v=np.abs(K[i][6])+np.abs(K[i][7])+np.abs(K[i][8])
                        if v>0.001:
                                print "please do not rotate any boxes in URDF -- not handled atm"
                                exit
                        [sx,sy,sz,x,y,z] = K[i][0:6]

                        p1 = [x+sx/2, y+sy/2, z+sz/2]
                        p2 = [x+sx/2, y+sy/2, z-sz/2]
                        p3 = [x+sx/2, y-sy/2, z+sz/2]
                        p4 = [x+sx/2, y-sy/2, z-sz/2]
                        p5 = [x-sx/2, y+sy/2, z+sz/2]
                        p6 = [x-sx/2, y+sy/2, z-sz/2]
                        p7 = [x-sx/2, y-sy/2, z+sz/2]
                        p8 = [x-sx/2, y-sy/2, z-sz/2]
                        hull = ConvexHull([p1,p2,p3,p4,p5,p6,p7,p8])
                        E=hull.equations[0::2]
                        Ah = np.array(E[0:,0:3])
                        bh = np.array(-E[0:,3])
                        ###normalize
                        for at in range(0,len(Ah)):
                                normA = np.linalg.norm(Ah[at])
                                Ah[at] = Ah[at]/normA
                                bh[at] = bh[at]/normA
                        self.A.append(Ah)
                        self.b.append(bh)
                        self.xyz.append([x,y,z])

                np.save("xyz.simcomplex",self.xyz)

if __name__ == "__main__":
        p = PolytopeSet()
        p.fromURDF("wall.urdf")

        #p.computeDistanceMatrix()
        p.getWalkableSurfaces()


        p.distanceWalkableSurfaceMatrix()
        #p.createWalkableSimplicialComplex()
        #p.computeProjectableObjectCandidates(1)
        #p.computeProjectableObjectCandidates(2)
        p.fromWalkableSurfaceComputeBoxElement(1)

