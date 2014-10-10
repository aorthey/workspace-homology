from BeautifulSoup import BeautifulSoup
from scipy.spatial import ConvexHull
import numpy as np
import os
import re
import math
from copy import copy 
from cvxpy import *
from pyhull.halfspace import Halfspace, HalfspaceIntersection

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
                
        def projectPointOntoHyperplane(self, v, a, b):
                a=a[0]
                return v - (v[0]*a[0]+v[1]*a[1]+v[2]*a[2] - b)*a

        def distancePolytopePolytope(self, Ai, bi, Aj, bj):
                xob = Variable(3)
                yob = Variable(3)
                objective = Minimize(sum_squares(xob  - yob ))
                constraints = [Ai*xob <= bi,Aj*yob <= bj]
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

                constraints.append(ApolyX*xob <= bpolyX)
                constraints.append(AsurfaceX*xob == bsurfaceX)

                constraints.append(Ai*yob <= bi)

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

                constraints.append(ApolyX*xob <= bpolyX)
                constraints.append(AsurfaceX*xob == bsurfaceX)

                constraints.append(ApolyY*yob <= bpolyY)
                constraints.append(AsurfaceY*yob == bsurfaceY)

                prob = Problem(objective, constraints)
                return sqrt(abs(prob.solve())).value


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


        def getWalkableSurfaces(self):
                self.W = []
                ## iterate over all objects and extract information if it is a
                ## walkable surface

                ##gravity vector
                vg = np.array((0,0,1))
                coneD = float(np.sqrt((2-2*math.cos(self.ROBOT_MAX_SLOPE*math.pi/180.0))))
                ctrW = 0
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
                                        print R.value, radius
                                        if radius >= self.ROBOT_FOOT_RADIUS:
                                                ##surface is walkable
                                                self.W.append([np.array(a),b[j],self.A[i],self.b[i],self.xyz[i]])
                                                ctrW = ctrW + 1


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
                DEBUG = 1
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
                        Ah = E[0:,0:3]
                        bh = -E[0:,3]
                        ###normalize
                        for at in range(0,len(Ah)):
                                normA = np.linalg.norm(Ah[at])
                                Ah[at] = Ah[at]/normA
                                bh[at] = bh[at]/normA
                        self.A.append(Ah)
                        self.b.append(bh)
                        self.xyz.append([x,y,z])

                np.save("xyz.simcomplex",self.xyz)

