from BeautifulSoup import BeautifulSoup
from scipy.spatial import ConvexHull
import numpy as np
import os
import re
import math
from copy import copy 
from cvxpy import *

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
                self.W=[]
                
        def computeDistanceMatrix(self):
                N = self.N
                self.M=np.zeros((N,N))
                self.D=np.zeros((N,N))
                for i in range(0,N):
                        for j in range(i+1,N):
                                xob = Variable(3)
                                yob = Variable(3)
                                objective = Minimize(sum_squares(xob  - yob ))
                                constraints = [self.A[i]*xob <= self.b[i],self.A[j]*yob <= self.b[j]]
                                prob = Problem(objective, constraints)

                                self.D[i,j] = sqrt(abs(prob.solve())).value
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
                                xob = Variable(3)
                                yob = Variable(3)
                                objective = Minimize(sum_squares(xob  - yob ))

                                AsurfaceX = self.W[i][0]
                                bsurfaceX = self.W[i][1]
                                ApolyX = self.W[i][2]
                                bpolyX = self.W[i][3]
                                AsurfaceY = self.W[j][0]
                                bsurfaceY = self.W[j][1]
                                ApolyY = self.W[j][2]
                                bpolyY = self.W[j][3]

                                constraints = []
                                constraints.append(AsurfaceX*xob <= bsurfaceX)
                                constraints.append(ApolyX*xob <= bpolyX)
                                constraints.append(AsurfaceX*yob <= bsurfaceY)
                                constraints.append(ApolyY*yob <= bpolyY)

                                prob = Problem(objective, constraints)
                                prob.solve()
                                self.WD[i,j]=self.WD[j,i]=prob.value
                                if prob.value < self.ROBOT_FOOT_RADIUS:
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
                                                self.W.append([a,b[j],self.A[i],self.b[i]])
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

                C3=[]
                for p in range(0,len(C3candidates)):
                        #[i,j,k,l]=C3candidates[p]
                        #iob = Variable(3)
                        #job = Variable(3)
                        #kob = Variable(3)
                        #lob = Variable(3)
                        #objective = Minimize(sum_squares(iob-job)+sum_squares(iob-kob)+sum_squares(iob-lob)+sum_squares(job-kob) + sum_squares(job-lob) +sum_squares(kob-lob))
                        #constraints = [A[i]*iob <= b[i],A[j]*job <= b[j],A[k]*kob <= b[k],A[l]*lob <= b[l]]
                        #prob = Problem(objective, constraints)

                        #dist = sqrt(abs(prob.solve())).value
                        #if dist < ROBOT_SPHERE_RADIUS:
                        #        C3.append([i,j,k,l])
                        C3.append([i,j,k,l])
                        print "3-cells ",p,"/",len(C3candidates)

                ###############################################################################
                ## delete subgraphs of size 4, which are fully connected, 
                ## and add two 2-cells instead
                ###############################################################################

                print C0
                print C1
                print C2
                print C3

                np.save("C0.simcomplex",C0)
                np.save("C1.simcomplex",C1)
                np.save("C2.simcomplex",C2)
                np.save("C3.simcomplex",C3)

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

