from BeautifulSoup import BeautifulSoup
import os
import re
import sys

DEBUG = 1

###############################################################################
## Parse URDF, extract box informations for each link
###############################################################################
urdf_fname = "wall.urdf"

ROBOT_SPHERE_RADIUS = 0.3
## header for URDF
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

N = len(K)
print "complexity ",N," over 2 for distance computation"

###############################################################################
## Take box values, convert to polytope representation
###############################################################################

from scipy.spatial import ConvexHull
A=[]
b=[]
for i in range(0,N):
        if abs(K[i][6])+abs(K[i][7])+abs(K[i][8]) > 0.001:
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
        A.append(Ah)
        b.append(bh)

###############################################################################
## Use convex optimization to compute distances between objects
###############################################################################
from cvxopt.base import *
from cvxpy import *
import numpy as np

M=np.zeros((N,N))
D=np.zeros((N,N))
for i in range(0,N):
        for j in range(i+1,N):
                xob = Variable(3)
                yob = Variable(3)
                objective = Minimize(sum_squares(xob  - yob ))
                constraints = [A[i]*xob <= b[i],A[j]*yob <= b[j]]
                prob = Problem(objective, constraints)

                D[i,j] = sqrt(abs(prob.solve())).value
                D[j,i] = D[i,j]
                if D[i,j] < ROBOT_SPHERE_RADIUS:
                        M[i,j] = 1
                        M[j,i] = 1
                else:
                        M[i,j] = 0
                        M[j,i] = 0

D=np.around(D,3)
print D
###############################################################################
## Use convex optimization to compute possible 2-cell candidates
###############################################################################
C2candidates=[]
for i in range(0,N):
        for j in range(i+1,N):
                if M[i,j]==1:
                        for k in range(j+1,N):
                                if M[i,k]+M[j,k]==2:
                                        #connection between 3 cells
                                        C2candidates.append([i,j,k])



###############################################################################
## Create simplicial complex from distance matrices
###############################################################################

C0=[]
C1=[]
for i in range(0,N):
        C0.append(i)
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

print C0
print C1
print C2
