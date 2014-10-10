## RUN in sage
from sage.all import *
import os
import re
import sys
import numpy as np

xyz = np.load("xyz.simcomplex.npy")
C0 = np.load("C0.simcomplex.npy")
C1 = np.load("C1.simcomplex.npy")
C2 = np.load("C2.simcomplex.npy")
#C3 = np.load("C3.simcomplex.npy")

C=[]
for i in range(0,len(C0)):
        C.append(C0[i].tolist())
for i in range(0,len(C1)):
        C.append(C1[i].tolist())
for i in range(0,len(C2)):
        C.append(C2[i].tolist())
#for i in range(0,len(C3)):
        #C.append(C3[i].tolist())

S = SimplicialComplex(C)

#print S.homology()

S2 = simplicial_complexes.Sphere(2)
S3 = simplicial_complexes.Sphere(3)
print S.betti()
print S.is_isomorphic(S2,certify=true)
print S.is_isomorphic(S3,certify=true)
