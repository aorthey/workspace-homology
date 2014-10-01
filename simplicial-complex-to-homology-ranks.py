from sage.all import *
import os
import re
import sys
import numpy as np

C0 = np.load("C0.simcomplex.npy")
C1 = np.load("C1.simcomplex.npy")
C2 = np.load("C2.simcomplex.npy")

C=[]
for i in range(0,len(C0)):
        C.append(C0[i].tolist())
for i in range(0,len(C1)):
        C.append(C1[i].tolist())
for i in range(0,len(C2)):
        C.append(C2[i].tolist())

S = SimplicialComplex(C)

print S.homology()


