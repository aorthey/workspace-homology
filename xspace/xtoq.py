import sys
sys.path.append("..")
from src.robotspecifications import *
from math import cos,sin,tan,pi,asin,acos,atan2,atan,sqrt
from plotter import rotFromRPY
import math
import numpy as np

def xtoq(X):

        dk=ROBOT_DIST_FOOT_SOLE
        d0=ROBOT_DIST_KNEE_FOOT
        d1=ROBOT_DIST_HIP_KNEE
        d2=ROBOT_DIST_WAIST_HIP 
        d3=ROBOT_DIST_NECK_WAIST
        d4=ROBOT_DIST_HEAD_NECK

        theta = np.array((0.0,0.0,0.0,0.0,0.0))
        q = np.array((0.0,0.0,0.0,0.0,0.0))

        N = len(X)
        print "dimension X:",N

        x0 = X[0]
        if len(x0)>1:
                print "wrong input",x0
                sys.exit(0)

        delta = VSTACK_DELTA

        zz = np.array((0,0,1))
        zzf = zz+x0
        vk = dk*zz
        while i*delta < dk:
                if d(x[i],zzf) > 0.001:
                        ##change to ankle knee segment
                        xp = x[i-1]


        return [q,theta]
