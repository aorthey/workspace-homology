import sys
sys.path.append("..")
from src.robotspecifications import *
from math import cos,sin,tan,pi,asin,acos,atan2,atan,sqrt
from plotter import rotFromRPY
import math
import numpy as np
DEBUG=1
def htoq(k,h1,h2,h3):

        if k>3 or k<0:
                print "[WARNING] k not in {0,1,2,3}"
                return [None,None]

        dk=ROBOT_DIST_FOOT_SOLE
        d0=ROBOT_DIST_KNEE_FOOT
        d1=ROBOT_DIST_HIP_KNEE
        d2=ROBOT_DIST_WAIST_HIP 
        d3=ROBOT_DIST_NECK_WAIST
        d4=ROBOT_DIST_HEAD_NECK

        l1 = h1 - dk
        if l1 < 0:
                return [None,None]
        l2 = h3 - d2 - h1
        if l2 < 0:
                return [None,None]

        d5 = sqrt(h2*h2+l2*l2)
        l3 = (d5*d5 - d4*d4 + d3*d3)/(2*d5)

        sqrtA = 4*l1*l1*d0*d0-pow(l1*l1-d1*d1+d0*d0,2)
        sqrtB = 4*d5*d5*d3*d3-pow(d5*d5-d4*d4+d3*d3,2)
        if sqrtA < 0:
                return [None,None]
        if sqrtB < 0:
                return [None,None]
        a=(1/(2*l1))*sqrt(sqrtA) 
        b=(1/(2*d5))*sqrt(sqrtB) 

        l0 = sqrt(d0*d0 - a*a)
        if math.isnan(l0):
                return [None,None]

        p1 = np.array((0,0,dk))
        p3 = np.array((0,0,h1))
        p4 = np.array((0,0,(h1+d2)))
        p6 = np.array((h2,0,h3))

        p3p1 = (p3-p1)/np.linalg.norm(p3-p1)
        p6p4 = (p6-p4)/np.linalg.norm(p6-p4)

        vk = np.array((0,0,dk))
        RYpos = rotFromRPY(0,pi/2,0)
        RYneg = rotFromRPY(0,-pi/2,0)

        v0 = l0*p3p1 + a*np.dot(RYpos,p3p1)
        v0p = l0*p3p1 + a*np.dot(RYneg,p3p1)

        v1 = p3-p1-v0
        v1p = p3-p1-v0p

        v2 = np.array((0,0,d2))

        v3 = l3*p6p4 +  b*np.dot(RYpos,p6p4)
        v3p = l3*p6p4 + b*np.dot(RYneg,p6p4)

        v4 = p6-p4-v3
        v4p = p6-p4-v3p

        v0n = v0/np.linalg.norm(v0)
        v0pn = v0p/np.linalg.norm(v0p)
        v1n = v1/np.linalg.norm(v1)
        v1pn = v1p/np.linalg.norm(v1p)
        v2n = v2/np.linalg.norm(v2)

        v3n = v3/np.linalg.norm(v3)
        v3pn = v3p/np.linalg.norm(v3p)
        v4n = v4/np.linalg.norm(v4)
        v4pn = v4p/np.linalg.norm(v4p)

        xx = np.array((1,0,0))
        yy = np.array((0,1,0))
        zz = np.array((0,0,1))

        theta = np.array((0.0,0.0,0.0,0.0,0.0)) ## real c-space values
        q = np.array((0.0,0.0,0.0,0.0,0.0)) ## deviation from main axes values

        ### given all v_i vectors, we can compute the angles with the z-axis (q)
        ### and also the relative angles between two consecutive vectors (theta)
        if k==0 or k==1:
                ## foot points towards (0,1,0), meaning path p0p1p2p3 is active
                g0 = acos(np.dot(v0n.T,zz))
                g1 = acos(np.dot(v1n.T,zz))
                g1rel = acos(np.dot(v0n.T,v1n.T))
                g2rel = acos(np.dot(v1n.T,v2n.T))

                if DEBUG:
                        ### check assumption cases
                        vtest = np.dot(RYneg,v0)
                        if np.dot(v1n,vtest) < 0:
                                print "Assumption that v1 lies on the right side of v0 is violated"
                                sys.exit(0)
                        vtest = np.dot(RYpos,v1)
                        if np.dot(v2,vtest) < 0:
                                print "Assumption that v2 lies on the left side of v1 is violated"
                                sys.exit(0)

                        vtest = np.dot(RYpos,zz)
                        if np.dot(v0,vtest) < 0:
                                print "Assumption that v0 lies on the left side of zaxis is violated"
                                sys.exit(0)

                q[0] = g0
                theta[0] = -g0

                q[1] = -g1
                theta[1] = g1rel

                q[2] = 0.0
                theta[2] = -g2rel

                if k==0:
                        ## take left path, such that p0p1p2p3p4p5p6 is active
                        g3 = acos(np.dot(v3n.T,zz))
                        g4 = acos(np.dot(v4n.T,zz))
                        g4rel = acos(np.dot(v3n.T,v4n.T))

                        if v3n[0]>0:
                                q[3] = g3
                                theta[3] = g3
                        else:
                                q[3] = -g3
                                theta[3] = -g3

                        ## due to the nature of circle intersections, v4 has to
                        ## lie on the right side of v3
                        theta[4] = -g4rel
                        if v4n[0] > 0:
                                q[4] = g4
                        else:
                                q[4] = -g4
                        if DEBUG:
                                ### check assumption cases
                                vtest = np.dot(RYneg,v3n)
                                if np.dot(v4n,vtest) < 0:
                                        print "Assumption that v4 lies on the right side of v3 is violated"
                                        sys.exit(0)
                else:
                        ## take right path, such that p0p1p2p3p4p5'p6 is active
                        g3 = acos(np.dot(v3pn.T,zz))
                        g4 = acos(np.dot(v4pn.T,zz))
                        g4rel = acos(np.dot(v3pn.T,v4pn.T))

                        if v3pn[0]>0:
                                q[3] = g3
                                theta[3] = g3
                        else:
                                q[3] = -g3
                                theta[3] = -g3

                        ## v4 lies on the left side of v3
                        theta[4] = g4rel
                        if v4pn[0] > 0:
                                q[4] = g4
                        else:
                                q[4] = -g4
                        if DEBUG:
                                ### check assumption cases
                                vtest = np.dot(RYpos,v3pn)
                                if np.dot(v4pn,vtest) < 0:
                                        print "Assumption that v4p lies on the left side of v3p is violated"
                                        sys.exit(0)
        else:
                ## foot points towards (0,-1,0), meaning path p0p1p2'p3 is active
                g0 = acos(np.dot(v0pn.T,zz))
                g1 = acos(np.dot(v1pn.T,zz))
                g1rel = acos(np.dot(v0pn.T,v1pn.T))
                g2rel = acos(np.dot(v1pn.T,v2n.T))

                if DEBUG:
                        vtest = np.dot(RYpos,v0pn)
                        if np.dot(v1pn,vtest) < 0:
                                print "Assumption that v1p lies on the left side of v0p is violated"
                                sys.exit(0)
                        vtest = np.dot(RYneg,v1pn)
                        if np.dot(v2n,vtest) < 0:
                                print "Assumption that v2 lies on the right side of v1p is violated"
                                sys.exit(0)
                        vtest = np.dot(RYneg,zz)
                        if np.dot(v0pn,vtest) < 0:
                                print "Assumption that v0 lies on the right side of zaxis is violated"
                                sys.exit(0)

                q[0] = -g0
                theta[0] = -g0

                q[1] = g1
                theta[1] = g1rel

                q[2] = 0.0
                theta[2] = -g2rel

                if k==3:
                        ## take left path, such that p0p1p2'p3p4p5p6 is active
                        g3 = acos(np.dot(v3n.T,zz))
                        g4 = acos(np.dot(v4n.T,zz))
                        g4rel = acos(np.dot(v3n.T,v4n.T))

                        if v3n[0]>0:
                                q[3] = g3
                                theta[3] = -g3
                        else:
                                q[3] = -g3
                                theta[3] = g3

                        ## due to the nature of circle intersections, v4 has to
                        ## lie on the right side of v3
                        theta[4] = g4rel
                        if v4n[0] > 0:
                                q[4] = g4
                        else:
                                q[4] = -g4
                        if DEBUG:
                                ### check assumption cases
                                vtest = np.dot(RYneg,v3n)
                                if np.dot(v4n,vtest) < 0:
                                        print "Assumption that v4 lies on the right side of v3 is violated"
                                        sys.exit(0)
                else:
                        ## take right path, such that p0p1p2'p3p4p5'p6 is active
                        g3 = acos(np.dot(v3pn.T,zz))
                        g4 = acos(np.dot(v4pn.T,zz))
                        g4rel = acos(np.dot(v3pn.T,v4pn.T))

                        if v3pn[0]>0:
                                q[3] = g3
                                theta[3] = -g3
                        else:
                                q[3] = -g3
                                theta[3] = g3

                        ## v4 lies on the left side of v3
                        theta[4] = -g4rel
                        if v4pn[0] > 0:
                                q[4] = g4
                        else:
                                q[4] = -g4
                        if DEBUG:
                                ### check assumption cases
                                vtest = np.dot(RYpos,v3pn)
                                if np.dot(v4pn,vtest) < 0:
                                        print "Assumption that v4p lies on the left side of v3p is violated"
                                        sys.exit(0)

        return [q,theta]

