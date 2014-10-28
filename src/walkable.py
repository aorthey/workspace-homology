import numpy as np
from numpy import dot
from scipy.spatial import ConvexHull
import math
from copy import copy 
from cvxpy import *
from src.polytope import Polytope
from src.robotspecifications import ROBOT_FOOT_RADIUS, ROBOT_MAX_SLOPE, ROBOT_SPHERE_RADIUS, ROBOT_FOOT_HEIGHT
from src.linalg import distanceWalkableSurfaceHyperplane, distancePolytopePolytope, getRotationMatrixAligningHyperplaneAndXYPlane
from src.linalg import projectPointOntoHyperplane
from src.linalg import distancePointHyperplane
from src.linalg import sortVertices2D
from src.linalg import getMeanFromVerticesNumpy
from src.linalg import distancePointWalkableSurface
from src.linalg import projectPointOntoWalkableSurface
from itertools import combinations
from math import acos,cos,sin,atan2
import Polygon, Polygon.IO

class WalkableSurface(Polytope):

        def __init__(self, ap, bp, A, b, iObject):
                self.ap = ap
                self.bp = bp
                self.A = A
                self.b = b
                self.iObject = iObject
                self.slope = np.linalg.norm(ap-np.array((0,0,1)))

        @classmethod
        def fromVertices(cls, ap, bp, V, iObject):
                ## assume vertices are in the plane
                zvalue = V[0,2]
                ap = ap
                bp = bp
                iObject = iObject
                A = np.zeros((len(V),3))
                b = np.zeros((len(V),1))

                vm = getMeanFromVerticesNumpy(V)
                vm = vm.flatten()
                for i in range(0,len(V)):
                        if i==len(V)-1:
                                vi = V[i,:]
                                vii = V[0,:]
                        else:
                                vi = V[i,:]
                                vii = V[i+1,:]

                        vvi = (vii-vi)/(np.linalg.norm(vii-vi))

                        av = (vi+dot(vm-vi,vvi)*(vvi))-vm
                        av = av/np.linalg.norm(av)
                        bv = dot(av,vi)
                        A[i,:]=av
                        b[i]=bv
                return cls(ap, bp, A, b, iObject)

        def __str__(self):
                out = ""
                out += "------------------------------\n"
                out += "Walkable Surface\n"
                out += "------------------------------\n"
                out += "A shape: "+str(self.A.shape)+"\n"
                out += "b shape: "+str(self.b.shape)+"\n"
                out += "ap, bp : "+str(self.ap)+" x == "+str(self.bp)+"\n"
                out += "object : "+str(self.iObject)+"\n"
                out += "slope  : "+str(self.slope)+"\n"
                out += "------------------------------\n"
                return out

        def getVertexRepresentation(self):
                An = np.vstack([self.A, self.ap])
                bn = np.vstack([self.b, self.bp])
                M = An.shape[0]
                N = An.shape[1]

                vertices = []
                for rowlist in combinations(range(M), N):
                        Ap = An[np.ix_(rowlist,range(0,N))]
                        bp = bn[np.ix_(rowlist)]
                        if np.linalg.det(Ap) != 0:
                                xp = np.linalg.solve(Ap,bp)
                                print xp
                                P = np.less_equal(dot(An,xp),bn)
                                d = distancePointHyperplane(xp,self.ap,self.bp)
                                if P.all() & d <= 0.001:
                                        vertices.append(xp)
                if len(vertices)==0:
                        #print "[WARNING] number of vertices for object is NULL"
                        return []

                V = np.zeros((len(vertices),3))
                theta = np.zeros((len(vertices),1))

                from src.linalg import getMeanFromVerticesList
                mean = getMeanFromVerticesList(vertices)

                for i in range(0,len(vertices)):
                        V[i,0]=vertices[i][0]
                        V[i,1]=vertices[i][1]
                        V[i,2]=vertices[i][2]
                        theta[i] = atan2(V[i,1]-mean[1],V[i,0]-mean[0])

                ## sort vertices clockwise order:
                Iv = np.argsort(theta.T)
                self.V = V[Iv][0]
                return self.V

        def createBox(self, DeltaL, DeltaU, DeltaSide=0.0):
                A_box =[]
                b_box =[]
                ##surface hyperplane, but opposite direction
                A_box.append(-self.ap)
                b_box.append(-self.bp-DeltaL)
                ##distance from surface hyperplane, pointing outside
                A_box.append(self.ap)
                b_box.append(self.bp+DeltaU)

                for j in range(0,len(self.A)):
                        aj = self.A[j]
                        bj = self.b[j]
                        if np.dot(self.ap,aj) >0.99: 
                                ##hard alignment, either
                                ##parallel or equal => discard
                                continue
                        [value, x0] = distanceWalkableSurfaceHyperplane(self,aj,bj)
                        if value < 0.0001:
                                #project hyperplane
                                ajp = aj - (dot(self.ap,aj))*self.ap
                                bjp = dot(x0.T,np.array(ajp).T)
                                A_box.append(ajp)
                                b_box.append(bjp+DeltaSide)


                A_clean = np.zeros((len(A_box),3))
                b_clean = np.zeros((len(b_box),1))
                for j in range(0,len(A_box)):
                        A_clean[j,:] = A_box[j]
                        b_clean[j] = b_box[j]

                A_box = A_clean
                b_box = b_clean

                return Polytope(A_box, b_box)


def getStartGoalWalkableSurfaces(wsurfaces, xstart, xgoal):
        N = len(wsurfaces)

        minStartI = 0
        minStartV = float("inf")
        for i in range(0,N):
                d = distancePointWalkableSurface(xstart, wsurfaces[i])
                if d < minStartV:
                        minStartV = d
                        minStartI = i
        minGoalI = 0
        minGoalV = float("inf")
        for i in range(0,N):
                d = distancePointWalkableSurface(xgoal, wsurfaces[i])
                if d < minGoalV:
                        minGoalV = d
                        minGoalI = i
        if minStartV > 0.5:
                print "[WARNING]: Distance of start contact to walkable surface seems too big"
        if minGoalV > 0.5:
                print "[WARNING]: Distance of goal contact to walkable surface seems too big"
        xstartProj = projectPointOntoWalkableSurface(xstart, wsurfaces[minStartI])
        xgoalProj = projectPointOntoWalkableSurface(xgoal, wsurfaces[minGoalI])

        return [minStartI, xstartProj, minGoalI, xgoalProj]


def WalkableSurfacesFromPolytopes(polytopes):
        W = []
        ## iterate over all objects and extract information if it is a
        ## walkable surface

        ##gravity vector
        vg = np.array((0,0,1))
        coneD = float(np.sqrt((2-2*math.cos(ROBOT_MAX_SLOPE*math.pi/180.0))))

        ctrW = 0

        print "-----------------------------------------------"
        print "Computing Walkable Surfaces from Polytopes"
        print "-----------------------------------------------"

        N = len(polytopes)
        ##iterate over all polytopes
        for i in range(0,N):
                ## iterate over all surface elements
                p = polytopes[i]
                for j in range(0,p.numberOfHalfspaces()):
                        ## iterate over all surface patches and check the
                        ## two conditions on walkability
                        A = copy(p.A)
                        b = copy(p.b)
                        K = len(A)
                        a = np.array(A[j])
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
                                prob = Problem(objective, constraints)
                                solver_output = prob.solve(solver=ECOS)
                                radius = prob.value
                                if radius >= ROBOT_FOOT_HEIGHT:
                                        #print ctrW,": radius on surface: ",radius
                                        ##surface is walkable
                                        ctrW = ctrW + 1
                                        W.append(WalkableSurface(a,b[j],p.A,p.b,i))

        print "Found",len(W),"initial walkable surfaces in environment"
        print "-----------------------------------------------"
        return W

def ProjectPolytopesDownInsideBox(polytopes, surface, box):
        ap = surface.ap
        bp = surface.bp
        A_box = box.A
        b_box = box.b

        apt = np.zeros((3,1))
        apt[0]=ap[0]
        apt[1]=ap[1]
        apt[2]=ap[2]
        Rxy = getRotationMatrixAligningHyperplaneAndXYPlane(apt,bp)

        #############################################################
        ## compute distance between box and objects in the scene
        #############################################################

        proj_objects=[]
        N = len(polytopes)
        for i in range(0,N):
                if i==surface.iObject:
                        continue
                p = polytopes[i]
                A_obj = p.A
                b_obj = p.b

                d=distancePolytopePolytope(A_obj,b_obj,A_box,b_box)
                if d < 0.001:
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

                        p_obj = Polytope(A_obj,b_obj)
                        p_iob = Polytope(A_iob,b_iob)
                        v_obj = p_obj.getVertexRepresentation()
                        v_iob = p_iob.getVertexRepresentation()
                        v_iob_prime = np.zeros((len(v_iob),3))
                        for j in range(0,len(v_iob)):
                                v_prime = projectPointOntoHyperplane(v_iob[j], ap, bp)
                                v_iob_prime[j] = dot(Rxy,v_prime)
                        proj_objects.append(v_iob_prime)

        #print "Found ",len(proj_objects)," objects to project down"

        #######################################################
        ## write to file for convex decomposition
        #######################################################

        p_box = Polytope(A_box,b_box)
        v_box = p_box.getVertexRepresentation()

        ## get the xy positions of the box vertices
        polygonBoxV = []
        for j in range(0,len(v_box)):
                x = np.around(v_box[j][0],4)
                y = np.around(v_box[j][1],4)
                z = np.around(v_box[j][2],4)
                vv = np.array((x,y,z))
                d = distancePointHyperplane(vv,ap,bp)
                if dot(ap,vv) <= bp:
                        vvp = dot(Rxy,(vv-(bp-d)*ap))
                else:
                        vvp = dot(Rxy,(vv-(bp+d)*ap))

                if np.absolute(vvp[2])>=0.001:
                        print "[WARNING1] z component has to be zero after rotation"
                        print " but is ",vvp[2]
                        print "vvp: ",vvp
                        print "vv: ",vv
                        print "ap,bp: ",ap,bp
                        print "d : ",d
                        print "Rxy: ",Rxy
                        raise "warning"

                polygonBoxV.append((vvp[0],vvp[1]))

        # get vertices of objects
        polygonObjArray = []
        for j in range(0,len(proj_objects)):
                vp = proj_objects[j]
                polygonObjV=[]
                for k in range(0,len(vp)):
                        xk = np.around(vp[k][0],2)
                        yk = np.around(vp[k][1],2)
                        zk = np.around(vp[k][2],2)

                        vv = np.array((xk,yk,zk))
                        d = distancePointHyperplane(vv,ap,bp)
                        #vvp = dot(Rxy,(vv-(bp+d)*ap))

                        if dot(ap,vv) <= bp:
                                vvp = dot(Rxy,(vv-(bp-d)*ap))
                        else:
                                vvp = dot(Rxy,(vv-(bp+d)*ap))

                        if np.absolute(vvp[2])>=0.001:
                                print "[WARNING2] z component has to be zero after rotation"
                                print " but is ",vvp[2]
                                print "vvp: ",vvp
                                print "vv: ",vv
                                print "ap,bp: ",ap,bp
                                print "d : ",d
                                print "Rxy: ",Rxy
                                raise "warning"

                        polygonObjV.append((vvp[0],vvp[1]))

                polygonObjArray.append(polygonObjV)

        #######################################################
        ## Create Polygons
        #######################################################

        pbox = Polygon.Polygon( polygonBoxV )
        qbox = pbox
        p = []

        for j in range(0,len(polygonObjArray)):
                if len(polygonObjArray[j])>0:
                        pobj = Polygon.Polygon( polygonObjArray[j] )
                        qbox = qbox - pobj

        qdecomp = qbox.triStrip()
        decompPolygons = []
        for j in range(0,len(qdecomp)):
                qdecomp[j]=sortVertices2D(qdecomp[j])
                projV = np.zeros((qdecomp[j].shape[0],qdecomp[j].shape[1]+1)); 
                projV[:,:-1] = qdecomp[j]

                for k in range(0,len(qdecomp[j])):
                        vvp = np.array((qdecomp[j][k][0],qdecomp[j][k][1],0))
                        vv = dot(np.linalg.inv(Rxy),vvp)+bp*ap
                        projV[k,:]=vv
                decompPolygons.append( projV )

        return decompPolygons
        #Polygon.IO.writeSVG("poly.img", qdecomp)
        #self.plot.polytopeFromVertices(v_box)

