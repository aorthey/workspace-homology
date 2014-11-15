from BeautifulSoup import BeautifulSoup
from scipy.spatial import ConvexHull
import re
import numpy as np
from src.polytope import Polytope

DEBUG = 1
def URDFtoPolytopes(urdf_fname):
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
                        if sx+sy > 0.7:
                                K.append([sx,sy,sz,x,y,z,ro,po,yo])
                else:
                        K.append([sx,sy,sz,x,y,z,ro,po,yo])

        N = len(K)
        A=[]
        b=[]
        xyz=[]

        objectsInURDF = []

        for i in range(0,N):
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
                bh = np.zeros((len(Ah),1))
                #bh = np.array((-E[0:,3]))
                for k in range(0,len(Ah)):
                        bh[k] = -E[k,3]
                #bh = -E[0:,3]
                ###normalize
                for at in range(0,len(Ah)):
                        normA = np.linalg.norm(Ah[at])
                        Ah[at] = Ah[at]/normA
                        bh[at] = bh[at]/normA
                p = Polytope(Ah,bh,np.array([x,y,z]))
                objectsInURDF.append(p)

        return objectsInURDF

        #np.save("xyz.simcomplex",self.xyz)

def URDFRobotToCapsules(urdf_fname):
        soup = BeautifulSoup(open(urdf_fname))
        links = soup.robot.findAll("collision")
        for i in range(0,len(links)):
                L=links[i]
                print L

        #for i in range(0,len(links)):
        #        L=links[i]
        #        st=L.geometry.box["size"]
        #        size = re.split(' ',st)
        #        sx = float(size[0])
        #        sy = float(size[1])
        #        sz = float(size[2])

        #        pos=L.origin["xyz"]
        #        pos = re.split(' ',pos)
        #        ori=L.origin["rpy"]
        #        ori = re.split(' ',ori)

        #        x = float(pos[0])
        #        y = float(pos[1])
        #        z = float(pos[2])
        #        ro = float(ori[0])
        #        po = float(ori[1])
        #        yo = float(ori[2])

        #        ## prune small boxes (DEBUG MODE)
        #        if DEBUG:
        #                if sx+sy > 0.5:
        #                        K.append([sx,sy,sz,x,y,z,ro,po,yo])
        #        else:
        #                K.append([sx,sy,sz,x,y,z,ro,po,yo])

        #N = len(K)
        #A=[]
        #b=[]
        #xyz=[]

        #objectsInURDF = []

        #for i in range(0,N):
        #        v=np.abs(K[i][6])+np.abs(K[i][7])+np.abs(K[i][8])
        #        if v>0.001:
        #                print "please do not rotate any boxes in URDF -- not handled atm" 
        #                exit 
        #        [sx,sy,sz,x,y,z] = K[i][0:6]

        #        p1 = [x+sx/2, y+sy/2, z+sz/2]
        #        p2 = [x+sx/2, y+sy/2, z-sz/2]
        #        p3 = [x+sx/2, y-sy/2, z+sz/2]
        #        p4 = [x+sx/2, y-sy/2, z-sz/2]
        #        p5 = [x-sx/2, y+sy/2, z+sz/2]
        #        p6 = [x-sx/2, y+sy/2, z-sz/2]
        #        p7 = [x-sx/2, y-sy/2, z+sz/2]
        #        p8 = [x-sx/2, y-sy/2, z-sz/2]
        #        hull = ConvexHull([p1,p2,p3,p4,p5,p6,p7,p8])
        #        E=hull.equations[0::2]
        #        Ah = np.array(E[0:,0:3])
        #        bh = np.array(-E[0:,3])
        #        ###normalize
        #        for at in range(0,len(Ah)):
        #                normA = np.linalg.norm(Ah[at])
        #                Ah[at] = Ah[at]/normA
        #                bh[at] = bh[at]/normA
        #        p = Polytope(Ah,bh,np.array([x,y,z]))
        #        objectsInURDF.append(p)

        #return objectsInURDF

        ##np.save("xyz.simcomplex",self.xyz)
