from src.polytope import Polytope
from src.walkable import WalkableSurface, WalkableSurfacesFromPolytopes, ProjectPolytopesDownInsideBox
from src.urdfparser import URDFtoPolytopes
from src.linalg import intersection
from src.robotspecifications import ROBOT_FOOT_RADIUS, ROBOT_MAX_SLOPE, ROBOT_SPHERE_RADIUS
from src.robotspecifications import ROBOT_FOOT_HEIGHT
from src.robotspecifications import ROBOT_VOLUME_MIN_HEIGHT, ROBOT_VOLUME_MAX_HEIGHT
from plotter import Plotter
plot = Plotter()

pobjects = URDFtoPolytopes("wall.urdf")
wsurfaces = WalkableSurfacesFromPolytopes(pobjects)

wboxesLow=[]
for i in range(0,len(wsurfaces)):
        wboxesLow.append(wsurfaces[i].createBox(0, ROBOT_FOOT_HEIGHT))

###############################################################################
# Visualize Complete Scene
###############################################################################
for i in range(0,len(pobjects)):
        plot.polytopeFromVertices(pobjects[i].getVertexRepresentation(),(0,0,0,0.05))

###############################################################################
# Project Objects in WboxesLow down
###############################################################################

Wsurfaces_decomposed = []
for i in range(0,len(wsurfaces)):

        ap = wsurfaces[i].ap
        bp = wsurfaces[i].bp
        iObject = wsurfaces[i].iObject
        p = ProjectPolytopesDownInsideBox(pobjects, wsurfaces[i], wboxesLow[i])

        ###############################################################################
        # Project Objects in Upper Box above splitted surfaces Down
        ###############################################################################

        for j in range(0,len(p)):
                Wsplit = WalkableSurface.fromVertices(ap,bp,p[j],iObject)
                plot.polytopeFromPolygonVertices3D(Wsplit.getVertexRepresentation(),0.1,(1,0,0,0.6))

                WsplitBox = Wsplit.createBox(ROBOT_FOOT_HEIGHT,ROBOT_VOLUME_MIN_HEIGHT)

                p2 = ProjectPolytopesDownInsideBox(pobjects, wsurfaces[i], WsplitBox)
                for k in range(0,len(p2)):
                        W2split = WalkableSurface.fromVertices(ap,bp,p2[k],iObject)
                        Wsurfaces_decomposed.append(W2split)

###############################################################################
# Visualize traversable surfaces
###############################################################################
for i in range(0,len(Wsurfaces_decomposed)):
        plot.polytopeFromPolygonVertices3D(Wsurfaces_decomposed[i].getVertexRepresentation(),0.2,(0,1,0,0.6))

###############################################################################
# Find topology of traversable surfaces
###############################################################################



plot.show()
