from src.polytope import Polytope
from src.walkable import WalkableSurface, WalkableSurfacesFromPolytopes, ProjectPolytopesDownInsideBox
from src.urdfparser import URDFtoPolytopes
from src.linalg import intersection
from src.robotspecifications import ROBOT_FOOT_RADIUS, ROBOT_MAX_SLOPE, ROBOT_SPHERE_RADIUS
from plotter import Plotter
plot = Plotter()

pobjects = URDFtoPolytopes("wall.urdf")
wsurfaces = WalkableSurfacesFromPolytopes(pobjects)

wboxesLow=[]
wboxesUp=[]
for i in range(0,len(wsurfaces)):
        wboxesLow.append(wsurfaces[i].createBox(0, 0.1))
        wboxesUp.append(wsurfaces[i].createBox(0.1, 1.0))


###############################################################################
# Visualize
###############################################################################
#for i in range(0,len(pobjects)):
#        plot.polytopeFromVertices(pobjects[i].getVertexRepresentation(),(1,0,0,0.05))

###############################################################################
# Project Objects in WboxesLow down
###############################################################################
plot.polytopeFromVertices(wboxesLow[1].getVertexRepresentation())

for i in range(0,1):
        plot.polytopeFromVertices(wsurfaces[i].getVertexRepresentation())
        plot.polytopeFromVertices(wboxesLow[i].getVertexRepresentation(),(0,1,0,0.2))
        #plot.polytopeFromVertices(wboxesUp[i].getVertexRepresentation(),(0,0,1,0.2))

        p = ProjectPolytopesDownInsideBox(pobjects, wsurfaces[i], wboxesLow[i])
        for i in range(0,len(p)):
                plot.polytopeFromPolygonVertices(p[i],0.1,(1,0,0,1))

plot.show()
