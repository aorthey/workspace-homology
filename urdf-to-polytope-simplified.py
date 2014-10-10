from polytopeset import PolytopeSet
p = PolytopeSet()
p.fromURDF("wall.urdf")
#p.computeDistanceMatrix()
p.getWalkableSurfaces()
p.distanceWalkableSurfaceMatrix()
#p.createWalkableSimplicialComplex()

