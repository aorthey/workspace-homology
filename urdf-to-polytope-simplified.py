from polytopeset import PolytopeSet
p = PolytopeSet()
p.fromURDF("wall.urdf")
p.computeDistanceMatrix()
p.getWalkableSurfaces()
p.distanceWalkableSurfaceMatrix()
#p.createWalkableSimplicialComplex()
#p.computeProjectableObjectCandidates(1)
#p.computeProjectableObjectCandidates(2)
p.fromWalkableSurfaceComputeBoxElement(1)

