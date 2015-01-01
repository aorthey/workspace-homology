from timeit import default_timer as timer
import numpy as np
import pickle

Xname = "../data/xspacemanifold-same-axes/xsamples.dat"
Hname = "../data/xspacemanifold-same-axes/hsamples.dat"

Xnamer = "../data/xspacemanifold-same-axes/xsamples-reduced.dat"
Hnamer = "../data/xspacemanifold-same-axes/hsamples-reduced.dat"

Xarray = pickle.load( open( Xname, "rb" ) )
Harray = pickle.load( open( Hname, "rb" ) )

N = len(Xarray)
Nr = 100

Xarray = Xarray[0:Nr]
Harray = Harray[0:Nr]

pickle.dump( Xarray, open( Xnamer, "wb" ) )
pickle.dump( Harray, open( Hnamer, "wb" ) )
print "reduced from",N,"samples to",Nr,"samples"
