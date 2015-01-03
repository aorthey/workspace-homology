from timeit import default_timer as timer
import numpy as np
import pickle

Xname = "../data/xspacemanifold-same-axes/xsamples.dat"
Hname = "../data/xspacemanifold-same-axes/hsamples.dat"


Xarray = pickle.load( open( Xname, "rb" ) )
Harray = pickle.load( open( Hname, "rb" ) )

N = len(Xarray)
Nr = 100
Xnamer = "../data/xspacemanifold-same-axes/xsamples-reduced-"+str(Nr)+".dat"
Hnamer = "../data/xspacemanifold-same-axes/hsamples-reduced-"+str(Nr)+".dat"

Xarray = Xarray[0:Nr]
Harray = Harray[0:Nr]

pickle.dump( Xarray, open( Xnamer, "wb" ) )
pickle.dump( Harray, open( Hnamer, "wb" ) )
print "reduced from",N,"samples to",Nr,"samples"
