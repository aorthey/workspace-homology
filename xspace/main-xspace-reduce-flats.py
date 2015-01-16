from timeit import default_timer as timer
import numpy as np
import pickle

XLname = "../data/xspacemanifold-same-axes/xsamplesL.dat"
XRname = "../data/xspacemanifold-same-axes/xsamplesR.dat"
XMname = "../data/xspacemanifold-same-axes/xsamplesM.dat"
Hname = "../data/xspacemanifold-same-axes/hsamples.dat"
HeadName = "../data/xspacemanifold-same-axes/headersamples.dat"

print "loading samples..."
XLarray = pickle.load( open( XLname, "rb" ) )
print "loaded",XLname
XRarray = pickle.load( open( XRname, "rb" ) )
print "loaded",XRname
XMarray = pickle.load( open( XMname, "rb" ) )
print "loaded",XMname
Harray = pickle.load( open( Hname, "rb" ) )
print "loaded",Hname

N = len(XLarray)
Nr = 100
XLnamer = "../data/xspacemanifold-same-axes/xsamplesL-reduced-"+str(Nr)+".dat"
XRnamer = "../data/xspacemanifold-same-axes/xsamplesR-reduced-"+str(Nr)+".dat"
XMnamer = "../data/xspacemanifold-same-axes/xsamplesM-reduced-"+str(Nr)+".dat"
Hnamer = "../data/xspacemanifold-same-axes/hsamples-reduced-"+str(Nr)+".dat"

XLarray = XLarray[0:Nr]
XRarray = XRarray[0:Nr]
XMarray = XMarray[0:Nr]
Harray = Harray[0:Nr]

print "writing reduced files..."
pickle.dump( XLarray, open( XLnamer, "wb" ) )
pickle.dump( XRarray, open( XRnamer, "wb" ) )
pickle.dump( XMarray, open( XMnamer, "wb" ) )
pickle.dump( Harray, open( Hnamer, "wb" ) )
print "reduced from",N,"samples to",Nr,"samples"
