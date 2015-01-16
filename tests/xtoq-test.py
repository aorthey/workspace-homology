import sys
sys.path.append("..")
from xspace.xtoq import *


X = np.zeros((40,1))

[q,theta]=xtoq(X)

print q

