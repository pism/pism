#!/usr/bin/env python

# see src/verif/tests/simpleP.c
simpleP_c_result="""$ ./simpleP
Enter  r  (in km; 0 <= r <= TESTP_L = 22.5):   20.0
Results from Test P:
    h = 180.0000 (m)  Po = 16.0687800 (bar)  |vb| = 46.26644 (m a-1)
    W_c = 0.58184968 (m)  W = 0.67507258 (m)  P = 2.0086731 (bar)"""

import numpy as np
from exactP import exactP_list

secpera = 31556926.0
EPS_ABS = 1.0e-12
EPS_REL = 1.0e-15

print "setting r = 20 km ..."
r = np.array([20.0e3])

#ierr = exactP(r*1000.0,&h,&magvb,&Wcrit,&W,&P,EPS_ABS[0],EPS_REL[0],1);
h, magvb, Wcrit, W, P = exactP_list(r, EPS_ABS, EPS_REL, 1)

j = 0
print "Results from Test P:"
print "    h = %.4f (m)  Po = %.7f (bar)  |vb| = %.5f (m a-1)" \
  % (h[j], 910.0 * 9.81 * h[j] / 1.0e5, magvb[j] * secpera)
print "    W_c = %.8f (m)  W = %.8f (m)  P = %.7f (bar)" \
  % (Wcrit[j], W[j], P[j]/1.0e5)

print ""
print "compare to SAVED simpleP.c result:"
print simpleP_c_result
