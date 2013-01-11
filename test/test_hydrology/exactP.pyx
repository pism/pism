# -*- mode: cython -*-
cimport numpy as np
cimport exactP_c

import numpy as np

ctypedef np.float64_t double_t
ctypedef np.int32_t int_t

def exactP_list(np.ndarray[dtype=double_t, ndim=1, mode='c'] r,
                double_t EPS_ABS, double_t EPS_REL, int_t ode_method):
    """
    h, magvb, Wcrit, W, P = exactP_list(r, EPS_ABS, EPS_REL, ode_method)
    """
    cdef np.ndarray[dtype=double_t, ndim=1, mode='c'] W, P, h, magvb, Wcrit
    cdef int N

    W = np.zeros_like(r)
    P = np.zeros_like(r)
    h = np.zeros_like(r)
    magvb = np.zeros_like(r)
    Wcrit = np.zeros_like(r)
    N = r.shape[0]

    errcode = exactP_c.exactP_list(<double*>r.data, N,
                                   <double*>h.data, <double*>magvb.data, <double*>Wcrit.data,
                                   <double*>W.data, <double*>P.data, EPS_ABS, EPS_REL, ode_method)

    if errcode != 0:
        raise Exception("exactP_list() call returned %d" % errcode)

    return h, magvb, Wcrit, W, P
