# -*- mode: cython -*-

cdef extern from *:
    int exactP_list(double *r, int N, double *h, double *magvb, double *Wcrit, double *W, double *P,
                    double EPS_ABS, double EPS_REL, int ode_method)

