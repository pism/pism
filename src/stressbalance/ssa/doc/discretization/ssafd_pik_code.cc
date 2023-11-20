const double dx2 = dx*dx, dy2 = dy*dy, d4 = 4*dx*dy, d2 = 2*dx*dy;

/* Coefficients of the discretization of the first equation; u first, then v. */
double eq1[] = {
 0,  -c_n*bPP/dy2,  0,
 -4*c_w*aMM/dx2,  (c_n*bPP+c_s*bMM)/dy2+(4*c_e*aPP+4*c_w*aMM)/dx2,  -4*c_e*aPP/dx2,
 0,  -c_s*bMM/dy2,  0,
 c_w*aMM*bPw/d2+c_n*aMn*bPP/d4,  (c_n*aPn*bPP-c_n*aMn*bPP)/d4+(c_w*aMM*bPP-c_e*aPP*bPP)/d2,  -c_e*aPP*bPe/d2-c_n*aPn*bPP/d4,
 (c_w*aMM*bMw-c_w*aMM*bPw)/d2+(c_n*aMM*bPP-c_s*aMM*bMM)/d4,  (c_n*aPP*bPP-c_n*aMM*bPP-c_s*aPP*bMM+c_s*aMM*bMM)/d4+(c_e*aPP*bPP-c_w*aMM*bPP-c_e*aPP*bMM+c_w*aMM*bMM)/d2,  (c_e*aPP*bPe-c_e*aPP*bMe)/d2+(c_s*aPP*bMM-c_n*aPP*bPP)/d4,
 -c_w*aMM*bMw/d2-c_s*aMs*bMM/d4,  (c_s*aMs*bMM-c_s*aPs*bMM)/d4+(c_e*aPP*bMM-c_w*aMM*bMM)/d2,  c_e*aPP*bMe/d2+c_s*aPs*bMM/d4,
};

/* Coefficients of the discretization of the second equation; u first, then v. */
double eq2[] = {
 c_w*aMM*bPw/d4+c_n*aMn*bPP/d2,  (c_n*aPn*bPP-c_n*aMn*bPP)/d2+(c_w*aMM*bPP-c_e*aPP*bPP)/d4,  -c_e*aPP*bPe/d4-c_n*aPn*bPP/d2,
 (c_w*aMM*bMw-c_w*aMM*bPw)/d4+(c_n*aMM*bPP-c_s*aMM*bMM)/d2,  (c_n*aPP*bPP-c_n*aMM*bPP-c_s*aPP*bMM+c_s*aMM*bMM)/d2+(c_e*aPP*bPP-c_w*aMM*bPP-c_e*aPP*bMM+c_w*aMM*bMM)/d4,  (c_e*aPP*bPe-c_e*aPP*bMe)/d4+(c_s*aPP*bMM-c_n*aPP*bPP)/d2,
 -c_w*aMM*bMw/d4-c_s*aMs*bMM/d2,  (c_s*aMs*bMM-c_s*aPs*bMM)/d2+(c_e*aPP*bMM-c_w*aMM*bMM)/d4,  c_e*aPP*bMe/d4+c_s*aPs*bMM/d2,
 0,  -4*c_n*bPP/dy2,  0,
 -c_w*aMM/dx2,  (4*c_n*bPP+4*c_s*bMM)/dy2+(c_e*aPP+c_w*aMM)/dx2,  -c_e*aPP/dx2,
 0,  -4*c_s*bMM/dy2,  0,
};

/* i indices */
const int I[] = {
 i-1,  i,  i+1,
 i-1,  i,  i+1,
 i-1,  i,  i+1,
 i-1,  i,  i+1,
 i-1,  i,  i+1,
 i-1,  i,  i+1,
};

/* j indices */
const int J[] = {
 j+1,  j+1,  j+1,
 j,  j,  j,
 j-1,  j-1,  j-1,
 j+1,  j+1,  j+1,
 j,  j,  j,
 j-1,  j-1,  j-1,
};

/* component indices */
const int C[] = {
 0,  0,  0,
 0,  0,  0,
 0,  0,  0,
 1,  1,  1,
 1,  1,  1,
 1,  1,  1,
};
