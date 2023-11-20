const double dx2 = dx*dx, dy2 = dy*dy, d4 = 4*dx*dy, d2 = 2*dx*dy;

/* Coefficients of the discretization of the first equation; u first, then v. */
double eq1[] = {
 0,  -c_n/dy2,  0,
 -4*c_w/dx2,  (c_s+c_n)/dy2+(4*c_w+4*c_e)/dx2,  -4*c_e/dx2,
 0,  -c_s/dy2,  0,
 c_n/d4+c_w/d2,  (c_w-c_e)/d2,  -c_n/d4-c_e/d2,
 (c_n-c_s)/d4,  0,  (c_s-c_n)/d4,
 -c_s/d4-c_w/d2,  (c_e-c_w)/d2,  c_s/d4+c_e/d2,
};

/* Coefficients of the discretization of the second equation; u first, then v. */
double eq2[] = {
 c_w/d4+c_n/d2,  (c_w-c_e)/d4,  -c_e/d4-c_n/d2,
 (c_n-c_s)/d2,  0,  (c_s-c_n)/d2,
 -c_w/d4-c_s/d2,  (c_e-c_w)/d4,  c_e/d4+c_s/d2,
 0,  -4*c_n/dy2,  0,
 -c_w/dx2,  (4*c_s+4*c_n)/dy2+(c_w+c_e)/dx2,  -c_e/dx2,
 0,  -4*c_s/dy2,  0,
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
