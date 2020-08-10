#include <cmath>

#include "pism/util/Vector2.hh"

namespace pism {

Vector2 exact_xy(double x, double y) {
  return {
    exp(x)*sin(2*M_PI*y),
    exp(x)*cos(2*M_PI*y)
  };
}

Vector2 source_xy(double x, double y, double B) {
  double x0 = 2*M_PI*y;
  double x1 = cos(x0);
  double x2 = pow(x1, 2);
  double x3 = pow(M_PI, 2);
  double x4 = pow(M_PI, 3);
  double x5 = pow(M_PI, 4);
  double x6 = sin(x0);
  double x7 = pow(x6, 2);
  double x8 = pow(2, 2.0/3.0)*B*exp((1.0/3.0)*x)/pow(x2*(4*x3 + 1 + 4*M_PI) + x7*(16*x3 - 8*M_PI + 4), 4.0/3.0);
  return {
    -1.0/3.0*x6*x8*(x2*(6*x3 + 20*x4 + 72*x5 - 2 + 3*M_PI) + x7*(-48*x3 + 32*x4 + 96*x5 - 8 + 36*M_PI)),
    -1.0/6.0*x1*x8*(x2*(-12*x3 + 136*x4 + 192*x5 - 18*M_PI - 1) + x7*(96*x3 - 128*x4 + 384*x5 - 24*M_PI - 4))
  };
}

Vector2 exact_xz(double x, double z, double B, double rhog, double s0, double alpha, double H) {
  return {
    2*H*alpha*rhog*x - 4*pow(alpha, 3)*pow(rhog, 3)*pow(x, 3)*(-pow(H, 4) + pow(-alpha*pow(x, 2) + s0 - z, 4))/pow(B, 3),
    0
  };
}

Vector2 source_xz(double x, double z, double B, double rhog, double s0, double alpha, double H) {
  double x0 = pow(x, 6);
  double x1 = pow(alpha, 6);
  double x2 = pow(x, 2);
  double x3 = -alpha*x2 + s0 - z;
  double x4 = pow(rhog, 6)/pow(B, 6);
  double x5 = x1*pow(x3, 6)*x4;
  double x6 = H*alpha*rhog;
  double x7 = pow(x, 4);
  double x8 = pow(alpha, 4);
  double x9 = pow(x3, 3);
  double x10 = pow(rhog, 3);
  double x11 = x10/pow(B, 3);
  double x12 = x11*x9;
  double x13 = x12*x8;
  double x14 = x13*x7;
  double x15 = pow(alpha, 3);
  double x16 = x11*x15*(-pow(H, 4) + pow(x3, 4));
  double x17 = x16*x2;
  double x18 = 32*x14 - 12*x17 + 2*x6;
  double x19 = 64*x0*x5 + pow(x18, 2);
  double x20 = pow(x19, -1.0/3.0);
  double x21 = pow(x3, 2);
  double x22 = pow(x, 3);
  double x23 = x10*x15*x22/pow(B, 2);
  double x24 = x13*x22;
  double x25 = pow(x, 5);
  double x26 = x11*x21;
  double x27 = pow(alpha, 5)*x25*x26;
  double x28 = x*x16;
  double x29 = (1.0/2.0)*B;
  double x30 = pow(x19, -4.0/3.0);
  double x31 = pow(x3, 5)*x4;
  double x32 = (1.0/3.0)*x18;
  return {
    -24*x20*x21*x23 + x20*x29*(896*x24 - 768*x27 - 96*x28) + 8*x23*x30*x9*(128*x0*x1*x31 - x32*(96*x12*x15*x2 - 192*x26*x7*x8)) + x29*x30*(128*x14 - 48*x17 + 8*x6)*(256*pow(alpha, 7)*pow(x, 7)*x31 - 128*x25*x5 - x32*(448*x24 - 384*x27 - 48*x28)),
    0
  };
}
} // end of namespace pism
