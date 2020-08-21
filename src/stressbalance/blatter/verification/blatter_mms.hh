#include <cmath>

#include "pism/util/Vector2.hh"

namespace pism {

Vector2 exact_xy(double x, double y);

Vector2 source_xy(double x, double y, double B);

Vector2 exact_xz(double x, double z, double B, double rhog, double s0, double alpha, double H);

Vector2 source_xz(double x, double z, double B, double rhog, double s0, double alpha, double H);

} // end of namespace pism
