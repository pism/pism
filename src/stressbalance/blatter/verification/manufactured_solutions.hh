// DO NOT EDIT. This code was generated by a Python script.
#include "pism/util/Vector2d.hh"

namespace pism {

Vector2d blatter_xy_exact(double x, double y);

Vector2d blatter_xy_source(double x, double y, double B);

Vector2d blatter_xz_exact(double x, double z, double A, double rho, double g, double s_0, double alpha, double H, double beta);

Vector2d blatter_xz_source(double x, double z, double A, double rho, double g, double s_0, double alpha, double H, double beta);

Vector2d blatter_xz_source_bed(double x, double z, double A, double rho, double g, double s_0, double alpha, double H, double beta);

Vector2d blatter_xz_source_surface(double x, double z, double A, double rho, double g, double s_0, double alpha, double H, double beta);

Vector2d blatter_xz_cfbc_exact(double x, double z, double B, double L, double rho_i, double rho_w, double g);

Vector2d blatter_xz_cfbc_source(double x, double z, double L, double rho_i, double rho_w, double g);

Vector2d blatter_xz_cfbc_surface(double x, double L, double rho_i, double rho_w, double g);

Vector2d blatter_xz_cfbc_base(double x, double L, double rho_i, double rho_w, double g);

Vector2d blatter_xz_halfar_exact(double x, double z, double H_0, double R_0, double rho_i, double g, double B);

Vector2d blatter_xz_halfar_source(double x, double z, double H_0, double R_0, double rho_i, double g, double B);

Vector2d blatter_xz_halfar_source_lateral(double x, double z, double H_0, double R_0, double rho_i, double g, double B);

Vector2d blatter_xz_halfar_source_surface(double x, double H_0, double R_0, double rho_i, double g, double B);

Vector2d blatter_xz_halfar_source_base(double x, double H_0, double R_0, double rho_i, double g, double B);

double blatter_xz_vanderveen_thickness(double x, double alpha, double H_0, double Q_0, double rho_i, double g, double B);

Vector2d blatter_xz_vanderveen_exact(double x, double alpha, double H_0, double Q_0, double rho_i, double g, double B);

Vector2d blatter_xz_vanderveen_source_lateral(double x, double alpha, double H_0, double Q_0, double rho_i, double g, double B);

Vector2d blatter_xz_vanderveen_source_surface(double x, double alpha, double H_0, double Q_0, double rho_i, double g, double B);

double blatter_xz_vanderveen_beta(double x, double alpha, double H_0, double Q_0, double rho_i, double g, double B);

} // end of namespace pism