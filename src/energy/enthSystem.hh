// Copyright (C) 2009-2011, 2013, 2014, 2015, 2016, 2017, 2018 Andreas Aschwanden and Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __enthSystem_hh
#define __enthSystem_hh

#include <vector>

#include "pism/util/ColumnSystem.hh"
#include "pism/util/EnthalpyConverter.hh"

namespace pism {

class Config;
class IceModelVec3;

namespace energy {

//! Tridiagonal linear system for conservation of energy in vertical column of ice enthalpy.
/*!
  See the page documenting \ref bombproofenth.  The top of
  the ice has a Dirichlet condition.
*/
class enthSystemCtx : public columnSystemCtx {

public:
  enthSystemCtx(const std::vector<double>& storage_grid,
                const std::string &prefix,
                double dx,  double dy, double dt,
                const Config &config,
                const IceModelVec3 &Enth3,
                const IceModelVec3 &u3,
                const IceModelVec3 &v3,
                const IceModelVec3 &w3,
                const IceModelVec3 &strain_heating3,
                EnthalpyConverter::Ptr EC);
  ~enthSystemCtx();

  void init(int i, int j, bool ismarginal, double ice_thickness);

  double k_from_T(double T) const;

  void set_surface_heat_flux(double hf);
  void set_surface_neumann_bc(double dE);
  void set_surface_dirichlet_bc(double E_surface);

  void set_basal_dirichlet_bc(double E_basal);
  void set_basal_heat_flux(double hf);
  void set_basal_neumann_bc(double dE);

  virtual void save_system(std::ostream &output, unsigned int M) const;

  void solve(std::vector<double> &result);

  double lambda() const {
    return m_lambda;
  }

  double Enth(size_t i) const {
    return m_Enth[i];
  }

  double Enth_s(size_t i) const {
    return m_Enth_s[i];
  }
protected:
  // enthalpy in ice at previous time step
  std::vector<double> m_Enth;
  // enthalpy level for CTS; function only of pressure
  std::vector<double> m_Enth_s;

  // temporary storage for ice enthalpy at (i,j), as well as north,
  // east, south, and west from (i,j)
  std::vector<double> m_E_ij, m_E_n, m_E_e, m_E_s, m_E_w;

  //! strain heating in the ice column
  std::vector<double> m_strain_heating;

  //! values of @f$ k \Delta t / (\rho c \Delta x^2) @f$
  std::vector<double> m_R;

  double m_ice_density, m_ice_c, m_ice_k, m_p_air,
    m_nu, m_R_cold, m_R_temp, m_R_factor;

  double m_ice_thickness,
    m_lambda;              //!< implicit FD method parameter
  double m_D0, m_U0, m_B0;   // coefficients of the first (basal) equation
  double m_L_ks, m_D_ks, m_U_ks, m_B_ks;   // coefficients of the last (surface) equation
  bool m_marginal, m_k_depends_on_T;

  bool m_exclude_horizontal_advection;
  bool m_exclude_vertical_advection;
  bool m_exclude_strain_heat;

  const IceModelVec3 &m_Enth3, &m_strain_heating3;
  EnthalpyConverter::Ptr m_EC;  // conductivity has known dependence on T, not enthalpy

  void compute_enthalpy_CTS();
  double compute_lambda();

  void assemble_R();
  void checkReadyToSolve();
};

} // end of namespace energy
} // end of namespace pism

#endif   //  ifndef __enthSystem_hh
