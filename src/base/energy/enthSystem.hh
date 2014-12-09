// Copyright (C) 2009-2011, 2013, 2014 Andreas Aschwanden and Ed Bueler
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

#include "columnSystem.hh"

namespace pism {

class Config;
class IceModelVec3;
class EnthalpyConverter;

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
                IceModelVec3 &Enth3,
                IceModelVec3 *u3,
                IceModelVec3 *v3,
                IceModelVec3 *w3,
                IceModelVec3 *strain_heating3,
                const EnthalpyConverter &EC);
  ~enthSystemCtx();

  void initThisColumn(int i, int j, bool my_ismarginal,
                      double ice_thickness);

  double k_from_T(double T);

  void setDirichletSurface(double my_Enth_surface);
  void setDirichletBasal(double Y);
  void setBasalHeatFlux(double hf);

  virtual void save_system(std::ostream &output, unsigned int M) const;

  void solveThisColumn(std::vector<double> &result);

  double lambda() {
    return m_lambda;
  }

  double Enth(size_t i) {
    return m_Enth[i];
  }

  double Enth_s(size_t i) {
    return m_Enth_s[i];
  }
protected:
  // enthalpy in ice at previous time step
  std::vector<double> m_Enth;
  // enthalpy level for CTS; function only of pressure
  std::vector<double> m_Enth_s;

  //! strain heating in the ice column
  std::vector<double> m_strain_heating;

  //! values of @f$ k \Delta t / (\rho c \Delta x^2) @f$
  std::vector<double> m_R;

  double m_ice_density, m_ice_c, m_ice_k, m_ice_K, m_ice_K0, m_p_air,
    m_nu, m_R_cold, m_R_temp, m_R_factor;

  double m_ice_thickness,
    m_lambda,              //!< implicit FD method parameter
    m_Enth_ks;             //!< top surface Dirichlet B.C.
  double m_D0, m_U0, m_B0;   // coefficients of the first (basal) equation
  bool m_ismarginal, m_c_depends_on_T, m_k_depends_on_T;

  IceModelVec3 *m_Enth3, *m_strain_heating3;
  const EnthalpyConverter &m_EC;  // conductivity has known dependence on T, not enthalpy

  void compute_enthalpy_CTS();
  double compute_lambda();

  void assemble_R();
  void checkReadyToSolve();
};

} // end of namespace pism

#endif   //  ifndef __enthSystem_hh
