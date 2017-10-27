// Copyright (C) 2009-2011, 2013, 2014, 2015, 2017 Ed Bueler
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

#ifndef __tempSystem_hh
#define __tempSystem_hh

#include "pism/util/ColumnSystem.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Mask.hh"

namespace pism {

class IceModelVec3;

namespace energy {
//! Tridiagonal linear system for vertical column of temperature-based conservation of energy problem.
/*!
  Call sequence like this:
  \code
  tempSystemCtx foo;
  foo.dx = ...  // set public constants
  foo.u = ...   // set public pointers
  foo.initAllColumns();
  for (j in ownership) {
  for (i in ownership) {
  ks = ...
  foo.setIndicesThisColumn(i,j,ks);
  [COMPUTE OTHER PARAMS]
  foo.setSchemeParamsThisColumn(mask,isMarginal,lambda);  
  foo.setSurfaceBoundaryValuesThisColumn(Ts);
  foo.setBasalBoundaryValuesThisColumn(Ghf,Tshelfbase,Rb);
  foo.solveThisColumn(x);
  }  
  }
  \endcode
*/
class tempSystemCtx : public columnSystemCtx {
public:
  tempSystemCtx(const std::vector<double>& storage_grid,
                const std::string &prefix,
                double dx, double dy, double dt,
                const Config &config,
                const IceModelVec3 &T3,
                const IceModelVec3 &u3,
                const IceModelVec3 &v3,
                const IceModelVec3 &w3,
                const IceModelVec3 &strain_heating3);

  void initThisColumn(int i, int j, bool is_marginal, MaskValue new_mask, double ice_thickness);

  void setSurfaceBoundaryValuesThisColumn(double my_Ts);
  void setBasalBoundaryValuesThisColumn(double my_G0, double my_Tshelfbase,
                                                  double my_Rb);

  void solveThisColumn(std::vector<double> &x);

  double lambda() {
    return m_lambda;
  }

  double w(int k) {
    return m_w[k];
  }
protected:
  double m_ice_density, m_ice_c, m_ice_k;
  const IceModelVec3 &m_T3, &m_strain_heating3;

  std::vector<double>  m_T, m_strain_heating;
  std::vector<double> m_T_n, m_T_e, m_T_s, m_T_w;

  double m_lambda, m_Ts, m_G0, m_Tshelfbase, m_Rb;
  MaskValue    m_mask;
  bool        m_is_marginal;
  double m_nu,
    m_rho_c_I,
    m_iceK,
    m_iceR;
private:
  bool
    m_surfBCsValid,
    m_basalBCsValid;

  double compute_lambda();
};

} // end of namespace energy
} // end of namespace pism

#endif  /* __tempSystem_hh */

