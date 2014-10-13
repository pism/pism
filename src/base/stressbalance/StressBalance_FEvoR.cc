/* Copyright (C) 2014 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "StressBalance_FEvoR.hh"
#include "PISMVars.hh"
#include "flowlaws.hh"
#include "ShallowStressBalance.hh"
#include "SSB_Modifier.hh"
#include "Mask.hh"
#include "enthalpyConverter.hh"

namespace pism {

StressBalance_FEvoR::StressBalance_FEvoR(IceGrid &g, ShallowStressBalance *sb, SSB_Modifier *ssb_mod,
                                         const Config &conf)
  : StressBalance(g, sb, ssb_mod, conf) {
  // empty
}

StressBalance_FEvoR::~StressBalance_FEvoR() {
  // empty
}

// FIXME: Code duplication.
static inline double D2(double u_x, double u_y, double u_z, double v_x, double v_y, double v_z) {
  return 0.5 * (PetscSqr(u_x + v_y) + u_x*u_x + v_y*v_y + 0.5 * (PetscSqr(u_y + v_x) + u_z*u_z + v_z*v_z));
}

PetscErrorCode StressBalance_FEvoR::compute_volumetric_strain_heating() {
  PetscErrorCode ierr;
  IceModelVec3 *u, *v, *enthalpy;
  IceModelVec2S *thickness;
  const IceFlowLaw *flow_law = m_stress_balance->get_flow_law();
  EnthalpyConverter &EC = m_stress_balance->get_enthalpy_converter();

  ierr = m_modifier->get_horizontal_3d_velocity(u, v); CHKERRQ(ierr);

  IceModelVec2Int *mask = dynamic_cast<IceModelVec2Int*>(m_variables->get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  MaskQuery m(*mask);

  thickness = dynamic_cast<IceModelVec2S*>(m_variables->get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(m_variables->get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(grid.com, 1, "enthalpy is not available");

  double
    n = flow_law->exponent(),
    exponent = 0.5 * (1.0 / n + 1.0);

  IceModelVec::AccessList list;
  // new code
  IceModelVec3* enhancement_factor = dynamic_cast<IceModelVec3*>(m_variables->get("enhancement_factor"));
  if (enhancement_factor == NULL) SETERRQ(grid.com, 1, "enhancement_factor is not available");

  double *EF_ij;
  list.add(*enhancement_factor);
  // end of new code

  list.add(*mask);
  list.add(*enthalpy);
  list.add(m_strain_heating);
  list.add(*thickness);
  list.add(*u);
  list.add(*v);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double H = (*thickness)(i,j);
    int ks = grid.kBelowHeight(H);
    double
      *u_ij, *u_w, *u_n, *u_e, *u_s,
      *v_ij, *v_w, *v_n, *v_e, *v_s;
    double *Sigma, *E_ij;

    double west = 1, east = 1, south = 1, north = 1,
      D_x = 0,                // 1/(dx), 1/(2dx), or 0
      D_y = 0;                // 1/(dy), 1/(2dy), or 0

    // x-derivative
    {
      if ((m.icy(i,j) && m.ice_free(i+1,j)) || (m.ice_free(i,j) && m.icy(i+1,j)))
        east = 0;
      if ((m.icy(i,j) && m.ice_free(i-1,j)) || (m.ice_free(i,j) && m.icy(i-1,j)))
        west = 0;

      if (east + west > 0)
        D_x = 1.0 / (grid.dx * (east + west));
      else
        D_x = 0.0;
    }

    // y-derivative
    {
      if ((m.icy(i,j) && m.ice_free(i,j+1)) || (m.ice_free(i,j) && m.icy(i,j+1)))
        north = 0;
      if ((m.icy(i,j) && m.ice_free(i,j-1)) || (m.ice_free(i,j) && m.icy(i,j-1)))
        south = 0;

      if (north + south > 0)
        D_y = 1.0 / (grid.dy * (north + south));
      else
        D_y = 0.0;
    }


    ierr = u->getInternalColumn(i,     j,     &u_ij); CHKERRQ(ierr);
    ierr = u->getInternalColumn(i - 1, j,     &u_w);  CHKERRQ(ierr);
    ierr = u->getInternalColumn(i + 1, j,     &u_e);  CHKERRQ(ierr);
    ierr = u->getInternalColumn(i,     j - 1, &u_s);  CHKERRQ(ierr);
    ierr = u->getInternalColumn(i,     j + 1, &u_n);  CHKERRQ(ierr);

    ierr = v->getInternalColumn(i,     j,     &v_ij); CHKERRQ(ierr);
    ierr = v->getInternalColumn(i - 1, j,     &v_w);  CHKERRQ(ierr);
    ierr = v->getInternalColumn(i + 1, j,     &v_e);  CHKERRQ(ierr);
    ierr = v->getInternalColumn(i,     j - 1, &v_s);  CHKERRQ(ierr);
    ierr = v->getInternalColumn(i,     j + 1, &v_n);  CHKERRQ(ierr);

    ierr =        enthalpy->getInternalColumn(i, j, &E_ij);  CHKERRQ(ierr);
    ierr = m_strain_heating.getInternalColumn(i, j, &Sigma); CHKERRQ(ierr);

    // new code
    ierr = enhancement_factor->getInternalColumn(i, j, &EF_ij); CHKERRQ(ierr);
    // end of new code

    for (int k = 0; k <= ks; ++k) {
      double dz,
        pressure = EC.getPressureFromDepth(H - grid.zlevels[k]),
        B        = flow_law->hardness_parameter(E_ij[k], pressure);

      double u_z = 0.0, v_z = 0.0,
        u_x = D_x * (west  * (u_ij[k] - u_w[k]) + east  * (u_e[k] - u_ij[k])),
        u_y = D_y * (south * (u_ij[k] - u_s[k]) + north * (u_n[k] - u_ij[k])),
        v_x = D_x * (west  * (v_ij[k] - v_w[k]) + east  * (v_e[k] - v_ij[k])),
        v_y = D_y * (south * (v_ij[k] - v_s[k]) + north * (v_n[k] - v_ij[k]));

      if (k > 0) {
        dz = grid.zlevels[k+1] - grid.zlevels[k-1];
        u_z = (u_ij[k+1] - u_ij[k-1]) / dz;
        v_z = (v_ij[k+1] - v_ij[k-1]) / dz;
      } else {
        // use one-sided differences for u_z and v_z on the bottom level
        dz = grid.zlevels[1] - grid.zlevels[0];
        u_z = (u_ij[1] - u_ij[0]) / dz;
        v_z = (v_ij[1] - v_ij[0]) / dz;
      }

      // new code
      Sigma[k] = 2.0 * pow(EF_ij[k], -1.0/n) * B * pow(D2(u_x, u_y, u_z, v_x, v_y, v_z), exponent);
      // end of new code
    } // k-loop

    int remaining_levels = grid.Mz - (ks + 1);
    if (remaining_levels > 0) {
      ierr = PetscMemzero(&Sigma[ks+1],
                          remaining_levels*sizeof(double)); CHKERRQ(ierr);
    }
  }

  return 0;
}

} // end of namespace pism
