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

#include "PISMStressBalance.hh"
#include "ShallowStressBalance.hh"
#include "SSB_Modifier.hh"
#include "coupler/PISMOcean.hh"
#include "base/enthalpyConverter.hh"
#include "base/rheology/flowlaws.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/Profiling.hh"
#include "base/stressbalance/StressBalance_FEvoR.hh"

namespace pism {
  namespace stressbalance{
    StressBalance_FEvoR::StressBalance_FEvoR(IceGrid::ConstPtr g, ShallowStressBalance *sb, SSB_Modifier *ssb_mod)
      : StressBalance(g, sb, ssb_mod) {
       
  // empty
}

StressBalance_FEvoR::~StressBalance_FEvoR() {
  // empty
}

// FIXME: Code duplication.
static inline double D2(double u_x, double u_y, double u_z, double v_x, double v_y, double v_z) {
  return 0.5 * (PetscSqr(u_x + v_y) + u_x*u_x + v_y*v_y + 0.5 * (PetscSqr(u_y + v_x) + u_z*u_z + v_z*v_z));
}

void StressBalance_FEvoR::compute_volumetric_strain_heating() {
  PetscErrorCode ierr;
  m_log->message(4,
                    "\n Beginning compute_volumetric_strain_heating()\n"); 

  const rheology::FlowLaw *flow_law = m_shallow_stress_balance->flow_law();
  EnthalpyConverter::Ptr EC = m_shallow_stress_balance->enthalpy_converter();

  const IceModelVec3
    &u = m_modifier->velocity_u(),
    &v = m_modifier->velocity_v();

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec3  *enthalpy  = m_grid->variables().get_3d_scalar("enthalpy");

  const IceModelVec2Int *mask = m_grid->variables().get_2d_mask("mask");
  MaskQuery m(*mask);

  double
    n = flow_law->exponent(),
    exponent = 0.5 * (1.0 / n + 1.0);

  IceModelVec::AccessList list;
  // new codes
  const IceModelVec3* enhancement_factor = m_grid->variables().get_3d_scalar("enhancement_factor");
  if (enhancement_factor == NULL) throw RuntimeError( "enhancement_factor is not available");

    const double *EF_ij;
  list.add(*enhancement_factor);
  // end of new code


  list.add(*mask);
  list.add(*enthalpy);
  list.add(m_strain_heating);
  list.add(*thickness);
  list.add(u);
  list.add(v);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double H = (*thickness)(i,j);
      int ks = m_grid->kBelowHeight(H);
      const double
        *u_ij, *u_w, *u_n, *u_e, *u_s,
        *v_ij, *v_w, *v_n, *v_e, *v_s;
      double *Sigma;
      const double *E_ij;

      double west = 1, east = 1, south = 1, north = 1,
        D_x = 0,                // 1/(dx), 1/(2dx), or 0
        D_y = 0;                // 1/(dy), 1/(2dy), or 0

      // x-derivative
      {
        if ((m.icy(i,j) && m.ice_free(i+1,j)) || (m.ice_free(i,j) && m.icy(i+1,j))) {
          east = 0;
        }
        if ((m.icy(i,j) && m.ice_free(i-1,j)) || (m.ice_free(i,j) && m.icy(i-1,j))) {
          west = 0;
        }

        if (east + west > 0) {
          D_x = 1.0 / (m_grid->dx() * (east + west));
        } else {
          D_x = 0.0;
        }
      }

      // y-derivative
      {
        if ((m.icy(i,j) && m.ice_free(i,j+1)) || (m.ice_free(i,j) && m.icy(i,j+1))) {
          north = 0;
        }
        if ((m.icy(i,j) && m.ice_free(i,j-1)) || (m.ice_free(i,j) && m.icy(i,j-1))) {
          south = 0;
        }

        if (north + south > 0) {
          D_y = 1.0 / (m_grid->dy() * (north + south));
        } else {
          D_y = 0.0;
        }
      }

      u_ij = u.get_column(i,     j);
      u_w  = u.get_column(i - 1, j);
      u_e  = u.get_column(i + 1, j);
      u_s  = u.get_column(i,     j - 1);
      u_n  = u.get_column(i,     j + 1);

      v_ij = v.get_column(i,     j);
      v_w  = v.get_column(i - 1, j);
      v_e  = v.get_column(i + 1, j);
      v_s  = v.get_column(i,     j - 1);
      v_n  = v.get_column(i,     j + 1);

      E_ij = enthalpy->get_column(i, j);
      Sigma = m_strain_heating.get_column(i, j);

    // new code
    EF_ij = enhancement_factor->get_column(i, j);
    // end of new code

    for (int k = 0; k <= ks; ++k) {
        double dz,
          pressure = EC->pressure(H - m_grid->z(k)),
          B        = flow_law->hardness_parameter(E_ij[k], pressure);

        double u_z = 0.0, v_z = 0.0,
          u_x = D_x * (west  * (u_ij[k] - u_w[k]) + east  * (u_e[k] - u_ij[k])),
          u_y = D_y * (south * (u_ij[k] - u_s[k]) + north * (u_n[k] - u_ij[k])),
          v_x = D_x * (west  * (v_ij[k] - v_w[k]) + east  * (v_e[k] - v_ij[k])),
          v_y = D_y * (south * (v_ij[k] - v_s[k]) + north * (v_n[k] - v_ij[k]));

        if (k > 0) {
          dz = m_grid->z(k+1) - m_grid->z(k-1);
          u_z = (u_ij[k+1] - u_ij[k-1]) / dz;
          v_z = (v_ij[k+1] - v_ij[k-1]) / dz;
        } else {
          // use one-sided differences for u_z and v_z on the bottom level
          dz = m_grid->z(1) - m_grid->z(0);
          u_z = (u_ij[1] - u_ij[0]) / dz;
          v_z = (v_ij[1] - v_ij[0]) / dz;
        }


      // new code
      Sigma[k] = 2.0 * pow(EF_ij[k], -1.0/n) * B * pow(D2(u_x, u_y, u_z, v_x, v_y, v_z), exponent);
      // end of new code
    } // k-loop


      int remaining_levels = m_grid->Mz() - (ks + 1);
      if (remaining_levels > 0) {
        ierr = PetscMemzero(&Sigma[ks+1],
                            remaining_levels*sizeof(double));
        PISM_CHK(ierr, "PetscMemzero");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
  m_log->message(4,
                    "\n Done with compute_volumetric_strain_heating()\n"); 

}
  } // end of namespace stressbalance
} // end of namespace pism
