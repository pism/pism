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

#include "PISMFEvoR.hh"
#include "PISMConfig.hh"
#include "PISMVars.hh"
#include "PISMStressBalance_diagnostics.hh"
#include "enthalpyConverter.hh"
#include <cassert>

namespace pism {

PISMFEvoR::PISMFEvoR(IceGrid &g, const Config &conf, EnthalpyConverter *EC,
                     StressBalance *stress_balance)
  : Component_TS(g, conf), m_stress_balance(stress_balance), m_EC(EC) {

  assert(m_EC != NULL);
  assert(m_stress_balance != NULL);

  PetscErrorCode ierr = allocate(); CHKERRCONTINUE(ierr);

  // will be allocated in init()
  m_pressure = NULL;
  m_tauxz    = NULL;
  m_tauyz    = NULL;

  // points to storage owned by IceModel
  m_enthalpy = NULL;
}

PISMFEvoR::~PISMFEvoR() {
  delete m_pressure;
  delete m_tauxz;
  delete m_tauyz;
}

PetscErrorCode PISMFEvoR::max_timestep(double t, double &dt, bool &restrict) {
  // FIXME: put real code here
  PetscErrorCode ierr = Component_TS::max_timestep(t, dt, restrict); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMFEvoR::update(double t, double dt) {
  m_t = t;
  m_dt = dt;
  PetscErrorCode ierr;

  IceModelVec3 *u = NULL, *v = NULL, *w = NULL;
  ierr = m_stress_balance->get_3d_velocity(u, v, w); CHKERRQ(ierr);

  assert(m_pressure != NULL);
  IceModelVec* pressure = NULL; // FIXME: use a smart pointer
  ierr = m_pressure->compute(pressure); CHKERRQ(ierr);
  IceModelVec3 *pressure3 = static_cast<IceModelVec3*>(pressure);

  assert(m_tauxz != NULL);
  IceModelVec* tauxz = NULL;    // FIXME: use a smart pointer
  ierr = m_tauxz->compute(tauxz); CHKERRQ(ierr);
  IceModelVec3 *tauxz3 = static_cast<IceModelVec3*>(tauxz);

  assert(m_tauyz != NULL);
  IceModelVec* tauyz = NULL;    // FIXME: use a smart pointer
  ierr = m_tauyz->compute(tauyz); CHKERRQ(ierr);
  IceModelVec3 *tauyz3 = static_cast<IceModelVec3*>(tauyz);

  {
    ierr = m_enhancement_factor.set(config.get("sia_enhancement_factor")); CHKERRQ(ierr);

    unsigned int n_particles = 1;
        /* terminology:
         *  particles exist in pism and contain one or more distrobutions 
         * of crystals that are tracked through time. They essentially are 
         * infinitesimely small. Distributions exist in FEvoR and contain 
         * sets of independent crystals (or in the case of NNI weakly 
         * dependant). In PISM-FEvoR you will likely not need to access the
         * crystals directly. Methods should be provided through FEvoR's 
         * distribution class fevor_distribution. 
         */
        
    for (unsigned int i = 0; i < n_particles; ++i) {
      double x = 0.0,
             y = 0.0,
             z = 0.0;
             
      double P   = 0.0, 
             txz = 0.0,
             tyz = 0.0;
             
      double E = 90e3, 
             T = 0.0;

      // check if the point (x,y,z) is within the domain:
      assert(0.0 <= z && z <= grid.Lz);
      assert(-grid.Lx <= x && x <= grid.Lx);
      assert(-grid.Ly <= y && y <= grid.Ly);

      ierr = interp_field_point( x, y, z, pressure3, P  );
      ierr = interp_field_point( x, y, z, tauxz3   , txy);
      ierr = interp_field_point( x, y, z, tauyz3   , txz);
      
      
      ierr = m_EC->getAbsTemp(E, P, T); CHKERRQ(ierr);
      
      // TODO put in FEvoR step.
    }
  }

  delete pressure;
  delete tauxz;
  delete tauyz;

  return 0;
}

virtual PetscErrorCode PISMFEvoR::interp_field_point( double &x, double &y, 
                                                      double &z, 
                                                      IceModelVec3 &field3, 
                                                      double &feildValue) {
    PetscErrorCode ierr;
    
    int I = 0, J = 0;
    grid.compute_point_neighbors(x, y, I, J);
    std::vector<double> weights = grid.compute_interp_weights(x, y);

    double *column0 = NULL, *column1 = NULL, *column2 = NULL, *column3 = NULL;
    ierr = field3->getInternalColumn(I,   J,   &column0); CHKERRQ(ierr);
    ierr = field3->getInternalColumn(I+1, J,   &column1); CHKERRQ(ierr);
    ierr = field3->getInternalColumn(I+1, J+1, &column2); CHKERRQ(ierr);
    ierr = field3->getInternalColumn(I,   J+1, &column3); CHKERRQ(ierr);

    unsigned int K = 0;
    // K + 1 (used below) should be at most Mz - 1
    while (K + 1 < grid.Mz - 1 && grid.zlevels[K + 1] < z) {
    K++;
    }

    double z_weight = (z - grid.zlevels[K]) / (grid.zlevels[K+1] - grid.zlevels[K]);

    double f0 = column0[K] + z_weight * (column0[K+1] - column0[K]);
    double f1 = column1[K] + z_weight * (column1[K+1] - column1[K]);
    double f2 = column2[K] + z_weight * (column2[K+1] - column2[K]);
    double f3 = column3[K] + z_weight * (column3[K+1] - column3[K]);

    feildValue = weights[0] * f0 + weights[1] * f1 + weights[2] * f2 + weights[3] * f3;
    
    return 0;
}

void PISMFEvoR::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  if (keyword != "none") {
    result.insert(m_enhancement_factor.metadata().get_string("short_name"));
  }
}

PetscErrorCode PISMFEvoR::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                           IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "enhancement_factor")) {
    ierr = m_enhancement_factor.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMFEvoR::write_variables(const std::set<std::string> &vars, const PIO& nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "enhancement_factor")) {
    ierr = m_enhancement_factor.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMFEvoR::allocate() {
  PetscErrorCode ierr;

  // SIAFD diffusive flux computation requires stencil width of 1
  const unsigned int stencil_width = 1;

  ierr = m_enhancement_factor.create(grid, "enhancement_factor", WITH_GHOSTS,
                                     stencil_width); CHKERRQ(ierr);
  ierr = m_enhancement_factor.set_attrs("diagnostic", // i.e. not needed to re-start the model
                                        "flow law enhancement factor",
                                        "1", // dimensionless
                                        "" /* no standard name */); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMFEvoR::init(Vars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
                    "* Initializing the Fabric Evolution with Recrystallization"
                    " (FEvoR) model...\n"); CHKERRQ(ierr);

  // make enhancement factor available to other PISM components:
  ierr = vars.add(m_enhancement_factor); CHKERRQ(ierr);

  m_enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (m_enthalpy == NULL) {
    SETERRQ(grid.com, 1, "enthalpy field is not available");
  }

  // It would be nice to be able to allocate these in
  // PISMFEvoR::allocate() or in the constructor, but pism::Vars is
  // not available there...
  if (m_pressure == NULL) {
    m_pressure = new PSB_pressure(m_stress_balance, grid, vars);
    assert(m_pressure != NULL);
  }

  if (m_tauxz == NULL) {
    m_tauxz = new PSB_tauxz(m_stress_balance, grid, vars);
    assert(m_tauxz != NULL);
  }

  if (m_tauyz == NULL) {
    m_tauyz = new PSB_tauyz(m_stress_balance, grid, vars);
    assert(m_tauyz != NULL);
  }

  return 0;
}

} // end of namespace pism
