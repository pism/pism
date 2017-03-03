// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2017 Constantine Khroulev
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

#include "PBLingleClark.hh"

#include <gsl/gsl_math.h>       // GSL_NAN

#include "base/util/io/PIO.hh"
#include "base/util/PISMTime.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_options.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/PISMVars.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"
#include "deformation.hh"

namespace pism {
namespace bed {

PBLingleClark::PBLingleClark(IceGrid::ConstPtr g)
  : BedDef(g) {

  m_work0 = m_topg.allocate_proc0_copy();

  // A work vector. This storage is used to put thickness change on rank 0 and to get the plate
  // displacement change back.
  m_work.create(m_grid, "work_vector", WITHOUT_GHOSTS);
  m_work.set_attrs("internal", "a work vector", "meters", "");

  m_topg_start.create(m_grid, "topg_start", WITHOUT_GHOSTS);
  m_topg_start.set_attrs("internal",
                         "bed elevation at the beginning of the run",
                         "meters", "");

  bool use_elastic_model = m_config->get_boolean("bed_deformation.lc.elastic_model");

  m_bdLC = NULL;

  const int
    Mx = m_grid->Mx(),
    My = m_grid->My(),
    Z  = 4,                     // use Z = 4 for now; to reduce global drift?
    Nx = Z*(Mx - 1) + 1,
    Ny = Z*(My - 1) + 1;

  const double
    Lx = Z * (m_grid->x0() - m_grid->x(0)),
    Ly = Z * (m_grid->y0() - m_grid->y(0));

  m_extended_grid = IceGrid::Shallow(m_grid->ctx(),
                                     Lx, Ly,
                                     m_grid->x0(), m_grid->y0(),
                                     Nx, Ny, NOT_PERIODIC);

  m_plate_displacement.create(m_extended_grid, "plate_displacement", WITHOUT_GHOSTS);
  m_plate_displacement.set_attrs("model state",
                                 "viscous plate displacement in the "
                                 "Lingle-Clark bed deformation model", "meters", "");
  m_plate_displacement.metadata().set_string("comment",
                                             "This field has no physical meaning"
                                             " and should not be used to interpret model results.");
  // coordinate variables of the extended grid should have different names
  m_plate_displacement.metadata().get_x().set_name("x_lc");
  m_plate_displacement.metadata().get_y().set_name("y_lc");

  // Set up scatters to and from processor 0 by allocating a processor 0 copy. (This copy is thrown
  // away immediately.)
  m_work0_extended = m_plate_displacement.allocate_proc0_copy();

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      m_bdLC = new BedDeformLC(*m_config, use_elastic_model,
                               Mx, My,
                               m_grid->dx(), m_grid->dy(),
                               Nx, Ny);
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}

PBLingleClark::~PBLingleClark() {
  if (m_bdLC != NULL) {
    delete m_bdLC;
    m_bdLC = NULL;
  }
}

void PBLingleClark::uplift_problem(const IceModelVec2S& ice_thickness,
                                   const IceModelVec2S& bed_uplift) {

  // Note: this method is not called every very often so it's OK to allocate a rank 0 copy of
  // bed_uplift and free it at the end of scope.

  petsc::Vec::Ptr
    thickness = m_work0,
    uplift    = bed_uplift.allocate_proc0_copy();

  ice_thickness.put_on_proc0(*thickness);
  bed_uplift.put_on_proc0(*uplift);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      m_bdLC->uplift_problem(*thickness, *uplift);
      PetscErrorCode ierr = VecCopy(m_bdLC->plate_displacement_change(), *m_work0);
      PISM_CHK(ierr, "VecCopy");
      ierr = VecCopy(m_bdLC->plate_displacement(), *m_work0_extended);
      PISM_CHK(ierr, "VecCopy");
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  m_topg.get_from_proc0(*m_work0);
  m_plate_displacement.get_from_proc0(*m_work0_extended);
}

void PBLingleClark::bootstrap_impl(const IceModelVec2S &bed,
                                   const IceModelVec2S &bed_uplift,
                                   const IceModelVec2S &ice_thickness) {
  m_t_beddef_last = m_grid->ctx()->time()->start();

  m_t  = GSL_NAN;
  m_dt = GSL_NAN;

  m_topg_start.copy_from(bed);
  m_topg_last.copy_from(m_topg);

  // initialize the plate displacement
  {
    bed_uplift.put_on_proc0(*m_work0);

    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {
        m_bdLC->bootstrap(*m_work0);
      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();
  }
}

/*! Initialize the Lingle-Clark bed deformation model using uplift.
 *
 * Inputs:
 *
 * - bed topography,
 * - ice thickness,
 * - plate displacement (either read from a file or bootstrapped using uplift) and
 *   possibly re-gridded.
 */
void PBLingleClark::init_impl(const InputOptions &opts) {
  m_log->message(2, "* Initializing the Lingle-Clark bed deformation model...\n");

  // Initialize bed topography and uplift maps.
  BedDef::init_impl(opts);

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  m_topg_start.copy_from(m_topg);
  m_topg_last.copy_from(m_topg);

  if (opts.type == INIT_RESTART) {
    // Set m_plate_displacement by reading from the input file.
    m_plate_displacement.read(opts.filename, opts.record);
  } else if (opts.type == INIT_BOOTSTRAP) {
    // Set m_plate_displacement by solving the "uplift problem".
    this->uplift_problem(ice_thickness, m_uplift);
    // Re-set m_topg because uplift_problem() modified it.
    m_topg.copy_from(m_topg_start);
  } else {
    // do nothing
  }

  // Try re-gridding plate_displacement.
  regrid("Lingle-Clark bed deformation model", m_plate_displacement, REGRID_WITHOUT_REGRID_VARS);

  // Now that m_plate_displacement is finally initialized, put it on rank 0 and initialize m_bdLC
  // itself.
  m_plate_displacement.put_on_proc0(*m_work0_extended);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {  // only processor zero does the step
      m_bdLC->init(*m_work0_extended);
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}

MaxTimestep PBLingleClark::max_timestep_impl(double t) const {
  (void) t;
  // no time-step restriction
  return MaxTimestep("bed_def lc");
}

//! Update the Lingle-Clark bed deformation model.
void PBLingleClark::update_with_thickness_impl(const IceModelVec2S &ice_thickness,
                                               double t, double dt) {

  if ((fabs(t - m_t)   < 1e-12) &&
      (fabs(dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = t;
  m_dt = dt;

  double t_final = m_t + m_dt;

  // Check if it's time to update:
  double dt_beddef = t_final - m_t_beddef_last; // in seconds
  if ((dt_beddef < m_config->get_double("bed_deformation.update_interval", "seconds") and
       t_final < m_grid->ctx()->time()->end()) or
      dt_beddef < 1e-12) {
    return;
  }

  m_t_beddef_last = t_final;

  ice_thickness.put_on_proc0(*m_work0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {  // only processor zero does the step
      m_bdLC->step(dt_beddef, *m_work0);
      PetscErrorCode ierr = VecCopy(m_bdLC->plate_displacement_change(), *m_work0);
      PISM_CHK(ierr, "VecCopy");
      ierr = VecCopy(m_bdLC->plate_displacement(), *m_work0_extended);
      PISM_CHK(ierr, "VecCopy");
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  // Update bed topography using the plate displacement change obtained from rank 0.
  {
    IceModelVec2S &dU = m_work;

    dU.get_from_proc0(*m_work0);

    IceModelVec::AccessList list{&m_topg, &m_topg_start, &dU};

    ParallelSection loop(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        m_topg(i, j) = m_topg_start(i, j) + dU(i, j);
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();
  }

  //! Finally, we need to update bed uplift and topg_last.
  compute_uplift(dt_beddef);
  m_topg_last.copy_from(m_topg);

  //! Increment the topg state counter. SIAFD relies on this!
  m_topg.inc_state_counter();

  // get plate displacement on the extended grid from processor 0
  m_plate_displacement.get_from_proc0(*m_work0_extended);
}

void PBLingleClark::define_model_state_impl(const PIO &output) const {
  BedDef::define_model_state_impl(output);
  m_plate_displacement.define(output);
}

void PBLingleClark::write_model_state_impl(const PIO &output) const {
  BedDef::write_model_state_impl(output);
  m_plate_displacement.write(output);
}

} // end of namespace bed
} // end of namespace pism
