// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev
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

namespace pism {
namespace bed {

PBLingleClark::PBLingleClark(IceGrid::ConstPtr g)
  : BedDef(g) {

  m_Hp0        = m_topg_initial.allocate_proc0_copy();
  m_bedp0      = m_topg_initial.allocate_proc0_copy();
  m_Hstartp0   = m_topg_initial.allocate_proc0_copy();
  m_bedstartp0 = m_topg_initial.allocate_proc0_copy();
  m_upliftp0   = m_topg_initial.allocate_proc0_copy();

  bool use_elastic_model = m_config->get_boolean("bed_deformation.lc.elastic_model");

  m_bdLC = NULL;

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      m_bdLC = new BedDeformLC(*m_config, use_elastic_model,
                               m_grid->Mx(), m_grid->My(), m_grid->dx(), m_grid->dy(),
                               4,     // use Z = 4 for now; to reduce global drift?
                               *m_Hstartp0, *m_bedstartp0, *m_upliftp0, *m_Hp0, *m_bedp0);
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}

PBLingleClark::~PBLingleClark() {
  if (m_bdLC != NULL) {
    delete m_bdLC;
  }
}

void PBLingleClark::init_with_inputs_impl(const IceModelVec2S &bed,
                                          const IceModelVec2S &bed_uplift,
                                          const IceModelVec2S &ice_thickness) {
  m_t_beddef_last = m_grid->ctx()->time()->start();

  m_t  = GSL_NAN;
  m_dt = GSL_NAN;

  ice_thickness.put_on_proc0(*m_Hstartp0);
  bed.put_on_proc0(*m_bedstartp0);
  bed_uplift.put_on_proc0(*m_upliftp0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      m_bdLC->uplift_init();
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  // this should be the last thing we do here
  m_topg_initial.copy_from(m_topg);
  m_topg_last.copy_from(m_topg);
}

//! Initialize the Lingle-Clark bed deformation model using uplift.
void PBLingleClark::init_impl() {
  m_log->message(2,
             "* Initializing the Lingle-Clark bed deformation model...\n");

  BedDef::init_impl();

  // Correct bed elevation.
  if (m_config->get_boolean("bed_deformation.lc.regrid_bed_special")) {
    correct_topg();
  }

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  this->init_with_inputs_impl(m_topg, m_uplift, *ice_thickness);
}

MaxTimestep PBLingleClark::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("bed_def lc");
}

void PBLingleClark::correct_topg() {

  options::String regrid_file("-regrid_file",
                              "Specifies the name of a file to regrid from");

  InputOptions opts = process_input_options(m_grid->com);

  // Stop if it was requested, but we're not bootstrapping *and* regridding.
  if (not (regrid_file.is_set() and opts.type == INIT_BOOTSTRAP)) {
    return;
  }

  PIO nc(m_grid->com, "guess_mode", regrid_file, PISM_READONLY);
  bool topg_initial_exists = nc.inq_var("topg_initial");
  bool topg_exists = nc.inq_var("topg");
  nc.close();

  // Stop if the regridding file does not have both topg and topg_initial.
  if (not (topg_initial_exists and topg_exists)) {
    return;
  }

  // Stop if the user asked to regrid topg (in this case no correction is necessary).
  options::StringSet regrid_vars("-regrid_vars", "Specifies regridding variables", "");

  if (regrid_vars.is_set()) {
    if (set_contains(regrid_vars, "topg")) {
      m_log->message(2,
                 "  Bed elevation correction requested, but -regrid_vars contains topg...\n");
      return;
    }
  }

  m_log->message(2,
             "  Correcting topg from the bootstrapping file '%s' by adding the effect\n"
             "  of the bed deformation from '%s'...\n",
             opts.filename.c_str(), regrid_file->c_str());

  IceModelVec2S topg_tmp;       // will be de-allocated at 'return 0' below.
  const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");
  topg_tmp.create(m_grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
  topg_tmp.set_attrs("model_state", "bedrock surface elevation (at the end of the previous run)",
                     "m", "bedrock_altitude");

  // Get topg and topg_initial from the regridding file.
  m_topg_initial.regrid(regrid_file, CRITICAL);
  topg_tmp.regrid(regrid_file, CRITICAL);

  // After bootstrapping, topg contains the bed elevation field from
  // -i file.

  topg_tmp.add(-1.0, m_topg_initial);
  // Now topg_tmp contains the change in bed elevation computed during the run
  // that produced -regrid_file.

  // Apply this change to topg from -i file:
  m_topg.add(1.0, topg_tmp);

  // Store the corrected topg as the new "topg_initial".
  m_topg_initial.copy_from(m_topg);
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

  ice_thickness.put_on_proc0(*m_Hp0);
  m_topg.put_on_proc0(*m_bedp0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {  // only processor zero does the step
      m_bdLC->step(dt_beddef, // time step, in seconds
                   t_final - m_grid->ctx()->time()->start()); // time since the start of the run, in seconds
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  m_topg.get_from_proc0(*m_bedp0);

  //! Finally, we need to update bed uplift and topg_last.
  compute_uplift(dt_beddef);
  m_topg_last.copy_from(m_topg);

  //! Increment the topg state counter. SIAFD relies on this!
  m_topg.inc_state_counter();
}

} // end of namespace bed
} // end of namespace pism
