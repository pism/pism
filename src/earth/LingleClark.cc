// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2017, 2018, 2019, 2020, 2021 Constantine Khroulev
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

#include "LingleClark.hh"

#include "pism/util/io/File.hh"
#include "pism/util/Time.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Vars.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/fftw_utilities.hh"
#include "LingleClarkSerial.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace bed {

LingleClark::LingleClark(IceGrid::ConstPtr grid)
  : BedDef(grid),
    m_total_displacement(m_grid, "bed_displacement"),
    m_relief(m_grid, "bed_relief"),
    m_load_thickness(grid, "load_thickness"),
    m_elastic_displacement(grid, "elastic_bed_displacement") {

  m_time_name = m_config->get_string("time.dimension_name") + "_lingle_clark";
  m_t_last = time().current();
  m_update_interval = m_config->get_number("bed_deformation.lc.update_interval", "seconds");
  m_t_eps = m_config->get_number("time_stepping.resolution", "seconds");

  if (m_update_interval < 1.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid bed_deformation.lc.update_interval = %f seconds",
                                  m_update_interval);
  }

  // A work vector. This storage is used to put thickness change on rank 0 and to get the plate
  // displacement change back.
  m_total_displacement.set_attrs("internal",
                                 "total (viscous and elastic) displacement "
                                 "in the Lingle-Clark bed deformation model",
                                 "meters", "meters", "", 0);

  m_work0 = m_total_displacement.allocate_proc0_copy();

  m_relief.set_attrs("internal",
                     "bed relief relative to the modeled bed displacement",
                     "meters", "meters", "", 0);

  bool use_elastic_model = m_config->get_flag("bed_deformation.lc.elastic_model");

  m_elastic_displacement.set_attrs("model state",
                                   "elastic part of the displacement in the "
                                   "Lingle-Clark bed deformation model; "
                                   "see :cite:`BLKfastearth`", "meters", "meters", "", 0);
  m_elastic_displacement0 = m_elastic_displacement.allocate_proc0_copy();

  const int
    Mx = m_grid->Mx(),
    My = m_grid->My(),
    Z  = m_config->get_number("bed_deformation.lc.grid_size_factor"),
    Nx = Z*(Mx - 1) + 1,
    Ny = Z*(My - 1) + 1;

  const double
    Lx = Z * (m_grid->x0() - m_grid->x(0)),
    Ly = Z * (m_grid->y0() - m_grid->y(0));

  m_extended_grid = IceGrid::Shallow(m_grid->ctx(),
                                     Lx, Ly,
                                     m_grid->x0(), m_grid->y0(),
                                     Nx, Ny, CELL_CORNER, NOT_PERIODIC);

  m_viscous_displacement.reset(new array::Scalar(m_extended_grid,
                                                 "viscous_bed_displacement"));
  m_viscous_displacement->set_attrs("model state",
                                    "bed displacement in the viscous half-space "
                                    "bed deformation model; "
                                    "see BuelerLingleBrown",
                                    "meters", "meters", "", 0);

  // coordinate variables of the extended grid should have different names
  m_viscous_displacement->metadata().x().set_name("x_lc");
  m_viscous_displacement->metadata().y().set_name("y_lc");

  // do not point to auxiliary coordinates "lon" and "lat".
  m_viscous_displacement->metadata()["coordinates"] = "";

  m_viscous_displacement0 = m_viscous_displacement->allocate_proc0_copy();

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      m_serial_model.reset(new LingleClarkSerial(m_log, *m_config, use_elastic_model,
                                                 Mx, My,
                                                 m_grid->dx(), m_grid->dy(),
                                                 Nx, Ny));
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}

LingleClark::~LingleClark() {
  // empty
}

/*!
 * Initialize the model by computing the viscous bed displacement using uplift and the elastic
 * response using ice thickness.
 *
 * Then compute the bed relief as the difference between bed elevation and total bed displacement.
 *
 * This method has to initialize m_viscous_displacement, m_elastic_displacement,
 * m_total_displacement, and m_relief.
 */
void LingleClark::bootstrap_impl(const array::Scalar &bed_elevation,
                                 const array::Scalar &bed_uplift,
                                 const array::Scalar &ice_thickness,
                                 const array::Scalar &sea_level_elevation) {
  m_t_last = time().current();

  m_topg_last.copy_from(bed_elevation);

  compute_load(bed_elevation, ice_thickness, sea_level_elevation,
               m_load_thickness);

  std::shared_ptr<petsc::Vec> thickness0 = m_load_thickness.allocate_proc0_copy();

  // initialize the plate displacement
  {
    bed_uplift.put_on_proc0(*m_work0);
    m_load_thickness.put_on_proc0(*thickness0);

    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {
        PetscErrorCode ierr = 0;

        m_serial_model->bootstrap(*thickness0, *m_work0);

        ierr = VecCopy(m_serial_model->total_displacement(), *m_work0);
        PISM_CHK(ierr, "VecCopy");

        ierr = VecCopy(m_serial_model->viscous_displacement(), *m_viscous_displacement0);
        PISM_CHK(ierr, "VecCopy");

        ierr = VecCopy(m_serial_model->elastic_displacement(), *m_elastic_displacement0);
        PISM_CHK(ierr, "VecCopy");
      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();
  }

  m_viscous_displacement->get_from_proc0(*m_viscous_displacement0);

  m_elastic_displacement.get_from_proc0(*m_elastic_displacement0);

  m_total_displacement.get_from_proc0(*m_work0);

  // compute bed relief
  m_topg.add(-1.0, m_total_displacement, m_relief);
}

/*!
 * Return the load response matrix for the elastic response.
 *
 * This method is used for testing only.
 */
array::Scalar::Ptr LingleClark::elastic_load_response_matrix() const {
  array::Scalar::Ptr result(new array::Scalar(m_extended_grid, "lrm"));

  int
    Nx = m_extended_grid->Mx(),
    Ny = m_extended_grid->My();

  auto lrm0 = result->allocate_proc0_copy();

  {
    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {
        std::vector<std::complex<double> > array(Nx * Ny);

        m_serial_model->compute_load_response_matrix((fftw_complex*)array.data());

        get_real_part((fftw_complex*)array.data(), 1.0, Nx, Ny, Nx, Ny, 0, 0, *lrm0);
      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();
  }

  result->get_from_proc0(*lrm0);

  return result;
}

/*! Initialize the Lingle-Clark bed deformation model.
 *
 * Inputs:
 *
 * - bed topography,
 * - ice thickness,
 * - plate displacement (either read from a file or bootstrapped using uplift) and
 *   possibly re-gridded.
 */
void LingleClark::init_impl(const InputOptions &opts, const array::Scalar &ice_thickness,
                            const array::Scalar &sea_level_elevation) {
  m_log->message(2, "* Initializing the Lingle-Clark bed deformation model...\n");

  if (opts.type == INIT_RESTART or opts.type == INIT_BOOTSTRAP) {
    File input_file(m_grid->com, opts.filename, PISM_NETCDF3, PISM_READONLY);

    if (input_file.find_variable(m_time_name)) {
      input_file.read_variable(m_time_name, {0}, {1}, &m_t_last);
    } else {
      m_t_last = time().current();
    }
  } else {
    m_t_last = time().current();
  }

  // Initialize bed topography and uplift maps.
  BedDef::init_impl(opts, ice_thickness, sea_level_elevation);

  m_topg_last.copy_from(m_topg);

  if (opts.type == INIT_RESTART) {
    // Set viscous displacement by reading from the input file.
    m_viscous_displacement->read(opts.filename, opts.record);
    // Set elastic displacement by reading from the input file.
    m_elastic_displacement.read(opts.filename, opts.record);
  } else if (opts.type == INIT_BOOTSTRAP) {
    this->bootstrap(m_topg, m_uplift, ice_thickness, sea_level_elevation);
  } else {
    // do nothing
  }

  // Try re-gridding plate_displacement.
  regrid("Lingle-Clark bed deformation model",
         *m_viscous_displacement, REGRID_WITHOUT_REGRID_VARS);
  regrid("Lingle-Clark bed deformation model",
         m_elastic_displacement, REGRID_WITHOUT_REGRID_VARS);

  compute_load(m_topg, ice_thickness, sea_level_elevation,
               m_load_thickness);

  // Now that viscous displacement and elastic displacement are finally initialized,
  // put them on rank 0 and initialize the serial model itself.
  {
    m_viscous_displacement->put_on_proc0(*m_viscous_displacement0);
    m_elastic_displacement.put_on_proc0(*m_work0);

    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {  // only processor zero does the work
        PetscErrorCode ierr = 0;

        m_serial_model->init(*m_viscous_displacement0, *m_work0);

        ierr = VecCopy(m_serial_model->total_displacement(), *m_work0);
        PISM_CHK(ierr, "VecCopy");
      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();
  }

  m_total_displacement.get_from_proc0(*m_work0);

  // compute bed relief
  m_topg.add(-1.0, m_total_displacement, m_relief);
}

MaxTimestep LingleClark::max_timestep_impl(double t) const {

  if (t < m_t_last) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "time %f is less than the previous time %f",
                                  t, m_t_last);
  }

  // Find the smallest time of the form m_t_last + k * m_update_interval that is greater
  // than t
  double k = ceil((t - m_t_last) / m_update_interval);

  double
    t_next = m_t_last + k * m_update_interval,
    dt_max = t_next - t;

  if (dt_max < m_t_eps) {
    dt_max = m_update_interval;
  }

  return MaxTimestep(dt_max, "bed_def lc");
}

/*!
 * Get total bed displacement on the PISM grid.
 */
const array::Scalar& LingleClark::total_displacement() const {
  return m_total_displacement;
}

const array::Scalar& LingleClark::viscous_displacement() const {
  return *m_viscous_displacement;
}

const array::Scalar& LingleClark::elastic_displacement() const {
  return m_elastic_displacement;
}

const array::Scalar& LingleClark::relief() const {
  return m_relief;
}

void LingleClark::step(const array::Scalar &ice_thickness,
                       const array::Scalar &sea_level_elevation,
                       double dt) {

  compute_load(m_topg, ice_thickness, sea_level_elevation,
               m_load_thickness);

  m_load_thickness.put_on_proc0(*m_work0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {  // only processor zero does the step
      PetscErrorCode ierr = 0;

      m_serial_model->step(dt, *m_work0);

      ierr = VecCopy(m_serial_model->total_displacement(), *m_work0);
      PISM_CHK(ierr, "VecCopy");

      ierr = VecCopy(m_serial_model->viscous_displacement(), *m_viscous_displacement0);
      PISM_CHK(ierr, "VecCopy");

      ierr = VecCopy(m_serial_model->elastic_displacement(), *m_elastic_displacement0);
      PISM_CHK(ierr, "VecCopy");
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  m_viscous_displacement->get_from_proc0(*m_viscous_displacement0);

  m_elastic_displacement.get_from_proc0(*m_elastic_displacement0);

  m_total_displacement.get_from_proc0(*m_work0);

  // Update bed elevation using bed displacement and relief.
  {
    m_total_displacement.add(1.0, m_relief, m_topg);
    // Increment the topg state counter. SIAFD relies on this!
    m_topg.inc_state_counter();
  }

  //! Finally, we need to update bed uplift and topg_last.
  compute_uplift(m_topg, m_topg_last, dt, m_uplift);
  m_topg_last.copy_from(m_topg);
}

//! Update the Lingle-Clark bed deformation model.
void LingleClark::update_impl(const array::Scalar &ice_thickness,
                              const array::Scalar &sea_level_elevation,
                              double t, double dt) {

  double
    t_next  = m_t_last + m_update_interval,
    t_final = t + dt;

  if (t_final < m_t_last) {
    throw RuntimeError(PISM_ERROR_LOCATION, "cannot go back in time");
  }

  if (std::abs(t_next - t_final) < m_t_eps) { // reached the next update time
    double dt_beddef = t_final - m_t_last;
    step(ice_thickness, sea_level_elevation, dt_beddef);
    m_t_last = t_final;
  }
}

void LingleClark::define_model_state_impl(const File &output) const {
  BedDef::define_model_state_impl(output);
  m_viscous_displacement->define(output);
  m_elastic_displacement.define(output);

  if (not output.find_variable(m_time_name)) {
    output.define_variable(m_time_name, PISM_DOUBLE, {});

    output.write_attribute(m_time_name, "long_name",
                        "time of the last update of the Lingle-Clark bed deformation model");
    output.write_attribute(m_time_name, "calendar", time().calendar());
    output.write_attribute(m_time_name, "units", time().units_string());
  }
}

void LingleClark::write_model_state_impl(const File &output) const {
  BedDef::write_model_state_impl(output);

  m_viscous_displacement->write(output);
  m_elastic_displacement.write(output);

  output.write_variable(m_time_name, {0}, {1}, &m_t_last);
}

DiagnosticList LingleClark::diagnostics_impl() const {
  DiagnosticList result = {
    {"viscous_bed_displacement", Diagnostic::wrap(*m_viscous_displacement)},
    {"elastic_bed_displacement", Diagnostic::wrap(m_elastic_displacement)},
  };
  return combine(result, BedDef::diagnostics_impl());
}

} // end of namespace bed
} // end of namespace pism
