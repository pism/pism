/* Copyright (C) 2016, 2017 PISM Authors
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

#include "AgeModel.hh"

#include "pism/age/AgeColumnSystem.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Vars.hh"
#include "pism/util/io/PIO.hh"

namespace pism {

AgeModelInputs::AgeModelInputs(const IceModelVec2S *thickness,
                               const IceModelVec3 *u,
                               const IceModelVec3 *v,
                               const IceModelVec3 *w)
  : ice_thickness(thickness), u3(u), v3(v), w3(w) {
  // empty
}

AgeModelInputs::AgeModelInputs() {
  ice_thickness = NULL;
  u3            = NULL;
  v3            = NULL;
  w3            = NULL;
}

static void check_input(const IceModelVec *ptr, const char *name) {
  if (ptr == NULL) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "ice age model input %s was not provided", name);
  }
}

void AgeModelInputs::check() const {
  check_input(ice_thickness, "ice_thickness");
  check_input(u3, "u3");
  check_input(v3, "v3");
  check_input(w3, "w3");
}

AgeModel::AgeModel(IceGrid::ConstPtr grid, stressbalance::StressBalance *stress_balance)
  : Component(grid), m_stress_balance(stress_balance) {

  // FIXME: should be able to use width=1...
  const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");

  m_ice_age.create(m_grid, "age", WITH_GHOSTS, WIDE_STENCIL);
  m_ice_age.set_attrs("model_state", "age of ice", "s", "" /* no standard name*/);
  m_ice_age.metadata().set_string("glaciological_units", "years");
  m_ice_age.metadata().set_double("valid_min", 0.0);

  m_work.create(m_grid,"work_vector",WITHOUT_GHOSTS);
  m_work.set_attrs("internal", "new values of age during time step", "s", "");
}

/*!
Let \f$\tau(t,x,y,z)\f$ be the age of the ice.  Denote the three-dimensional
velocity field within the ice fluid as \f$(u,v,w)\f$.  The age equation
is \f$d\tau/dt = 1\f$, that is, ice may move but it gets one year older in one
year.  Thus
    \f[ \frac{\partial \tau}{\partial t} + u \frac{\partial \tau}{\partial x}
        + v \frac{\partial \tau}{\partial y} + w \frac{\partial \tau}{\partial z} = 1 \f]
This equation is purely advective and hyperbolic.  The right-hand side is "1" as
long as age \f$\tau\f$ and time \f$t\f$ are measured in the same units.
Because the velocity field is incompressible, \f$\nabla \cdot (u,v,w) = 0\f$,
we can rewrite the equation as
    \f[ \frac{\partial \tau}{\partial t} + \nabla \left( (u,v,w) \tau \right) = 1 \f]
There is a conservative first-order numerical method; see AgeColumnSystem::solveThisColumn().

The boundary condition is that when the ice falls as snow it has age zero.
That is, \f$\tau(t,x,y,h(t,x,y)) = 0\f$ in accumulation areas.  There is no
boundary condition elsewhere on the ice upper surface, as the characteristics
go outward in the ablation zone.  If the velocity in the bottom cell of ice
is upward (\f$w>0\f$) then we also apply a zero age boundary condition,
\f$\tau(t,x,y,0) = 0\f$.  This is the case where ice freezes on at the base,
either grounded basal ice freezing on stored water in till, or marine basal ice.
(Note that the water that is frozen-on as ice might be quite "old" in the sense
that its most recent time in the atmosphere was long ago; this comment is
relevant to any analysis which relates isotope ratios to modeled age.)

The numerical method is a conservative form of first-order upwinding, but the
vertical advection term is computed implicitly.  Thus there is no CFL-type
stability condition from the vertical velocity; CFL is only for the horizontal
velocity.  We use a finely-spaced, equally-spaced vertical grid in the
calculation.  Note that the columnSystemCtx methods coarse_to_fine() and
fine_to_coarse() interpolate back and forth between this fine grid and
the storage grid.  The storage grid may or may not be equally-spaced.  See
AgeColumnSystem::solve() for the actual method.
 */
void AgeModel::update(double t, double dt, const AgeModelInputs &inputs) {

  // fix a compiler warning
  (void) t;

  inputs.check();

  const IceModelVec2S &ice_thickness = *inputs.ice_thickness;

  const IceModelVec3
    &u3 = *inputs.u3,
    &v3 = *inputs.v3,
    &w3 = *inputs.w3;

  AgeColumnSystem system(m_grid->z(), "age",
                         m_grid->dx(), m_grid->dy(), dt,
                         m_ice_age, u3, v3, w3); // linear system to solve in each column

  size_t Mz_fine = system.z().size();
  std::vector<double> x(Mz_fine);   // space for solution

  IceModelVec::AccessList list{&ice_thickness, &u3, &v3, &w3, &m_ice_age, &m_work};

  unsigned int Mz = m_grid->Mz();

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      system.init(i, j, ice_thickness(i, j));

      if (system.ks() == 0) {
        // if no ice, set the entire column to zero age
        m_work.set_column(i, j, 0.0);
      } else {
        // general case: solve advection PDE

        // solve the system for this column; call checks that params set
        system.solve(x);

        // put solution in IceModelVec3
        system.fine_to_coarse(x, i, j, m_work);

        // Ensure that the age of the ice is non-negative.
        //
        // FIXME: this is a kludge. We need to ensure that our numerical method has the maximum
        // principle instead. (We may still need this for correctness, though.)
        double *column = m_work.get_column(i, j);
        for (unsigned int k = 0; k < Mz; ++k) {
          if (column[k] < 0.0) {
            column[k] = 0.0;
          }
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_work.update_ghosts(m_ice_age);
}

const IceModelVec3 & AgeModel::age() const {
  return m_ice_age;
}

MaxTimestep AgeModel::max_timestep_impl(double t) const {
  // fix a compiler warning
  (void) t;

  if (m_stress_balance == NULL) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "AgeModel: no stress balance provided."
                                  " Cannot compute max. time step.");
  }

  return MaxTimestep(m_stress_balance->max_timestep_cfl_3d().dt_max.value(), "age model");
}

void AgeModel::init(const InputOptions &opts) {

  m_log->message(2, "* Initializing the age model...\n");


  double initial_age_years = m_config->get_double("age.initial_value", "years");

  if (opts.type == INIT_RESTART) {
    PIO input_file(m_grid->com, "guess_mode", opts.filename, PISM_READONLY);

    if (input_file.inq_var("age")) {
      m_ice_age.read(input_file, opts.record);
    } else {
      m_log->message(2,
                     "PISM WARNING: input file '%s' does not have the 'age' variable.\n"
                     "  Setting it to %f years...\n",
                     opts.filename.c_str(), initial_age_years);
      m_ice_age.set(m_config->get_double("age.initial_value", "seconds"));
    }
  } else {
    m_log->message(2, " - setting initial age to %.4f years\n", initial_age_years);
    m_ice_age.set(m_config->get_double("age.initial_value", "seconds"));
  }

  regrid("Age Model", m_ice_age, REGRID_WITHOUT_REGRID_VARS);
}

void AgeModel::define_model_state_impl(const PIO &output) const {
  m_ice_age.define(output);
}

void AgeModel::write_model_state_impl(const PIO &output) const {
  m_ice_age.write(output);
}

} // end of namespace pism
