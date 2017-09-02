// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "PO_runoff_SMB.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace ocean {

#define T0_RELTOL   1.0e-14
#define ITER_MAXED_OUT 999

/* parameters needed for root problem: */
struct sgl_params {
  double A, B, a, b, alpha, beta;
};
  
/* the root problem is to make this function zero: */
double temp_to_sgl(double T0, void *params) {
  struct sgl_params *p = (struct sgl_params *) params;
  return (p->A * pow(p->a * T0 + p->b, p->alpha) + p->B) * pow(T0, p->beta) - 1;
};

  
Runoff_SMB::Runoff_SMB(IceGrid::ConstPtr g, OceanModel* in)
  : PScalarForcing<OceanModel,OceanModifier>(g, in) {

  m_option_prefix = "-ocean_runoff_smb";
  m_offset_name = "delta_T";
  
  m_offset = new Timeseries(*m_grid, m_offset_name, m_config->get_string("time.dimension_name"));
  m_offset->variable().set_string("units", "Kelvin");
  m_offset->variable().set_string("long_name", "air temperature offsets");
  m_offset->dimension().set_string("units", m_grid->ctx()->time()->units_string());

  m_temp_to_runoff_a = m_config->get_double("surface.temp_to_runoff_a");
  m_temp_to_runoff_b = m_config->get_double("surface.temp_to_runoff_b");

  m_runoff_to_ocean_melt_a = m_config->get_double("ocean.runoff_to_ocean_melt_a");
  m_runoff_to_ocean_melt_b = m_config->get_double("ocean.runoff_to_ocean_melt_b");

  m_runoff_to_ocean_melt_power_alpha = m_config->get_double("ocean.runoff_to_ocean_melt_power_alpha");
  m_runoff_to_ocean_melt_power_beta = m_config->get_double("ocean.runoff_to_ocean_melt_power_beta");

  int status = 0, iter = 0,  max_iter = 200;
  double T0 = 2.0, T0_lo = 0.0, T0_hi = 5.0;
  
  const gsl_root_fsolver_type *solvT;
  gsl_root_fsolver *solv;
  gsl_function F;
  struct sgl_params params;

  params.A =  m_runoff_to_ocean_melt_a;
  params.B =  m_runoff_to_ocean_melt_b;
  params.a =  m_temp_to_runoff_a;
  params.b =  m_temp_to_runoff_b;
  params.alpha =  m_runoff_to_ocean_melt_power_alpha;
  params.beta =  m_runoff_to_ocean_melt_power_beta;
  
  F.function = &temp_to_sgl;
  F.params   = &params;
  
  solvT      = gsl_root_fsolver_brent; /* faster than bisection but still bracketing */
  solv       = gsl_root_fsolver_alloc(solvT);
  gsl_root_fsolver_set(solv, &F, T0_lo, T0_hi);
  
  iter = 0;
  do {
    iter++;
    
    status = gsl_root_fsolver_iterate(solv);
    
    if (status != GSL_SUCCESS) {
      goto cleanup;
    }

    T0    = gsl_root_fsolver_root(solv);
    T0_lo = gsl_root_fsolver_x_lower(solv);
    T0_hi = gsl_root_fsolver_x_upper(solv);
    
    status = gsl_root_test_interval(T0_lo, T0_hi, 0, T0_RELTOL);
  } while ((status == GSL_CONTINUE) && (iter < max_iter));
  
  if (iter >= max_iter) {
    printf("!!!ERROR: root finding iteration reached maximum iterations; QUITTING!\n");
    goto cleanup;
  }
  
  m_log->message(2,
             "  - root solver found offset T_0 = %2.4f\n", T0);

 cleanup:
  gsl_root_fsolver_free(solv);
  m_current_forcing_0 = T0;
  
}

Runoff_SMB::~Runoff_SMB()
{
  // empty
}

void Runoff_SMB::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_input_model->init();

  m_log->message(2,
             "* Initializing ice shelf base mass flux forcing using scalar multiplier...\n"
             "*   derived from delta_T air temperature modifier\n");

  init_internal();

}

MaxTimestep Runoff_SMB::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("ocean runoff_SMB");
}
  
void Runoff_SMB::shelf_base_mass_flux_impl(IceModelVec2S &result) const {
  m_input_model->shelf_base_mass_flux(result);

  double a = m_temp_to_runoff_a,
    b = m_temp_to_runoff_b,
    A = m_runoff_to_ocean_melt_a,
    B = m_runoff_to_ocean_melt_b,
    alpha = m_runoff_to_ocean_melt_power_alpha,
    beta = m_runoff_to_ocean_melt_power_beta,
    m_scale = ((A * pow(a * m_current_forcing + b, alpha) + B) * pow(m_current_forcing, beta)) / (((A * pow(a * m_current_forcing_0 + b, alpha) + B)) * pow(m_current_forcing_0, beta));

  if (isnan(m_scale)) {
      m_scale = 0.0;
    }

  m_log->message(5,
                 "  - T_a = %2.1f, ocs = %2.4f\n", m_current_forcing, m_scale + 1);
  
  result.scale(1 + m_scale);
}

} // end of namespace ocean
} // end of namespace pism
