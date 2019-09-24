#include "Example.hh"

#include "base/util/PISMConfigInterface.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_options.hh"
#include "base/util/MaxTimestep.hh"

namespace pism {
namespace ocean {

Example::Example(IceGrid::ConstPtr g)
  : OceanModel(g) {

  // assume that climate_forcing.buffer_size is big enough
  m_shelf_melt_rate.create(m_grid, "shelf_base_melt_rate",
                           m_config->get_double("climate_forcing.buffer_size"));
  m_shelf_melt_rate.set_attrs("internal", "shelf base melt rate", "m / second", "");
}

Example::~Example() {
  // empty
}

void Example::update_impl(double t, double dt) {
  m_t  = t;
  m_dt = dt;

  m_shelf_melt_rate.update(t, dt);

  // Use mid-point of the interval. (We restricted the time step, so
  // the choice of the point within the time step does not matter.)
  m_shelf_melt_rate.interp(t + 0.5 * dt);

  // Alternatively one could call. This does not require a time step restriction.
  // m_shelf_melt_rate.average(t, dt);
}

void Example::init_impl() {
  m_log->message(2, "* Initializing the example ocean model...\n");

  options::String input_file("-ocean_example_file", "Shelf melt rate input file.");

  if (input_file.is_set()) {
    m_log->message(2, "  Reading shelf base melt rate from %s...\n",
                   input_file->c_str());

    m_shelf_melt_rate.init(input_file, 0.0, 0.0);
  } else {
    m_shelf_melt_rate.init_constant(0.0);
  }
}

MaxTimestep Example::max_timestep_impl(double t) const {
  // Assume that temporal variations in the melt rate have to be resolved.
  return m_shelf_melt_rate.max_timestep(t);

  // Use this to disable the time step restriction
  return MaxTimestep("example ocean model");
}

void Example::shelf_base_temperature_impl(IceModelVec2S &result) const {
  // PISM uses MKS. This is obviously wrong, but this just an example.
  result.set(273.15);
}

void Example::sea_level_elevation_impl(double &result) const {
  // Also wrong.
  result = 0.0;
}

//! @brief Computes mass flux in [kg m-2 s-1], from assumption that
//! basal heat flux rate converts to mass flux.
void Example::shelf_base_mass_flux_impl(IceModelVec2S &result) const {
  result.copy_from(m_shelf_melt_rate);
}

} // end of namespape ocean
} // end of namespace pism
