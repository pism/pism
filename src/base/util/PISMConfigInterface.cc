/* Copyright (C) 2015 PISM Authors
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

#include <cassert>
#include <mpi.h>

#include "base/util/io/PIO.hh"
#include "PISMConfigInterface.hh"
#include "PISMUnits.hh"
#include "pism_const.hh"
#include "pism_options.hh"
#include "error_handling.hh"

namespace pism {

struct Config::Impl {
  Impl(units::System::Ptr sys)
    : unit_system(sys) {
    // empty
  }
  //! Unit system. @fixme: this should be moved to the Context class.
  units::System::Ptr unit_system;

  std::string filename;

  //! @brief Set of parameters set by the user. Used to warn about parameters that were set but were
  //! not used.
  std::set<std::string> parameters_set_by_user;
  //! @brief Set of parameters used in a run. Used to warn about parameters that were set but were
  //! not used.
  std::set<std::string> parameters_used;
};

Config::Config(units::System::Ptr system)
  : m_impl(new Impl(system)) {
  // empty
}

Config::~Config() {
  delete m_impl;
}

units::System::Ptr Config::unit_system() const {
  return m_impl->unit_system;
}

void Config::read(MPI_Comm com, const std::string &file) {

  PIO nc(com, "netcdf3"); // OK to use netcdf3

  nc.open(file, PISM_READONLY);

  this->read(nc);

  nc.close();
}

void Config::read(const PIO &nc) {
  this->read_impl(nc);

  m_impl->filename = nc.inq_filename();
}

void Config::write(const PIO &nc) const {
  this->write_impl(nc);
}

void Config::write(MPI_Comm com, const std::string &file, bool append) const {

  PIO nc(com, "netcdf3"); // OK to use netcdf3

  IO_Mode mode = append ? PISM_READWRITE : PISM_READWRITE_MOVE;

  nc.open(file, mode);

  this->write(nc);

  nc.close();
}

//! \brief Returns the name of the file used to initialize the database.
std::string Config::filename() const {
  return m_impl->filename;
}

void Config::import_from(const Config &other) {
  Doubles doubles = other.all_doubles();
  Strings strings = other.all_strings();
  Booleans booleans = other.all_booleans();

  Doubles::const_iterator i;
  for (i = doubles.begin(); i != doubles.end(); ++i) {
    this->set_double(i->first, i->second, USER);
  }

  Strings::const_iterator j;
  for (j = strings.begin(); j != strings.end(); ++j) {
    this->set_string(j->first, j->second, USER);
  }

  Booleans::const_iterator k;
  for (k = booleans.begin(); k != booleans.end(); ++k) {
    this->set_boolean(k->first, k->second, USER);
  }
}

//! @brief Update values from the other config variable, overwriting present values but avoiding
//! adding new ones.
void Config::update_from(const Config &other) {
  Doubles doubles = other.all_doubles();
  Strings strings = other.all_strings();
  Booleans booleans = other.all_booleans();

  Doubles::const_iterator i;
  for (i = doubles.begin(); i != doubles.end(); ++i) {
    if (this->is_set(i->first)) {
      this->set_double(i->first, i->second, USER);
    }
  }

  Strings::const_iterator j;
  for (j = strings.begin(); j != strings.end(); ++j) {
    if (this->is_set(j->first)) {
      this->set_string(j->first, j->second, USER);
    }
  }

  Booleans::const_iterator k;
  for (k = booleans.begin(); k != booleans.end(); ++k) {
    if (this->is_set(k->first)) {
      this->set_boolean(k->first, k->second, USER);
    }
  }
}

const std::set<std::string>& Config::parameters_set_by_user() const {
  return m_impl->parameters_set_by_user;
}

const std::set<std::string>& Config::parameters_used() const {
  return m_impl->parameters_used;
}

bool Config::is_set(const std::string &name) const {
  return this->is_set_impl(name);
}

Config::Doubles Config::all_doubles() const {
  return this->all_doubles_impl();
}

double Config::get_double(const std::string &name) const {
  m_impl->parameters_used.insert(name);
  return this->get_double_impl(name);
}

double Config::get_double(const std::string &name,
                          const std::string &u1, const std::string &u2) const {
  double value = this->get_double(name);
  return units::convert(m_impl->unit_system, value, u1, u2);
}

void Config::set_double(const std::string &name, double value,
                        Config::SettingFlag flag) {
  if (flag == USER) {
    m_impl->parameters_set_by_user.insert(name);
  }

  this->set_double_impl(name, value);
}

Config::Strings Config::all_strings() const {
  return this->all_strings_impl();
}

std::string Config::get_string(const std::string &name) const {
  m_impl->parameters_used.insert(name);
  return this->get_string_impl(name);
}

void Config::set_string(const std::string &name,
                        const std::string &value,
                        Config::SettingFlag flag) {
  if (flag == USER) {
    m_impl->parameters_set_by_user.insert(name);
  }

  this->set_string_impl(name, value);
}

Config::Booleans Config::all_booleans() const {
  return this->all_booleans_impl();
}

bool Config::get_boolean(const std::string& name) const {
  m_impl->parameters_used.insert(name);
  return this->get_boolean_impl(name);
}

void Config::set_boolean(const std::string& name, bool value,
                         Config::SettingFlag flag) {
  if (flag == USER) {
    m_impl->parameters_set_by_user.insert(name);
  }

  this->set_boolean_impl(name, value);
}

void print_config(int verbosity_threshhold, MPI_Comm com, const Config &config) {
  const int v = verbosity_threshhold;

  verbPrintf(v, com,
             "### Strings:\n"
             "###\n");

  Config::Strings strings = config.all_strings();
  Config::Strings::const_iterator j;
  for (j = strings.begin(); j != strings.end(); ++j) {
    std::string name  = j->first;
    std::string value = j->second;

    if (value.empty() or ends_with(name, "_doc") or ends_with(name, "_units")) {
      continue;
    }

    verbPrintf(v, com, "  %s = \"%s\"\n", name.c_str(), value.c_str());
  }

  verbPrintf(v, com,
             "### Doubles:\n"
             "###\n");

  Config::Doubles doubles = config.all_doubles();
  Config::Doubles::const_iterator k;
  for (k = doubles.begin(); k != doubles.end(); ++k) {
    std::string name  = k->first;
    double value = k->second;
    std::string units = strings[name + "_units"]; // will be empty if not set

    if (fabs(value) >= 1.0e7 or fabs(value) <= 1.0e-4) {
      // use scientific notation if a number is big or small
      verbPrintf(v, com, "  %s = %12.3e (%s)\n", name.c_str(), value, units.c_str());
    } else {
      verbPrintf(v, com, "  %s = %12.5f (%s)\n", name.c_str(), value, units.c_str());
    }
  }

  verbPrintf(v, com,
             "### Booleans:\n"
             "###\n");

  Config::Booleans booleans = config.all_booleans();
  Config::Booleans::const_iterator p;
  for (p = booleans.begin(); p != booleans.end(); ++p) {
    std::string name  = p->first;
    std::string value = p->second ? "true" : "false";

    verbPrintf(v, com, "  %s = %s\n", name.c_str(), value.c_str());
  }

  verbPrintf(v, com,
             "### List of configuration parameters ends here.\n"
             "###\n");
}

void print_unused_parameters(int verbosity_threshhold, MPI_Comm com,
                             const Config &config) {
  std::set<std::string> parameters_set = config.parameters_set_by_user();
  std::set<std::string> parameters_used = config.parameters_used();

  if (options::Bool("-options_left", "report unused options")) {
    verbosity_threshhold = getVerbosityLevel();
  }

  std::set<std::string>::const_iterator k;
  for (k = parameters_set.begin(); k != parameters_set.end(); ++k) {

    if (ends_with(*k, "_doc")) {
      continue;
    }

    if (parameters_used.find(*k) == parameters_used.end()) {
      verbPrintf(verbosity_threshhold, com,
                 "PISM WARNING: flag or parameter \"%s\" was set but was not used!\n",
                 k->c_str());

    }
  }
}

// command-line options

//! Get a flag from a command-line option.
/*!
  If called as `boolean_from_option("foo", "foo")`, checks both `-foo` and `-no_foo`.

  - if `-foo` is set, calls `set_boolean("foo", true)`,

  - if `-no_foo` is set, calls `set_boolean("foo", false)`,

  - if *both* are set, prints an error message and stops,

  - if none, does nothing.

*/
void Config::boolean_from_option(const std::string &name, const std::string &flag) {

  bool foo    = options::Bool("-" + name, get_string_impl(flag + "_doc"));
  bool no_foo = options::Bool("-no_" + name, get_string_impl(flag + "_doc"));

  if (foo and no_foo) {
    throw RuntimeError::formatted("Inconsistent command-line options: both -%s and -no_%s are set.\n",
                                  name.c_str(), name.c_str());
  }

  if (foo) {
    set_boolean(flag, true, USER);
  }

  if (no_foo) {
    set_boolean(flag, false, USER);
  }
}

//! Sets a configuration parameter from a command-line option.
/*!
  If called as scalar_from_option("foo", "foo"), checks -foo and calls set("foo", value).

  Does nothing if -foo was not set.

  Note that no unit conversion is performed; parameters should be stored in
  input units and converted as needed. (This allows saving parameters without
  converting again.)
*/
void Config::scalar_from_option(const std::string &name, const std::string &parameter) {
  options::Real option("-" + name,
                       get_string_impl(parameter + "_doc"),
                       get_double_impl(parameter));
  if (option.is_set()) {
    this->set_double(parameter, option, USER);
  }
}

void Config::string_from_option(const std::string &name, const std::string &parameter) {

  options::String value("-" + name,
                        get_string_impl(parameter + "_doc"),
                        get_string_impl(parameter));
  if (value.is_set()) {
    this->set_string(parameter, value, USER);
  }
}

//! \brief Set a keyword parameter from a command-line option.
/*!
 * This sets the parameter "parameter" after checking the "-name" command-line
 * option. This option requires an argument, which has to match one of the
 * keyword given in a comma-separated list "choices_list".
 */
void Config::keyword_from_option(const std::string &name,
                                 const std::string &parameter,
                                 const std::string &choices) {

  options::Keyword keyword("-" + name,
                           this->get_string_impl(parameter + "_doc"),
                           choices,
                           this->get_string_impl(parameter));

  if (keyword.is_set()) {
    this->set_string(parameter, keyword, USER);
  }
}

void Config::set_from_options() {

  this->keyword_from_option("periodicity", "grid_periodicity", "none,x,y,xy");
  this->keyword_from_option("z_spacing", "grid_ice_vertical_spacing", "quadratic,equal");

  // Energy modeling
  this->boolean_from_option("use_Kirchhoff_law", "use_Kirchhoff_law");
  this->boolean_from_option("varc", "use_linear_in_temperature_heat_capacity");
  this->boolean_from_option("vark",
                            "use_temperature_dependent_thermal_conductivity");

  this->boolean_from_option("bmr_in_cont", "include_bmr_in_continuity");

  {
    options::Keyword energy("-energy",
                            "choose the energy model (one of 'none', 'cold', 'enthalpy')",
                            "none,cold,enthalpy", "enthalpy");

    if (energy.is_set()) {
      if (energy == "none") {
        this->set_boolean("do_energy", false, USER);
        // Allow selecting cold ice flow laws in isothermal mode.
        this->set_boolean("do_cold_ice_methods", true, USER);
      } else if (energy == "cold") {
        this->set_boolean("do_energy", true, USER);
        this->set_boolean("do_cold_ice_methods", true, USER);
      } else if (energy == "enthalpy") {
        this->set_boolean("do_energy", true, USER);
        this->set_boolean("do_cold_ice_methods", false, USER);
      } else {
        assert(false and "this can't happen: options::Keyword validates its input");
      }
    }
  }

  // at bootstrapping, choose whether the method uses smb as upper boundary for
  // vertical velocity
  this->keyword_from_option("boot_temperature_heuristic",
                            "bootstrapping_temperature_heuristic", "smb,quartic_guess");

  this->scalar_from_option("low_temp", "global_min_allowed_temp");
  this->scalar_from_option("max_low_temps", "max_low_temp_count");

  // Sub-models
  this->boolean_from_option("age", "do_age");
  this->boolean_from_option("mass", "do_mass_conserve");

  // hydrology
  this->keyword_from_option("hydrology", "hydrology_model",
                            "null,routing,distributed");
  this->boolean_from_option("hydrology_use_const_bmelt",
                            "hydrology_use_const_bmelt");
  this->scalar_from_option("hydrology_const_bmelt",
                           "hydrology_const_bmelt");
  this->scalar_from_option("hydrology_tillwat_max",
                           "hydrology_tillwat_max");
  this->scalar_from_option("hydrology_tillwat_decay_rate",
                           "hydrology_tillwat_decay_rate");
  this->scalar_from_option("hydrology_hydraulic_conductivity",
                           "hydrology_hydraulic_conductivity");
  this->scalar_from_option("hydrology_thickness_power_in_flux",
                           "hydrology_thickness_power_in_flux");
  this->scalar_from_option("hydrology_gradient_power_in_flux",
                           "hydrology_gradient_power_in_flux");
  // additional to hydrology::Routing, these apply to hydrology::Distributed:
  this->scalar_from_option("hydrology_roughness_scale",
                           "hydrology_roughness_scale");
  this->scalar_from_option("hydrology_cavitation_opening_coefficient",
                           "hydrology_cavitation_opening_coefficient");
  this->scalar_from_option("hydrology_creep_closure_coefficient",
                           "hydrology_creep_closure_coefficient");
  this->scalar_from_option("hydrology_regularizing_porosity",
                           "hydrology_regularizing_porosity");

  // Time-stepping
  this->keyword_from_option("calendar", "calendar",
                            "standard,gregorian,proleptic_gregorian,noleap,365_day,360_day,julian,none");

  this->scalar_from_option("adapt_ratio",
                           "adaptive_timestepping_ratio");

  this->scalar_from_option("timestep_hit_multiples",
                           "timestep_hit_multiples");

  this->boolean_from_option("count_steps", "count_time_steps");
  this->scalar_from_option("max_dt", "maximum_time_step_years");


  // SIA-related
  this->scalar_from_option("bed_smoother_range", "bed_smoother_range");

  this->keyword_from_option("gradient", "surface_gradient_method",
                            "eta,haseloff,mahaffy");

  // rheology-related
  this->scalar_from_option("sia_n", "sia_Glen_exponent");
  this->scalar_from_option("ssa_n", "ssa_Glen_exponent");

  this->keyword_from_option("sia_flow_law", "sia_flow_law",
                            "arr,arrwarm,gk,gpbld,hooke,isothermal_glen,pb");

  this->keyword_from_option("ssa_flow_law", "ssa_flow_law",
                            "arr,arrwarm,gpbld,hooke,isothermal_glen,pb");

  this->scalar_from_option("sia_e", "sia_enhancement_factor");
  this->scalar_from_option("ssa_e", "ssa_enhancement_factor");

  this->boolean_from_option("e_age_coupling", "e_age_coupling");

  // This parameter is used by the Goldsby-Kohlstedt flow law.
  this->scalar_from_option("ice_grain_size", "ice_grain_size");

  this->boolean_from_option("grain_size_age_coupling",
                            "compute_grain_size_using_age");

  // SSA
  // Decide on the algorithm for solving the SSA
  this->keyword_from_option("ssa_method", "ssa_method", "fd,fem");

  this->scalar_from_option("ssa_eps",  "epsilon_ssa");
  this->scalar_from_option("ssa_maxi", "max_iterations_ssafd");
  this->scalar_from_option("ssa_rtol", "ssafd_relative_convergence");

  this->scalar_from_option("ssafd_nuH_iter_failure_underrelaxation", "ssafd_nuH_iter_failure_underrelaxation");

  this->boolean_from_option("ssa_dirichlet_bc", "ssa_dirichlet_bc");
  this->boolean_from_option("cfbc", "calving_front_stress_boundary_condition");

  // Basal sliding fiddles
  this->boolean_from_option("brutal_sliding", "brutal_sliding");
  this->scalar_from_option("brutal_sliding_scale","brutal_sliding_scale");

  this->scalar_from_option("sliding_scale_factor_reduces_tauc",
                           "sliding_scale_factor_reduces_tauc");

  // SSA Inversion

  this->keyword_from_option("inv_method","inv_ssa_method",
                            "sd,nlcg,ign,tikhonov_lmvm,tikhonov_cg,tikhonov_blmvm,tikhonov_lcl,tikhonov_gn");

  this->keyword_from_option("inv_design_param",
                            "inv_design_param","ident,trunc,square,exp");

  this->scalar_from_option("inv_target_misfit","inv_target_misfit");

  this->scalar_from_option("tikhonov_penalty","tikhonov_penalty_weight");
  this->scalar_from_option("tikhonov_atol","tikhonov_atol");
  this->scalar_from_option("tikhonov_rtol","tikhonov_rtol");
  this->scalar_from_option("tikhonov_ptol","tikhonov_ptol");

  this->keyword_from_option("inv_state_func",
                            "inv_state_func",
                            "meansquare,log_ratio,log_relative");
  this->keyword_from_option("inv_design_func",
                            "inv_design_func","sobolevH1,tv");

  this->scalar_from_option("inv_design_cL2","inv_design_cL2");
  this->scalar_from_option("inv_design_cH1","inv_design_cH1");
  this->scalar_from_option("inv_ssa_tv_exponent","inv_ssa_tv_exponent");
  this->scalar_from_option("inv_log_ratio_scale","inv_log_ratio_scale");

  // Basal strength
  this->scalar_from_option("till_cohesion", "till_cohesion");
  this->scalar_from_option("till_reference_void_ratio",
                           "till_reference_void_ratio");
  this->scalar_from_option("till_compressibility_coefficient",
                           "till_compressibility_coefficient");
  this->scalar_from_option("till_effective_fraction_overburden",
                           "till_effective_fraction_overburden");
  this->scalar_from_option("till_log_factor_transportable_water",
                           "till_log_factor_transportable_water");

  // read the comma-separated list of four values
  options::RealList topg_to_phi("-topg_to_phi", "phi_min, phi_max, topg_min, topg_max");
  if (topg_to_phi.is_set()) {
    if (topg_to_phi->size() != 4) {
      throw RuntimeError::formatted("option -topg_to_phi requires a comma-separated list with 4 numbers; got %d",
                                    (int)topg_to_phi->size());
    }
    this->set_boolean("till_use_topg_to_phi", true);
    this->set_double("till_topg_to_phi_phi_min", topg_to_phi[0]);
    this->set_double("till_topg_to_phi_phi_max", topg_to_phi[1]);
    this->set_double("till_topg_to_phi_topg_min", topg_to_phi[2]);
    this->set_double("till_topg_to_phi_topg_max", topg_to_phi[3]);
  }

  this->boolean_from_option("tauc_slippery_grounding_lines",
                            "tauc_slippery_grounding_lines");
  this->boolean_from_option("tauc_add_transportable_water",
                            "tauc_add_transportable_water");

  this->keyword_from_option("yield_stress", "yield_stress_model",
                            "constant,mohr_coulomb");

  // all basal strength models use this in ice-free areas
  this->scalar_from_option("high_tauc", "high_tauc");

  // controls regularization of plastic basal sliding law
  this->scalar_from_option("plastic_reg", "plastic_regularization");

  // "friction angle" in degrees. We allow -plastic_phi without an
  // argument: MohrCoulombYieldStress interprets that as "set
  // constant till friction angle using the default read from a config
  // file or an override file".
  bool plastic_phi_set = options::Bool("-plastic_phi", "use constant till_phi");
  if (plastic_phi_set) {
    this->scalar_from_option("plastic_phi", "default_till_phi");
  }

  // use pseudo plastic instead of pure plastic; see iMbasal.cc
  this->boolean_from_option("pseudo_plastic", "do_pseudo_plastic_till");

  // power in denominator on pseudo_plastic_uthreshold; typical is q=0.25; q=0 is pure plastic
  this->scalar_from_option("pseudo_plastic_q", "pseudo_plastic_q");

  // threshold; at this velocity tau_c is basal shear stress
  this->scalar_from_option("pseudo_plastic_uthreshold",
                           "pseudo_plastic_uthreshold");

  this->boolean_from_option("subgl", "sub_groundingline");

  // Ice shelves
  this->boolean_from_option("part_grid", "part_grid");

  this->boolean_from_option("part_grid_reduce_frontal_thickness",
                            "part_grid_reduce_frontal_thickness");

  this->boolean_from_option("part_redist", "part_redist");

  this->scalar_from_option("nu_bedrock", "nu_bedrock");
  bool nu_bedrock = options::Bool("-nu_bedrock", "constant viscosity near margins");
  if (nu_bedrock) {
    this->set_boolean("nu_bedrock_set", true, USER);
  }

  // fracture density
  this->boolean_from_option("fractures", "do_fracture_density");
  this->boolean_from_option("write_fd_fields", "write_fd_fields");
  this->scalar_from_option("fracture_softening",
                           "fracture_density_softening_lower_limit");

  // Calving
  this->string_from_option("calving", "calving_methods");

  this->scalar_from_option("thickness_calving_threshold", "thickness_calving_threshold");

  // evaluates the adaptive timestep based on a CFL criterion with respect to the eigenCalving rate
  this->boolean_from_option("cfl_eigen_calving", "cfl_eigen_calving");
  this->scalar_from_option("eigen_calving_K", "eigen_calving_K");

  this->boolean_from_option("kill_icebergs", "kill_icebergs");

  // Output
  this->keyword_from_option("o_order", "output_variable_order",
                            "xyz,yxz,zyx");

  this->keyword_from_option("o_format", "output_format",
                            "netcdf3,quilt,netcdf4_parallel,pnetcdf,hdf5");

  this->scalar_from_option("summary_vol_scale_factor_log10",
                           "summary_vol_scale_factor_log10");
  this->scalar_from_option("summary_area_scale_factor_log10",
                           "summary_area_scale_factor_log10");

  // Metadata
  this->string_from_option("title", "run_title");
  this->string_from_option("institution", "institution");

  // Skipping
  this->boolean_from_option("skip", "do_skip");
  this->scalar_from_option("skip_max", "skip_max");

  // Shortcuts

  // option "-pik" turns on a suite of PISMPIK effects (but NOT a calving choice,
  // and in particular NOT  "-calving eigen_calving")
  bool pik = options::Bool("-pik", "enable suite of PISM-PIK mechanisms");
  if (pik) {
    this->set_boolean("calving_front_stress_boundary_condition", true, USER);
    this->set_boolean("part_grid", true, USER);
    this->set_boolean("part_redist", true, USER);
    this->set_boolean("kill_icebergs", true, USER);
    this->set_boolean("sub_groundingline", true, USER);
  }

  if (this->get_string("calving_methods").find("eigen_calving") != std::string::npos) {
    this->set_boolean("part_grid", true, USER);
    // eigen-calving requires a wider stencil:
    this->set_double("grid_max_stencil_width", 3);
  }

  // all calving mechanisms require iceberg removal
  if (this->get_string("calving_methods").empty() == false) {
    this->set_boolean("kill_icebergs", true, USER);
  }

  // kill_icebergs requires part_grid
  if (this->get_boolean("kill_icebergs")) {
    this->set_boolean("part_grid", true, USER);
  }

  this->keyword_from_option("stress_balance", "stress_balance_model",
                            "none,prescribed_sliding,sia,ssa,prescribed_sliding+sia,ssa+sia");

  bool test_climate_models = options::Bool("-test_climate_models",
                                           "Disable ice dynamics to test climate models");
  if (test_climate_models) {
    this->set_string("stress_balance_model", "none", USER);
    this->set_boolean("do_energy", false, USER);
    this->set_boolean("do_age", false, USER);
    // let the user decide if they want to use "-no_mass" or not
  }

  this->keyword_from_option("bed_def",
                            "bed_deformation_model", "none,iso,lc");
  this->boolean_from_option("bed_def_lc_elastic_model", "bed_def_lc_elastic_model");

  this->boolean_from_option("dry", "is_dry_simulation");

  this->boolean_from_option("clip_shelf_base_salinity",
                            "ocean_three_equation_model_clip_salinity");

  this->scalar_from_option("meltfactor_pik", "ocean_pik_melt_factor");

  // old options
  options::deprecated("-sliding_scale_brutal",
                      "-brutal_sliding' and '-brutal_sliding_scale");
  options::deprecated("-ssa_sliding", "-stress_balance ...");
  options::deprecated("-ssa_floating_only", "-stress_balance ...");
  options::deprecated("-sia", "-stress_balance ...");
  options::deprecated("-no_sia", "-stress_balance ...");
  options::deprecated("-hold_tauc", "-yield_stress constant");
  options::deprecated("-ocean_kill", "-calving ocean_kill -ocean_kill_file foo.nc");
  options::deprecated("-eigen_calving", "-calving eigen_calving -eigen_calving_K XXX");
  options::deprecated("-calving_at_thickness",
                      "-calving thickness_calving -thickness_calving_threshold XXX");
  options::deprecated("-float_kill", "-calving float_kill");
  options::deprecated("-no_energy", "-energy none");
  options::deprecated("-cold", "-energy cold");
}

} // end of namespace pism
