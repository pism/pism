// Copyright (C) 2004-2025 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef PISM_ICEMODEL_H
#define PISM_ICEMODEL_H

//! \file IceModel.hh Definition of class IceModel.
/*! \file IceModel.hh
  IceModel is a big class which is an ice flow model.  It contains all parts that
  are not well-defined, separated components.  Such components are better places
  to put sub-models that have a clear, general interface to the rest of an ice
  sheet model.

  IceModel has pointers to well-defined components, when they exist.

  IceModel generally interprets user options, and initializes components based on
  such options.  It manages the initialization sequences (%e.g. a restart from a
  file containing a complete model state, versus bootstrapping).
*/

#include <map>
#include <set>
#include <string>
#include <vector>
#include <memory>

// IceModel owns a bunch of fields, so we have to include this.
#include "pism/util/array/Vector.hh"
#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Time.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/geometry/GeometryEvolution.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/basalstrength/YieldStress.hh"
#include "pism/util/ScalarForcing.hh" // for use with std::unique_ptr

namespace pism {

namespace ocean {
class OceanModel;
class PyOceanModel;
namespace sea_level {
class SeaLevel;
}
}

namespace surface {
class SurfaceModel;
}

namespace hydrology {
class Hydrology;
}

namespace calving {
class EigenCalving;
class vonMisesCalving;
class FloatKill;
class HayhurstCalving;
class CalvingAtThickness;
class IcebergRemover;
}

class FractureDensity;

namespace energy {
class BedThermalUnit;
class Inputs;
class EnergyModelStats;
class EnergyModel;
}

namespace frontalmelt {
  class FrontalMelt;
}

namespace bed {
class BedDef;
}

namespace array {
class Forcing;
class CellType;
}

class Grid;
class AgeModel;
class Isochrones;
class Component;
class FrontRetreat;
class PrescribedRetreat;
class ScalarForcing;

class OutputWriter;

enum IceModelTerminationReason {PISM_DONE, PISM_CHEKPOINT, PISM_SIGNAL};

enum DiagnosticReport {DIAG_NONE, DIAG_ALL, DIAG_SPATIAL, DIAG_SCALAR, DIAG_JSON};

//! The base class for PISM. Contains all essential variables, parameters, and flags for modelling
//! an ice sheet.
class IceModel {
public:
  IceModel(std::shared_ptr<Grid> grid, const std::shared_ptr<Context> &context);

  // the destructor must be virtual merely because some members are virtual
  virtual ~IceModel();

  std::shared_ptr<Grid> grid() const;
  std::shared_ptr<Context> ctx() const;

  void init(DiagnosticReport report_type);

  /** Run PISM in the "standalone" mode. */
  IceModelTerminationReason run();

  /** Advance the current PISM run to a specific time */
  IceModelTerminationReason run_to(double run_end);

  void list_diagnostics(DiagnosticReport report_type) const;

  void write_final_output();

  const array::Scalar &calving() const;
  const array::Scalar &frontal_melt() const;
  const array::Scalar &forced_retreat() const;

  const stressbalance::StressBalance* stress_balance() const;
  const ocean::OceanModel* ocean_model() const;
  const energy::BedThermalUnit* bedrock_thermal_model() const;
  const energy::EnergyModel* energy_balance_model() const;
  const YieldStress* basal_yield_stress_model() const;
  const bed::BedDef* bed_deformation_model() const;

  /*!
   * Replace the ocean model with an implementation in Python.
   */
  void set_python_ocean_model(std::shared_ptr<ocean::PyOceanModel> model);

  const Geometry& geometry() const;
  const GeometryEvolution& geometry_evolution() const;

  double dt() const;

protected:
  virtual void allocate_submodels();
  virtual void allocate_stressbalance();
  virtual void allocate_age_model();
  virtual void allocate_isochrones();
  virtual void allocate_bed_deformation();
  virtual void allocate_bedrock_thermal_unit();
  virtual void allocate_energy_model();
  virtual void allocate_subglacial_hydrology();
  virtual void allocate_basal_yield_stress();
  virtual void allocate_couplers();
  virtual void allocate_geometry_evolution();
  virtual void allocate_iceberg_remover();

  virtual stressbalance::Inputs stress_balance_inputs();

  virtual energy::Inputs energy_model_inputs();

  virtual YieldStressInputs yield_stress_inputs();

  virtual void model_state_setup(InputOptions input_options);
  virtual void misc_setup(InputOptions input_options, DiagnosticReport report_type);
  virtual void init_calving();
  virtual void init_frontal_melt();
  virtual void init_front_retreat();
  virtual void update_diagnostics(double t, double dt);

  virtual std::map<std::string, Diagnostic::Ptr> allocate_spatial_diagnostics();
  virtual std::map<std::string, TSDiagnostic::Ptr> allocate_scalar_diagnostics();
  void deallocate_unused_diagnostics();

  void init_outputs(InputOptions options, DiagnosticReport report_type);

  /*!
   * Return the set of names of diagnostic quantities corresponding to a `keyword` (none,
   * small, medium, big_2d, big).
   */
  virtual std::set<std::string> output_variables(const std::string &keyword);

  /*!
   * Return the set of all the diagnostic variables corresponding to the set of requested
   * diagnostic *quantities* in `variable_names` (one quantity may map to two or more
   * variables, e.g. velocities (2 or 3 vector components) of quantities on the staggered
   * grid (2 grid offsets)).
   */
  virtual std::set<VariableMetadata>
  diagnostic_variables(const std::set<std::string> &variable_names) const;

  /*!
   * Return the set of "state" variables, i.e. variables that describe the state of
   * IceModel and of all its sub-models.
   */
  virtual std::set<VariableMetadata> state_variables() const;

  /*!
   * Return the set of "state" variables needed by the provided set of 2D and 3D
   * diagnostic quantities.
   */
  std::set<VariableMetadata>
  diagnostic_state_variables(const std::set<std::string> &variable_names) const;

  /*!
   * Return the set of "common" variables (i.e. ones written to most files): step counter,
   * model-years-per-processor-hour, wall clock time since start, configuration
   * parameters, grid mapping (if known).
   */
  std::set<VariableMetadata> common_metadata() const;

  //! Name of the output file
  std::string m_output_filename;

  // Set of diagnostic quantities to put in the output file:
  std::set<std::string> m_output_vars;

  //! Set of variables that will be written to the output file
  std::set<VariableMetadata> m_output_file_contents;

  void init_final_output();

  virtual double step(bool do_mass_continuity, bool do_skip);
  virtual void pre_step_hook();
  virtual void post_step_hook();

  void reset_counters();

  virtual void bootstrap_2d(const File &input_file);

  virtual void compute_lat_lon();

  virtual void restart_2d(const File &input_file, unsigned int record);
  virtual void initialize_2d();

  void define_time(const OutputFile &file, bool with_bounds = false) const;
  void define_variables(const OutputFile &file, const std::set<VariableMetadata> &variables) const;

  virtual void write_state(const OutputFile &file) const;

  virtual void write_run_stats(const OutputFile &file) const;

  virtual void write_diagnostics(const OutputFile &file,
                                 const std::set<std::string> &variable_names) const;

  std::string save_state_on_error(const std::string &suffix,
                                  const std::set<std::string> &additional_variables);

  //! Computational grid
  const std::shared_ptr<Grid> m_grid;
  //! Configuration flags and parameters
  std::shared_ptr<Config> m_config;
  //! Execution context
  std::shared_ptr<Context> m_ctx;
  //! Unit system
  const units::System::Ptr m_sys;
  //! Logger
  std::shared_ptr<Logger> m_log;
  //! Time manager
  std::shared_ptr<Time> m_time;

  std::shared_ptr<OutputWriter> m_output_writer;

  //! stores global attributes saved in a PISM output file
  VariableMetadata m_output_global_attributes;
  std::string m_output_history;

  //! the list of sub-models, for writing model states and obtaining diagnostics
  std::map<std::string,const Component*> m_submodels;

  std::unique_ptr<hydrology::Hydrology> m_subglacial_hydrology;
  std::shared_ptr<YieldStress> m_basal_yield_stress_model;

  std::shared_ptr<array::Forcing> m_surface_input_for_hydrology;

  std::shared_ptr<energy::BedThermalUnit> m_btu;
  std::shared_ptr<energy::EnergyModel> m_energy_model;

  std::shared_ptr<AgeModel> m_age_model;
  std::shared_ptr<Isochrones> m_isochrones;

  std::shared_ptr<calving::IcebergRemover>     m_iceberg_remover;
  std::shared_ptr<calving::FloatKill>          m_float_kill_calving;
  std::shared_ptr<calving::CalvingAtThickness> m_thickness_threshold_calving;
  std::shared_ptr<calving::EigenCalving>       m_eigen_calving;
  std::shared_ptr<calving::HayhurstCalving>    m_hayhurst_calving;
  std::shared_ptr<calving::vonMisesCalving>    m_vonmises_calving;
  std::shared_ptr<PrescribedRetreat>           m_prescribed_retreat;

  // scalar time-dependent scaling for retreat rates coming from eigen calving, von Mises
  // calving, or Hayhurst calving
  std::unique_ptr<ScalarForcing> m_calving_rate_factor;

  std::shared_ptr<FrontRetreat> m_front_retreat;

  std::shared_ptr<surface::SurfaceModel>      m_surface;
  std::shared_ptr<ocean::OceanModel>          m_ocean;
  std::shared_ptr<frontalmelt::FrontalMelt>   m_frontal_melt;
  std::shared_ptr<ocean::sea_level::SeaLevel> m_sea_level;

  std::shared_ptr<bed::BedDef> m_beddef;

  // state variables and some diagnostics/internals

  Geometry m_geometry;
  std::unique_ptr<GeometryEvolution> m_geometry_evolution;
  bool m_new_bed_elevation;

  //! ghosted
  array::Scalar2 m_basal_yield_stress;
  //! rate of production of basal meltwater (ice-equivalent); no ghosts
  array::Scalar m_basal_melt_rate;
  //! temperature at the top surface of the bedrock thermal layer
  array::Scalar m_bedtoptemp;

  std::shared_ptr<FractureDensity> m_fracture;

  //! mask to determine Dirichlet boundary locations for the sliding velocity
  array::Scalar2 m_velocity_bc_mask;
  //! Dirichlet boundary velocities
  array::Vector2 m_velocity_bc_values;

  //! Mask prescribing locations where ice thickness is held constant
  array::Scalar1 m_ice_thickness_bc_mask;

  // parameters
  //! mass continuity time step, s
  double m_dt;
  //! time of last update for enthalpy/temperature
  double m_t_TempAge;
  //! enthalpy/temperature and age time-steps
  double m_dt_TempAge;

  unsigned int m_skip_countdown;

  std::string m_adaptive_timestep_reason;

  std::string m_stdout_flags;

  unsigned int m_step_counter;

  // see iceModel.cc
  virtual void allocate_storage();

  struct TimesteppingInfo {
    double dt;
    std::string reason;
    unsigned int skip_counter;
  };
  virtual TimesteppingInfo max_timestep(unsigned int counter);

  virtual MaxTimestep max_timestep_diffusivity();
  virtual unsigned int skip_counter(double input_dt, double input_dt_diffusivity);

  // see energy.cc
  virtual void bedrock_thermal_model_step(double t, double dt);
  virtual void energy_step(double t, double dt);

  virtual void hydrology_step(double t, double dt);

  virtual void combine_basal_melt_rate(const Geometry &geometry,
                                       const array::Scalar &shelf_base_mass_flux,
                                       const array::Scalar &grounded_basal_melt_rate,
                                       array::Scalar &result);

  enum ConsistencyFlag {REMOVE_ICEBERGS, DONT_REMOVE_ICEBERGS};
  void enforce_consistency_of_geometry(ConsistencyFlag flag);

  /*!
   * Compute the mask (`result`) identifying "open ocean" cells, i.e. "ice free ocean"
   * cells that are reachable from open ocean cells at an edge of the domain.
   *
   * Here `result` is ghosted so that we can pass it to the connected component labeling
   * algorithm.
   */
  void identify_open_ocean(const array::CellType &cell_type, array::Scalar1 &result);

  virtual void front_retreat_step(double t, double dt);

  void compute_geometry_change(const array::Scalar &thickness,
                               const array::Scalar &Href,
                               const array::Scalar &thickness_old,
                               const array::Scalar &Href_old,
                               bool add_values,
                               array::Scalar &output);

  // see iMIO.cc
  virtual void regrid();

  // see iMfractures.cc
  virtual void update_fracture_density(double dt);

  // see iMreport.cc
  virtual double compute_temperate_base_fraction(double ice_area);
  virtual double compute_original_ice_fraction(double ice_volume);
  virtual void print_summary(bool tempAndAge, double dt);
  virtual void print_summary_line(bool printPrototype, bool tempAndAge,
                                  double delta_t,
                                  double volume, double area,
                                  double meltfrac, double max_diffusivity);


  // see iMutil.cc
  virtual int process_signals();
  virtual void append_history(const std::string &string);

  // working space (a convenience)
  static const int m_n_work2d = 4;
  mutable std::vector<std::shared_ptr<array::Scalar2>> m_work2d;

  std::shared_ptr<stressbalance::StressBalance> m_stress_balance;

  struct ThicknessChanges {
    ThicknessChanges(const std::shared_ptr<const Grid> &grid);

    // calving during the last time step
    array::Scalar calving;

    // frontal melt during the last time step
    array::Scalar frontal_melt;

    // parameterized retreat
    array::Scalar forced_retreat;
  };

  ThicknessChanges m_thickness_change;

  /*!
   * The set of variables that the "state" of IceModel consists of.
   */
  std::set<array::Array*> m_model_state;
  //! Requested spatially-variable diagnostics.
  std::map<std::string,Diagnostic::Ptr> m_diagnostics;
  //! Requested scalar diagnostics.
  std::map<std::string,TSDiagnostic::Ptr> m_ts_diagnostics;

  // This is related to the snapshot saving feature
  std::string m_snapshots_filename;
  std::shared_ptr<OutputFile> m_snapshot_file;
  bool m_split_snapshots;
  std::vector<double> m_snapshot_times;
  unsigned int m_current_snapshot;
  //! set of diagnostics to write to snapshot files (one diagnostic may correspond to
  //! multiple NetCDF variables)
  std::set<std::string> m_snapshot_vars;
  //! set of variables that will be written to snapshot files
  std::set<VariableMetadata> m_snapshot_file_contents;
  void init_snapshots();
  void write_snapshot();
  MaxTimestep snapshots_max_timestep(double my_t);

  //! file to write scalar time-series to
  std::shared_ptr<OutputFile> m_ts_file;
  //! requested times for scalar time-series
  std::shared_ptr<std::vector<double>> m_ts_times;
  std::set<std::string> m_ts_vars;
  void init_timeseries();
  void flush_timeseries();
  MaxTimestep ts_max_timestep(double my_t);

  // spatially-varying time-series
  std::string m_extra_filename;
  std::vector<double> m_extra_times;
  unsigned int m_next_extra;
  double m_last_extra;
  std::set<std::string> m_extra_vars;

  //! set of variables that will be written to extra files
  std::set<VariableMetadata> m_extra_file_contents;

  std::shared_ptr<OutputFile> m_extra_file;
  void init_extras();
  void write_extras();
  MaxTimestep extras_max_timestep(double my_t);

  // automatic checkpoints
  std::string m_checkpoint_filename;
  double m_last_checkpoint_time;
  std::set<std::string> m_checkpoint_vars;
  //! set of variables that will be written to checkpoint files
  std::set<VariableMetadata> m_checkpoint_file_contents;
  void init_checkpoints();
  bool write_checkpoint();

  // last time at which PISM hit a multiple of X years, see the configuration parameter
  // time_stepping.hit_multiples
  double m_timestep_hit_multiples_last_time;

  // diagnostic viewers; see iMviewers.cc
  virtual void update_viewers();
  virtual void view_field(const array::Array *field);
  std::map<std::string,
           std::vector<std::shared_ptr<petsc::Viewer> > > m_viewers;

private:
  double m_start_time;    // this is used in the wall-clock-time checkpoint code
};

MaxTimestep reporting_max_timestep(const std::vector<double> &times,
                                   double t,
                                   double eps,
                                   const std::string &description);

void bedrock_surface_temperature(const array::Scalar &sea_level,
                                 const array::CellType &cell_type,
                                 const array::Scalar &bed_topography,
                                 const array::Scalar &ice_thickness,
                                 const array::Scalar &basal_enthalpy,
                                 const array::Scalar &ice_surface_temperature,
                                 array::Scalar &result);

} // end of namespace pism

#endif /* PISM_ICEMODEL_H */
