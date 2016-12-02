// Copyright (C) 2004-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __iceModel_hh
#define __iceModel_hh

//! \file iceModel.hh Definition of class IceModel.
/*! \file iceModel.hh
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

// IceModel owns a bunch of fields, so we have to include this.
#include "base/util/iceModelVec.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/Context.hh"
#include "base/util/Logger.hh"
#include "base/util/PISMTime.hh"
#include "base/util/IceModelVec2CellType.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/MaxTimestep.hh"

namespace pism {

namespace ocean {
class OceanModel;
}

namespace surface {
class SurfaceModel;
}

namespace stressbalance {
class StressBalance;
}

namespace hydrology {
class Hydrology;
}

namespace calving {
class EigenCalving;
class vonMisesCalving;
class OceanKill;
class FloatKill;
class CalvingAtThickness;
class IcebergRemover;
}

class FrontalMelt;

namespace energy {
class BedThermalUnit;
class EnergyModelInputs;
class EnergyModelStats;
class EnergyModel;
}

namespace bed {
class BedDef;
}

// forward declarations
class IceGrid;
class YieldStress;
class AgeModel;
class IceModelVec2CellType;
class Component;

struct FractureFields {
  FractureFields(IceGrid::ConstPtr grid);

  IceModelVec2S density;
  IceModelVec2S growth_rate;
  IceModelVec2S healing_rate;
  IceModelVec2S flow_enhancement;
  IceModelVec2S age;
  IceModelVec2S toughness;

  //! major and minor principal components of horizontal strain-rate tensor (temporary storage)
  IceModelVec2 strain_rates;

  //! components of horizontal stress tensor along axes and shear stress (temporary storage)
  IceModelVec2 deviatoric_stresses;
};


//! The base class for PISM.  Contains all essential variables, parameters, and flags for modelling an ice sheet.
class IceModel {
public:
  // see iceModel.cc for implementation of constructor and destructor:
  IceModel(IceGrid::Ptr g, Context::Ptr context);
  virtual ~IceModel(); // must be virtual merely because some members are virtual

  // see iMinit.cc
  virtual void time_setup();

  IceGrid::Ptr grid() const;

  Context::Ptr ctx() const;

  virtual void allocate_submodels();
  virtual void allocate_stressbalance();
  virtual void allocate_age_model();
  virtual void allocate_bed_deformation();
  virtual void allocate_bedrock_thermal_unit();
  virtual void allocate_energy_model();
  virtual void allocate_subglacial_hydrology();
  virtual void allocate_basal_yield_stress();
  virtual void allocate_couplers();
  virtual void allocate_iceberg_remover();

  virtual void model_state_setup();
  virtual void misc_setup();
  virtual void init_diagnostics();
  virtual void init_calving();

  virtual void list_diagnostics();

  // see iceModel.cc
  void init();


  /** Run PISM in the "standalone" mode. */
  virtual void run();
  /** Advance the current PISM run to a specific time */
  virtual void run_to(double time);
  virtual void step(bool do_mass_continuity, bool do_skip);
  void reset_counters();

  // see iMbootstrap.cc
  virtual void bootstrap_2d(const PIO &input_file);

  // see iMoptions.cc
  virtual void setFromOptions();
  virtual std::set<std::string> output_size_from_option(const std::string &option,
                                                        const std::string &description,
                                                        const std::string &default_value);
  virtual std::set<std::string> set_output_size(const std::string &keyword);
  virtual std::string get_output_size(const std::string &option);

  // see iMutil.cc
  virtual void additionalAtStartTimestep();
  virtual void additionalAtEndTimestep();
  virtual void compute_cell_areas(); // is an initialization step; should go there

  // see iMIO.cc
  virtual void restart_2d(const PIO &input_file, unsigned int record);

  virtual void initialize_2d() __attribute__((noreturn));

  void initialize_cumulative_fluxes(const PIO &input_file);
  void reset_cumulative_fluxes();

  virtual void writeFiles(const std::string &default_filename);
  virtual void write_model_state(const PIO &nc);

  enum MetadataFlag {WRITE_MAPPING                         = 1,
                     WRITE_MAPPING_AND_RUN_STATS           = 1 | 2,
                     WRITE_MAPPING_AND_GLOBAL_ATTRIBUTES   = 1 | 4,
                     WRITE_RUN_STATS                       = 2,
                     WRITE_RUN_STATS_AND_GLOBAL_ATTRIBUTES = 2 | 4,
                     WRITE_GLOBAL_ATTRIBUTES               = 4,
                     WRITE_ALL                             = 1 | 2 | 4};

  virtual void write_metadata(const PIO &nc, MetadataFlag flag);
  virtual void write_diagnostics(const PIO &nc, const std::set<std::string> &vars,
                                 IO_Type nctype);
protected:

  //! Computational grid
  const IceGrid::Ptr m_grid;
  //! Configuration flags and parameters
  const Config::Ptr m_config;
  //! Execution context
  const Context::Ptr m_ctx;
  //! Unit system
  const units::System::Ptr m_sys;
  //! Logger
  const Logger::Ptr m_log;
  //! Time manager
  const Time::Ptr m_time;

  //! stores global attributes saved in a PISM output file
  VariableMetadata m_output_global_attributes;

  //! run statistics
  VariableMetadata m_run_stats;

  //! the list of sub-models, for writing model states and obtaining diagnostics
  std::map<std::string,const Component*> m_submodels;

  hydrology::Hydrology   *m_subglacial_hydrology;
  YieldStress *m_basal_yield_stress_model;

  energy::BedThermalUnit *m_btu;
  energy::EnergyModel *m_energy_model;

  AgeModel *m_age_model;

  calving::IcebergRemover     *m_iceberg_remover;
  calving::OceanKill          *m_ocean_kill_calving;
  calving::FloatKill          *m_float_kill_calving;
  calving::CalvingAtThickness *m_thickness_threshold_calving;
  calving::EigenCalving       *m_eigen_calving;
  calving::vonMisesCalving    *m_vonmises_calving;
  FrontalMelt                 *m_frontal_melt;

  surface::SurfaceModel *m_surface;
  ocean::OceanModel     *m_ocean;
  bed::BedDef           *m_beddef;

  // state variables and some diagnostics/internals

  //! ice surface elevation; ghosted
  IceModelVec2S m_ice_surface_elevation;
  //! ghosted
  IceModelVec2S m_ice_thickness;
  //! ghosted
  IceModelVec2S m_basal_yield_stress;
  //! rate of production of basal meltwater (ice-equivalent); no ghosts
  IceModelVec2S m_basal_melt_rate;
  //! Longitude; ghosted to compute cell areas
  IceModelVec2S m_longitude;
  //! Latitude; ghosted to compute cell areas
  IceModelVec2S m_latitude;
  //! accumulated mass advected to a partially filled grid cell
  IceModelVec2S m_Href;
  //! cell areas (computed using the WGS84 datum)
  IceModelVec2S m_cell_area;
  //! flux divergence
  IceModelVec2S m_flux_divergence;

  FractureFields *m_fracture;

protected:

  //! \brief mask for flow type with values ice_free_bedrock, grounded_ice, floating_ice,
  //! ice_free_ocean
  IceModelVec2CellType m_cell_type;

  //! mask to determine Dirichlet boundary locations
  IceModelVec2Int m_ssa_dirichlet_bc_mask;
  //! Dirichlet boundary velocities
  IceModelVec2V m_ssa_dirichlet_bc_values;
  
  //! mask to determine grounding line position
  IceModelVec2S m_gl_mask;

  // parameters
  //! mass continuity time step, s
  double m_dt;
  //! time of last update for enthalpy/temperature
  double t_TempAge;
  //! enthalpy/temperature and age time-steps
  double dt_TempAge;

  struct FluxCounters {
    FluxCounters();

    double H_to_Href;
    double Href_to_H;
    double discharge;
    double grounded_basal;
    double nonneg_rule;
    double sub_shelf;
    double sum_divQ_SIA;
    double sum_divQ_SSA;
    double surface;
  };

  struct FluxFields {
    FluxFields(IceGrid::ConstPtr grid);
    void reset();
    void regrid(const PIO &input_file);

    //! climatic_mass_balance
    IceModelVec2S climatic_mass_balance;
    //! grounded basal (melt/freeze-on) cumulative flux
    IceModelVec2S basal_grounded;
    //! floating (sub-shelf) basal (melt/freeze-on) cumulative flux
    IceModelVec2S basal_floating;
    //! cumulative nonnegative-rule flux
    IceModelVec2S nonneg;
    //! cumulative discharge (calving) flux
    IceModelVec2S discharge;
  };

protected:
  FluxCounters m_cumulative_fluxes;
  FluxFields m_cumulative_flux_fields;
  unsigned int m_skip_countdown;

  std::string m_adaptive_timestep_reason;

  std::string m_stdout_flags;

  // see iceModel.cc
  virtual void createVecs();

  virtual MaxTimestep max_timestep_diffusivity();
  virtual void max_timestep(double &dt_result, unsigned int &skip_counter);
  virtual unsigned int skip_counter(double input_dt, double input_dt_diffusivity);

  // see iMenergy.cc
  virtual void energyStep();

  virtual void combine_basal_melt_rate();

  // see iMgeometry.cc
  virtual void enforce_consistency_of_geometry();
  virtual void cell_interface_fluxes(bool dirichlet_bc,
                                     int i, int j,
                                     StarStencil<Vector2> input_velocity,
                                     StarStencil<double> input_flux,
                                     StarStencil<double> &output_velocity,
                                     StarStencil<double> &output_flux);
  virtual void adjust_flow(StarStencil<int> mask,
                           StarStencil<double> &SSA_velocity,
                           StarStencil<double> &SIA_flux);
  virtual void massContExplicitStep();
  virtual void update_grounded_cell_fraction();
  virtual void do_calving();
  virtual void Href_cleanup();
  virtual void update_cumulative_discharge(const IceModelVec2S &thickness,
                                           const IceModelVec2S &thickness_old,
                                           const IceModelVec2S &Href,
                                           const IceModelVec2S &Href_old);


  // see iMIO.cc
  virtual void dumpToFile(const std::string &filename);
  virtual void regrid(int dimensions);
  virtual void regrid_variables(const PIO &regrid_file,
                                const std::set<std::string> &regrid_vars,
                                unsigned int ndims);

  // see iMfractures.cc
  virtual void calculateFractureDensity();

  // see iMpartgrid.cc
  virtual void residual_redistribution(IceModelVec2S &residual);
  virtual void residual_redistribution_iteration(IceModelVec2S &residual, bool &done);

  // see iMreport.cc
  virtual double compute_temperate_base_fraction(double ice_area);
  virtual double compute_original_ice_fraction(double ice_volume);
  virtual void summary(bool tempAndAge);
  virtual void summaryPrintLine(bool printPrototype, bool tempAndAge,
                                double delta_t,
                                double volume, double area,
                                double meltfrac, double max_diffusivity);

public:

  // see iMreport.cc;  methods for computing diagnostic quantities:
  // scalar:
  double ice_volume(double thickness_threshold) const;
  double ice_volume_not_displacing_seawater(double thickness_threshold) const;
  double sealevel_volume(double thickness_threshold) const;
  double ice_volume_temperate(double thickness_threshold) const;
  double ice_volume_cold(double thickness_threshold) const;
  double ice_area(double thickness_threshold) const;
  double ice_area_grounded(double thickness_threshold) const;
  double ice_area_floating(double thickness_threshold) const;
  double ice_area_temperate(double thickness_threshold) const;
  double ice_area_cold(double thickness_threshold) const;

protected:
  // see iMutil.cc
  virtual int endOfTimeStepHook();
  virtual void stampHistoryCommand();
  virtual void stampHistoryEnd();
  virtual void stampHistory(const std::string &);
  virtual void update_run_stats();

protected:
  // working space (a convenience)
  static const int m_n_work2d = 4;
  mutable IceModelVec2S m_work2d[m_n_work2d];

  stressbalance::StressBalance *m_stress_balance;

public:
  const stressbalance::StressBalance* stress_balance() const;
  const ocean::OceanModel* ocean_model() const;
  const bed::BedDef* bed_model() const;
  const energy::BedThermalUnit* bedrock_thermal_model() const;
  const energy::EnergyModel* energy_balance_model() const;

  const IceModelVec2S& ice_thickness() const;
  const IceModelVec2S& ice_surface_elevation() const;
  const IceModelVec2CellType& cell_type() const;
  const IceModelVec2S &cell_area() const;

  FluxCounters cumulative_fluxes() const;
  const IceModelVec2S& flux_divergence() const;
  const FluxFields& cumulative_fluxes_2d() const;
  double dt() const;

protected:

  std::map<std::string,Diagnostic::Ptr> m_diagnostics;
  std::map<std::string,TSDiagnostic::Ptr> m_ts_diagnostics;

  // Set of variables to put in the output file:
  std::set<std::string> m_output_vars;

  // This is related to the snapshot saving feature
  std::string m_snapshots_filename;
  bool m_save_snapshots, m_snapshots_file_is_ready, m_split_snapshots;
  std::vector<double> m_snapshot_times;
  std::set<std::string> m_snapshot_vars;
  unsigned int m_current_snapshot;
  void init_snapshots();
  void write_snapshot();
  MaxTimestep save_max_timestep(double my_t);

  // scalar time-series
  bool m_save_ts;                 //! true if the user requested time-series output
  //! file to write time-series to
  std::string m_ts_filename;
  //! The history attribute in the -ts_file. Read from -ts_file if -ts_append is set, otherwise
  //! empty.
  std::string m_old_ts_file_history;
  std::vector<double> m_ts_times; //! times requested
  unsigned int m_current_ts;      //! index of the current time
  std::set<std::string> m_ts_vars;                //! variables requested
  void init_timeseries();
  void flush_timeseries();
  void write_timeseries();
  MaxTimestep ts_max_timestep(double my_t);

  // spatially-varying time-series
  bool m_save_extra, m_extra_file_is_ready, m_split_extra;
  std::string m_extra_filename;
  //! The history attribute in the -extra_file. Read from -extra_file if -extra_append is set,
  //! otherwise empty.
  std::string m_old_extra_file_history;
  std::vector<double> m_extra_times;
  unsigned int m_next_extra;
  double m_last_extra;
  std::set<std::string> m_extra_vars;
  TimeBoundsMetadata m_extra_bounds;
  TimeseriesMetadata m_timestamp;
  void init_extras();
  void write_extras();
  MaxTimestep extras_max_timestep(double my_t);

  // automatic backups
  double m_backup_interval;
  std::string m_backup_filename;
  double m_last_backup_time;
  std::set<std::string> m_backup_vars;
  void init_backups();
  void write_backup();

  // last time at which PISM hit a multiple of X years, see the
  // timestep_hit_multiples configuration parameter
  double m_timestep_hit_multiples_last_time;

  // diagnostic viewers; see iMviewers.cc
  virtual void init_viewers();
  virtual void update_viewers();
  virtual void view_field(const IceModelVec *field);
  std::set<std::string> m_map_viewers, m_slice_viewers;
  std::map<std::string,petsc::Viewer::Ptr> viewers;

private:
  double m_start_time;    // this is used in the wall-clock-time backup code
};

void check_minimum_ice_thickness(const IceModelVec2S &ice_thickness);
void check_maximum_ice_thickness(const IceModelVec2S &ice_thickness);

void bedrock_surface_temperature(double sea_level,
                                 const IceModelVec2CellType &cell_type,
                                 const IceModelVec2S &bed_topography,
                                 const IceModelVec2S &ice_thickness,
                                 const IceModelVec2S &basal_enthalpy,
                                 const IceModelVec2S &ice_surface_temperature,
                                 IceModelVec2S &result);

} // end of namespace pism

#endif /* __iceModel_hh */

