// Copyright (C) 2004-2015 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <signal.h>
#include <gsl/gsl_rng.h>
#include <petscsnes.h>
#include <petsctime.h>          // PetscGetTime()

#include "flowlaws.hh"


#include "pism_const.hh"
#include "iceModelVec.hh"
#include "PISMConfig.hh"
#include "PISMVars.hh"

namespace pism {

// forward declarations
class IceGrid;
class EnthalpyConverter;
class Hydrology;
class YieldStress;
class StressBalance;
class SurfaceModel;
class OceanModel;
class BedDef;
class BedThermalUnit;
class Diagnostic;
class TSDiagnostic;
class IcebergRemover;
class OceanKill;
class FloatKill;
class CalvingAtThickness;
class EigenCalving;


//! The base class for PISM.  Contains all essential variables, parameters, and flags for modelling an ice sheet.
class IceModel {
  // The following classes implement various diagnostic computations.
  // 2D and 3D:
  friend class IceModel_hardav;
  friend class IceModel_bwp;
  friend class IceModel_cts;
  friend class IceModel_dhdt;
  friend class IceModel_temp;
  friend class IceModel_temp_pa;
  friend class IceModel_temppabase;
  friend class IceModel_enthalpybase;
  friend class IceModel_enthalpysurf;
  friend class IceModel_tempbase;
  friend class IceModel_tempsurf;
  friend class IceModel_liqfrac;
  friend class IceModel_tempicethk;
  friend class IceModel_tempicethk_basal;
  friend class IceModel_new_mask;
  friend class IceModel_climatic_mass_balance_cumulative;
  friend class IceModel_dHdt;
  friend class IceModel_flux_divergence;
  // scalar:
  friend class IceModel_ivol;
  friend class IceModel_slvol;
  friend class IceModel_divoldt;
  friend class IceModel_iarea;
  friend class IceModel_imass;
  friend class IceModel_dimassdt;
  friend class IceModel_ivoltemp;
  friend class IceModel_ivolcold;
  friend class IceModel_ivolg;
  friend class IceModel_ivolf;
  friend class IceModel_iareatemp;
  friend class IceModel_iareacold;
  friend class IceModel_ienthalpy;
  friend class IceModel_iareag;
  friend class IceModel_iareaf;
  friend class IceModel_dt;
  friend class IceModel_max_diffusivity;
  friend class IceModel_surface_flux;
  friend class IceModel_surface_flux_cumulative;
  friend class IceModel_grounded_basal_flux;
  friend class IceModel_grounded_basal_flux_cumulative;
  friend class IceModel_sub_shelf_flux;
  friend class IceModel_sub_shelf_flux_cumulative;
  friend class IceModel_nonneg_flux;
  friend class IceModel_nonneg_flux_cumulative;
  friend class IceModel_discharge_flux;
  friend class IceModel_discharge_flux_cumulative;
  friend class IceModel_nonneg_flux_2D_cumulative;
  friend class IceModel_grounded_basal_flux_2D_cumulative;
  friend class IceModel_floating_basal_flux_2D_cumulative;
  friend class IceModel_discharge_flux_2D_cumulative;
  friend class IceModel_max_hor_vel;
  friend class IceModel_sum_divQ_flux;
  friend class IceModel_H_to_Href_flux;
  friend class IceModel_Href_to_H_flux;
public:
  // see iceModel.cc for implementation of constructor and destructor:
  IceModel(IceGrid &g, Config &config, Config &overrides);
  virtual ~IceModel(); // must be virtual merely because some members are virtual

  // see iMinit.cc
  virtual PetscErrorCode grid_setup();

  virtual PetscErrorCode allocate_submodels();
  virtual PetscErrorCode allocate_enthalpy_converter();
  virtual PetscErrorCode allocate_stressbalance();
  virtual PetscErrorCode allocate_bed_deformation();
  virtual PetscErrorCode allocate_bedrock_thermal_unit();
  virtual PetscErrorCode allocate_subglacial_hydrology();
  virtual PetscErrorCode allocate_basal_yield_stress();
  virtual PetscErrorCode allocate_couplers();
  virtual PetscErrorCode allocate_iceberg_remover();

  virtual PetscErrorCode init_couplers();
  virtual PetscErrorCode init_step_couplers();
  virtual PetscErrorCode set_grid_from_options();
  virtual PetscErrorCode set_grid_defaults();
  virtual PetscErrorCode model_state_setup();
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode allocate_internal_objects();
  virtual PetscErrorCode misc_setup();
  virtual PetscErrorCode init_diagnostics();
  virtual PetscErrorCode init_calving();

  virtual PetscErrorCode list_diagnostics();

  // see iceModel.cc
  PetscErrorCode init();


  /** Run PISM in the "standalone" mode. */
  virtual PetscErrorCode run();
  /** Advance the current PISM run to a specific time */
  virtual PetscErrorCode run_to(double time);
  virtual PetscErrorCode step(bool do_mass_continuity, bool do_energy, bool do_age, bool do_skip);
  virtual PetscErrorCode setExecName(const std::string &my_executable_short_name);
  virtual void reset_counters();

  // see iMbootstrap.cc 
  virtual PetscErrorCode bootstrapFromFile(const std::string &fname);
  virtual PetscErrorCode bootstrap_2d(const std::string &fname);
  virtual PetscErrorCode putTempAtDepth();

  // see iMoptions.cc
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode output_size_from_option(const std::string &option,
                                                 const std::string &description,
                                                 const std::string &default_value,
                                                 std::set<std::string> &result);
  virtual PetscErrorCode set_output_size(const std::string &keyword,
                                         std::set<std::string> &result);
  virtual std::string get_output_size(const std::string &option);

  // see iMutil.cc
  virtual PetscErrorCode additionalAtStartTimestep();
  virtual PetscErrorCode additionalAtEndTimestep();
  virtual PetscErrorCode compute_cell_areas(); // is an initialization step; should go there

  // see iMIO.cc
  virtual PetscErrorCode initFromFile(const std::string &name);
  virtual PetscErrorCode writeFiles(const std::string &default_filename);
  virtual PetscErrorCode write_model_state(const PIO &nc);
  virtual PetscErrorCode write_metadata(const PIO &nc,
                                        bool write_mapping,
                                        bool write_run_stats);
  virtual PetscErrorCode write_variables(const PIO &nc, const std::set<std::string> &vars,
                                         IO_Type nctype);
protected:

  IceGrid &grid;

  Config &config,           //!< configuration flags and parameters
    &overrides;                 //!< flags and parameters overriding config, see -config_override

  NCVariable global_attributes, //!< stores global attributes saved in a PISM output file
    mapping,                    //!< grid projection (mapping) parameters
    run_stats;                  //!< run statistics

  Hydrology   *subglacial_hydrology;
  YieldStress *basal_yield_stress_model;

  EnthalpyConverter *EC;
  BedThermalUnit *btu;

  IcebergRemover     *iceberg_remover;
  OceanKill          *ocean_kill_calving;
  FloatKill          *float_kill_calving;
  CalvingAtThickness *thickness_threshold_calving;
  EigenCalving       *eigen_calving;

  SurfaceModel *surface;
  OceanModel   *ocean;
  BedDef       *beddef;
  bool external_surface_model, external_ocean_model;

  //! \brief A dictionary with pointers to IceModelVecs below, for passing them
  //! from the IceModel core to other components (such as surface and ocean models)
  Vars variables;

  // state variables and some diagnostics/internals
  IceModelVec2S ice_surface_elevation,          //!< ice surface elevation; ghosted
    ice_thickness,              //!< ghosted
    basal_yield_stress,         //!< ghosted
    basal_melt_rate,           //!< rate of production of basal meltwater (ice-equivalent); no ghosts
    vLongitude, //!< Longitude; ghosted to compute cell areas
    vLatitude,  //!< Latitude; ghosted to compute cell areas
    bed_topography,             //!< bed topography; ghosted
    bed_uplift_rate,    //!< bed uplift rate; no ghosts
    geothermal_flux,   //!< geothermal flux; no ghosts
    vFD,    //!< fracture density
    vFG,    //!< fracture growth rate
    vFH,    //!< fracture healing rate
    vFE,    //!< fracture flow enhancement
    vFA,    //!< fracture age
    vFT,    //!< fracture toughness
    bedtoptemp,     //!< temperature seen by bedrock thermal layer, if present; no ghosts
    vHref,          //!< accumulated mass advected to a partially filled grid cell
    climatic_mass_balance,              //!< accumulation/ablation rate; no ghosts
    climatic_mass_balance_cumulative,    //!< cumulative climatic_mass_balance
    grounded_basal_flux_2D_cumulative, //!< grounded basal (melt/freeze-on) cumulative flux
    floating_basal_flux_2D_cumulative, //!< floating (sub-shelf) basal (melt/freeze-on) cumulative flux
    nonneg_flux_2D_cumulative,         //!< cumulative nonnegative-rule flux
    discharge_flux_2D_cumulative,      //!< cumulative discharge (calving) flux (2D field)
    ice_surface_temp,           //!< ice temperature at the ice surface but below firn; no ghosts
    liqfrac_surface,    //!< ice liquid water fraction at the top surface of the ice
    shelfbtemp,         //!< ice temperature at the shelf base; no ghosts
    shelfbmassflux,     //!< ice mass flux into the ocean at the shelf base; no ghosts
    cell_area,          //!< cell areas (computed using the WGS84 datum)
    flux_divergence;    //!< flux divergence

public:
  IceModelVec2S* get_geothermal_flux();
protected:

  IceModelVec2 strain_rates; //!< major and minor principal components of horizontal strain-rate tensor
  
  IceModelVec2 deviatoric_stresses; //!< components of horizontal stress tensor along axes and shear stress

  IceModelVec2Int vMask, //!< \brief mask for flow type with values ice_free_bedrock,
  //!< grounded_ice, floating_ice, ice_free_ocean
    vBCMask; //!< mask to determine Dirichlet boundary locations
 
  IceModelVec2V vBCvel; //!< Dirichlet boundary velocities
  
  IceModelVec2S gl_mask, //!< mask to determine grounding line position
    gl_mask_x, //!< mask to determine grounding line position in x-direction
    gl_mask_y; //!< mask to determine grounding line position in y-direction

  IceModelVec3
  T3,             //!< absolute temperature of ice; K (ghosted)
    Enth3,          //!< enthalpy; J / kg (ghosted)
    age3;           //!< age of ice; s (ghosted because it is averaged onto the staggered-grid)

  // parameters
  double   dt,     //!< mass continuity time step, s
    t_TempAge,  //!< time of last update for enthalpy/temperature
    dt_TempAge,  //!< enthalpy/temperature and age time-steps
    maxdt_temporary, dt_force,
    CFLviolcount,    //!< really is just a count, but GlobalSum requires this type
    dt_from_cfl, CFLmaxdt, CFLmaxdt2D,
    gmaxu, gmaxv, gmaxw,  // global maximums on 3D grid of abs value of vel components
    grounded_basal_ice_flux_cumulative,
    nonneg_rule_flux_cumulative,
    sub_shelf_ice_flux_cumulative,
    surface_ice_flux_cumulative,
    sum_divQ_SIA_cumulative,
    sum_divQ_SSA_cumulative,
    Href_to_H_flux_cumulative,
    H_to_Href_flux_cumulative,
    discharge_flux_cumulative;      //!< cumulative discharge (calving) flux

  unsigned int skipCountDown;

  // flags
  std::string m_adaptive_timestep_reason;

  std::string stdout_flags;

  std::string executable_short_name;
  
protected:
  // see iceModel.cc
  virtual PetscErrorCode createVecs();

  // see iMadaptive.cc
  virtual PetscErrorCode max_timestep_cfl_3d(double &dt_result);
  virtual PetscErrorCode max_timestep_cfl_2d(double &dt_result);
  virtual PetscErrorCode max_timestep_diffusivity(double &dt_result);
  virtual PetscErrorCode max_timestep(double &dt_result, unsigned int &skip_counter);
  virtual PetscErrorCode countCFLViolations(double* CFLviol);
  virtual unsigned int skip_counter(double input_dt, double input_dt_diffusivity);

  // see iMage.cc
  virtual PetscErrorCode ageStep();

  // see iMenergy.cc
  virtual PetscErrorCode energyStep();
  virtual PetscErrorCode get_bed_top_temp(IceModelVec2S &result);
  virtual bool checkThinNeigh(IceModelVec2S &thickness, int i, int j, const double threshold);

  virtual PetscErrorCode combine_basal_melt_rate();

  // see iMenthalpy.cc
  virtual PetscErrorCode compute_enthalpy_cold(IceModelVec3 &temperature, IceModelVec3 &result);
  virtual PetscErrorCode compute_enthalpy(IceModelVec3 &temperature, IceModelVec3 &liquid_water_fraction,
                                          IceModelVec3 &result);
  virtual PetscErrorCode compute_liquid_water_fraction(IceModelVec3 &enthalpy, IceModelVec3 &result);

  virtual PetscErrorCode setCTSFromEnthalpy(IceModelVec3 &result);

  virtual PetscErrorCode enthalpyAndDrainageStep(double* vertSacrCount,
                                                 double* liquifiedVol, double* bulgeCount);

  // see iMgeometry.cc
  virtual PetscErrorCode updateSurfaceElevationAndMask();
  virtual PetscErrorCode update_mask(IceModelVec2S &bed,
                                     IceModelVec2S &ice_thickness,
                                     IceModelVec2Int &mask);
  virtual PetscErrorCode update_surface_elevation(IceModelVec2S &bed,
                                                  IceModelVec2S &ice_thickness,
                                                  IceModelVec2S &result);
  virtual void cell_interface_fluxes(bool dirichlet_bc,
                                     int i, int j,
                                     planeStar<Vector2> input_velocity,
                                     planeStar<double> input_flux,
                                     planeStar<double> &output_velocity,
                                     planeStar<double> &output_flux);
  virtual void adjust_flow(planeStar<int> mask,
                           planeStar<double> &SSA_velocity,
                           planeStar<double> &SIA_flux);
  virtual PetscErrorCode massContExplicitStep();
  virtual PetscErrorCode update_floatation_mask();
  virtual PetscErrorCode do_calving();
  virtual PetscErrorCode Href_cleanup();
  virtual PetscErrorCode update_cumulative_discharge(IceModelVec2S &thickness,
                                                     IceModelVec2S &thickness_old,
                                                     IceModelVec2S &Href,
                                                     IceModelVec2S &Href_old);


  // see iMIO.cc
  virtual PetscErrorCode dumpToFile(const std::string &filename);
  virtual PetscErrorCode regrid(int dimensions);
  virtual PetscErrorCode regrid_variables(const std::string &filename, const std::set<std::string> &regrid_vars, unsigned int ndims);
  virtual PetscErrorCode init_enthalpy(const std::string &filename, bool regrid, int last_record);

  // see iMfractures.cc
  virtual PetscErrorCode calculateFractureDensity();

  // see iMpartgrid.cc
  double get_threshold_thickness(planeStar<int> Mask,
                                 planeStar<double> thickness,
                                 planeStar<double> surface_elevation,
                                 double bed_elevation,
                                 bool reduce_frontal_thickness);
  virtual PetscErrorCode residual_redistribution(IceModelVec2S &residual);
  virtual PetscErrorCode residual_redistribution_iteration(IceModelVec2S &residual, bool &done);

  // see iMreport.cc
  virtual PetscErrorCode energyStats(double iarea,double &gmeltfrac);
  virtual PetscErrorCode ageStats(double ivol, double &gorigfrac);
  virtual PetscErrorCode summary(bool tempAndAge);
  virtual PetscErrorCode summaryPrintLine(PetscBool printPrototype, bool tempAndAge,
                                          double delta_t,
                                          double volume, double area,
                                          double meltfrac, double max_diffusivity);

  // see iMreport.cc;  methods for computing diagnostic quantities:
  // scalar:
  virtual PetscErrorCode compute_ice_volume(double &result);
  virtual PetscErrorCode compute_sealevel_volume(double &result);
  virtual PetscErrorCode compute_ice_volume_temperate(double &result);
  virtual PetscErrorCode compute_ice_volume_cold(double &result);
  virtual PetscErrorCode compute_ice_area(double &result);
  virtual PetscErrorCode compute_ice_area_temperate(double &result);
  virtual PetscErrorCode compute_ice_area_cold(double &result);
  virtual PetscErrorCode compute_ice_area_grounded(double &result);
  virtual PetscErrorCode compute_ice_area_floating(double &result);
  virtual PetscErrorCode compute_ice_enthalpy(double &result);

  // see iMtemp.cc
  virtual PetscErrorCode excessToFromBasalMeltLayer(const double rho, const double c, const double L,
                                                    const double z, const double dz,
                                                    double *Texcess, double *bwat);
  virtual PetscErrorCode temperatureStep(double* vertSacrCount, double* bulgeCount);

  // see iMutil.cc
  virtual int            endOfTimeStepHook();
  virtual PetscErrorCode stampHistoryCommand();
  virtual PetscErrorCode stampHistoryEnd();
  virtual PetscErrorCode stampHistory(const std::string &);
  virtual PetscErrorCode update_run_stats();
  virtual PetscErrorCode check_maximum_thickness();
  virtual PetscErrorCode check_maximum_thickness_hook(const int old_Mz);

protected:
  // working space (a convenience)
  static const int nWork2d=3;
  IceModelVec2S vWork2d[nWork2d];
  IceModelVec2V vWork2dV;

  // 3D working space
  IceModelVec3 vWork3d;

  StressBalance *stress_balance;

public:
  StressBalance* get_stress_balance();
protected:

  std::map<std::string,Diagnostic*> diagnostics;
  std::map<std::string,TSDiagnostic*> ts_diagnostics;

  // Set of variables to put in the output file:
  std::set<std::string> output_vars;

  // This is related to the snapshot saving feature
  std::string snapshots_filename;
  bool save_snapshots, snapshots_file_is_ready, split_snapshots;
  std::vector<double> snapshot_times;
  std::set<std::string> snapshot_vars;
  unsigned int current_snapshot;
  PetscErrorCode init_snapshots();
  PetscErrorCode write_snapshot();

  // scalar time-series
  bool save_ts;                 //! true if the user requested time-series output
  std::string ts_filename;              //! file to write time-series to
  std::vector<double> ts_times; //! times requested
  unsigned int current_ts;      //! index of the current time
  std::set<std::string> ts_vars;                //! variables requested
  PetscErrorCode init_timeseries();
  PetscErrorCode flush_timeseries();
  PetscErrorCode write_timeseries();
  PetscErrorCode ts_max_timestep(double my_t, double& my_dt, bool &restrict);

  // spatially-varying time-series
  bool save_extra, extra_file_is_ready, split_extra;
  std::string extra_filename;
  std::vector<double> extra_times;
  unsigned int next_extra;
  double last_extra;
  std::set<std::string> extra_vars;
  NCTimeBounds extra_bounds;
  NCTimeseries timestamp;
  PetscErrorCode init_extras();
  PetscErrorCode write_extras();
  PetscErrorCode extras_max_timestep(double my_t, double& my_dt, bool &restrict);

  // automatic backups
  double backup_interval;
  std::string backup_filename;
  double last_backup_time;
  std::set<std::string> backup_vars;
  PetscErrorCode init_backups();
  PetscErrorCode write_backup();

  // last time at which PISM hit a multiple of X years, see the
  // timestep_hit_multiples configuration parameter
  double timestep_hit_multiples_last_time;

  // diagnostic viewers; see iMviewers.cc
  virtual PetscErrorCode init_viewers();
  virtual PetscErrorCode update_viewers();
  std::set<std::string> map_viewers, slice_viewers;
  int     id, jd;            // sounding indexes
  std::map<std::string,PetscViewer> viewers;

private:
  PetscLogDouble start_time;    // this is used in the wall-clock-time backup code
};

} // end of namespace pism

#endif /* __iceModel_hh */

