// Copyright (C) 2004-2012 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
#include <petsctime.h>		// PetscGetTime()

#include "flowlaws.hh"


#include "pism_const.hh"
#include "iceModelVec.hh"
#include "NCVariable.hh"
#include "PISMVars.hh"

// forward declarations
class IceGrid;
class EnthalpyConverter;
class PISMYieldStress;
class IceBasalResistancePlasticLaw;
class PISMStressBalance;
class PISMSurfaceModel;
class PISMOceanModel;
class PISMBedDef;
class PISMBedThermalUnit;
class PISMDiagnostic;
class PISMTSDiagnostic;

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond


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
  friend class IceModel_acab_cumulative;
  friend class IceModel_dHdt;
  // scalar:
  friend class IceModel_ivol;
  friend class IceModel_slvol;
  friend class IceModel_divoldt;
  friend class IceModel_iarea;
  friend class IceModel_imass;
  friend class IceModel_dimassdt;
  friend class IceModel_ivoltemp;
  friend class IceModel_ivoltempf;
  friend class IceModel_ivolcold;
  friend class IceModel_ivolcoldf;
  friend class IceModel_ivolg;
  friend class IceModel_ivolf;
  friend class IceModel_iareatemp;
  friend class IceModel_iareatempf;
  friend class IceModel_iareacold;
  friend class IceModel_iareacoldf;
  friend class IceModel_ienthalpy;
  friend class IceModel_iareag;
  friend class IceModel_iareaf;
  friend class IceModel_dt;
  friend class IceModel_max_diffusivity;
  friend class IceModel_surface_flux;
  friend class IceModel_cumulative_surface_flux;
  friend class IceModel_basal_flux;
  friend class IceModel_cumulative_basal_flux;
  friend class IceModel_sub_shelf_flux;
  friend class IceModel_cumulative_sub_shelf_flux;
  friend class IceModel_nonneg_flux;
  friend class IceModel_cumulative_nonneg_flux;
  friend class IceModel_ocean_kill_flux;
  friend class IceModel_cumulative_ocean_kill_flux;
  friend class IceModel_float_kill_flux;
  friend class IceModel_cumulative_float_kill_flux;
  friend class IceModel_discharge_flux;
  friend class IceModel_cumulative_discharge_flux;
  friend class IceModel_max_hor_vel;
public:
  // see iceModel.cc for implementation of constructor and destructor:
  IceModel(IceGrid &g, NCConfigVariable &config, NCConfigVariable &overrides);
  virtual ~IceModel(); // must be virtual merely because some members are virtual

  // see iMinit.cc
  virtual PetscErrorCode grid_setup();

  virtual PetscErrorCode allocate_submodels();
  virtual PetscErrorCode allocate_flowlaw();
  virtual PetscErrorCode allocate_enthalpy_converter();
  virtual PetscErrorCode allocate_basal_resistance_law();
  virtual PetscErrorCode allocate_stressbalance();
  virtual PetscErrorCode allocate_bed_deformation();
  virtual PetscErrorCode allocate_bedrock_thermal_unit();
  virtual PetscErrorCode allocate_basal_yield_stress();

  virtual PetscErrorCode init_couplers();
  virtual PetscErrorCode set_grid_from_options();
  virtual PetscErrorCode set_grid_defaults();
  virtual PetscErrorCode model_state_setup();
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode allocate_internal_objects();
  virtual PetscErrorCode misc_setup();
  virtual PetscErrorCode init_diagnostics();
  virtual PetscErrorCode init_ocean_kill();

  virtual PetscErrorCode list_diagnostics();

  // see iceModel.cc
  PetscErrorCode init();
  virtual PetscErrorCode run();
  virtual PetscErrorCode step(bool do_mass_continuity, bool do_energy, bool do_age, bool do_skip);
  virtual PetscErrorCode setExecName(string my_executable_short_name);
  virtual void reset_counters();

  // see iMbootstrap.cc 
  virtual PetscErrorCode bootstrapFromFile(string fname);
  virtual PetscErrorCode bootstrap_2d(string fname);
  virtual PetscErrorCode bootstrap_3d();
  virtual PetscErrorCode putTempAtDepth();

  // see iMoptions.cc
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode set_output_size(string option, string description,
					 string default_value, set<string> &result);
  virtual string         get_output_size(string option);

  // see iMutil.cc
  virtual void attach_surface_model(PISMSurfaceModel *surf);
  virtual void attach_ocean_model(PISMOceanModel *ocean);
  virtual PetscErrorCode additionalAtStartTimestep();
  virtual PetscErrorCode additionalAtEndTimestep();
  virtual PetscErrorCode compute_cell_areas(); // is an initialization step; should go there

  // see iMIO.cc
  virtual PetscErrorCode initFromFile(string);
  virtual PetscErrorCode writeFiles(string default_filename);
  virtual PetscErrorCode write_model_state(string filename);
  virtual PetscErrorCode write_metadata(string filename, bool write_mapping = true);
  virtual PetscErrorCode write_variables(string filename, set<string> vars,
					 PISM_IO_Type nctype);
protected:

  IceGrid               &grid;

  NCConfigVariable      mapping, //!< grid projection (mapping) parameters
    &config,			 //!< configuration flags and parameters
    &overrides;			 //!< flags and parameters overriding config, see -config_override
  NCGlobalAttributes    global_attributes;

  PISMYieldStress *basal_yield_stress;
  IceBasalResistancePlasticLaw *basal;

  EnthalpyConverter *EC;
  PISMBedThermalUnit *btu;

  PISMSurfaceModel *surface;
  PISMOceanModel   *ocean;
  PISMBedDef       *beddef;

  //! \brief A dictionary with pointers to IceModelVecs below, for passing them
  //! from the IceModel core to other components (such as surface and ocean models)
  PISMVars variables;

  // state variables and some diagnostics/internals
  IceModelVec2S
        vh,		//!< ice surface elevation; ghosted
        vH,		//!< ice thickness; ghosted
        vtauc,		//!< yield stress for basal till (plastic or pseudo-plastic model); ghosted
        vbwat,		//!< thickness of the basal meltwater; ghosted
        vbmr,           //!< rate of production of basal meltwater (ice-equivalent); no ghosts
        vLongitude,	//!< Longitude; ghosted to compute cell areas
        vLatitude,	//!< Latitude; ghosted to compute cell areas
        vbed,		//!< bed topography; ghosted
        vuplift,	//!< bed uplift rate; no ghosts
        vGhf,		//!< geothermal flux; no ghosts
        bedtoptemp,     //!< temperature seen by bedrock thermal layer, if present; no ghosts
                        //!< ghosted to be able to compute tauc "redundantly"

        vHref,          //!< accumulated mass advected to a partially filled grid cell
        vHresidual,     //!< residual ice mass of a not any longer partially (fully) filled grid cell
        vPrinStrain1,   //!< major principal component of horizontal strain-rate tensor
        vPrinStrain2,   //!< minor principal component of horizontal strain-rate tensor

    acab,		//!< accumulation/ablation rate; no ghosts
    acab_cumulative,    //!< cumulative acab
    artm,		//!< ice temperature at the ice surface but below firn; no ghosts
    liqfrac_surface,    //!< ice liquid water fraction at the top surface of the ice
    shelfbtemp,		//!< ice temperature at the shelf base; no ghosts
    shelfbmassflux,	//!< ice mass flux into the ocean at the shelf base; no ghosts
    cell_area;		//!< cell areas (computed using the WGS84 datum)

	
 
  IceModelVec2Int vMask, //!< \brief mask for flow type with values ice_free_bedrock,
                         //!< grounded_ice, floating_ice, ice_free_ocean
    ocean_kill_mask,     //!< mask used by the -ocean_kill code 
    vIcebergMask, //!< mask for iceberg identification

    vBCMask; //!< mask to determine Dirichlet boundary locations
 
  IceModelVec2V vBCvel; //!< Dirichlet boundary velocities


  IceModelVec3
        T3,		//!< absolute temperature of ice; K (ghosted)
        Enth3,          //!< enthalpy; J / kg (ghosted)
        tau3;		//!< age of ice; s (ghosted because it is averaged onto the staggered-grid)

  // parameters
  PetscReal   dt,     //!< mass continuity time step, s
              t_TempAge,  //!< time of last update for enthalpy/temperature
              dt_TempAge,  //!< enthalpy/temperature and age time-steps
              maxdt_temporary, dt_force,
              CFLviolcount,    //!< really is just a count, but PISMGlobalSum requires this type
              dt_from_diffus, dt_from_cfl, CFLmaxdt, CFLmaxdt2D, dt_from_eigencalving,
              gDmax,		// global max of the diffusivity
              gmaxu, gmaxv, gmaxw,  // global maximums on 3D grid of abs value of vel components
    cumulative_basal_ice_flux,
    cumulative_float_kill_flux,
    cumulative_discharge_flux,
    cumulative_nonneg_rule_flux,
    cumulative_ocean_kill_flux,
    cumulative_sub_shelf_ice_flux,
    cumulative_surface_ice_flux;
  PetscInt    skipCountDown;

  // physical parameters used frequently enough to make looking up via
  // config.get() a hassle; initialized in the IceModel constructor from the
  // configuration file; SHOULD NOT be hard-wired.
  PetscScalar standard_gravity;

  // flags
  PetscBool  shelvesDragToo, allowAboveMelting;
  PetscBool  repeatRedist, putOnTop;
  char        adaptReasonFlag;

  string      stdout_flags, stdout_ssa;

  string executable_short_name;
  
protected:
  // see iceModel.cc
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode deallocate_internal_objects();

  // see iMadaptive.cc
  virtual PetscErrorCode computeMax3DVelocities();
  virtual PetscErrorCode computeMax2DSlidingSpeed();
  virtual PetscErrorCode adaptTimeStepDiffusivity();
  virtual PetscErrorCode determineTimeStep(const bool doTemperatureCFL);
  virtual PetscErrorCode countCFLViolations(PetscScalar* CFLviol);

  // see iMage.cc
  virtual PetscErrorCode ageStep();

  // see iMcalving.cc
  virtual PetscErrorCode eigenCalving();
  virtual PetscErrorCode calvingAtThickness();
  virtual PetscErrorCode dt_from_eigenCalving();

  // see iMenergy.cc
  virtual PetscErrorCode energyStep();
  virtual PetscErrorCode get_bed_top_temp(IceModelVec2S &result);
  virtual bool checkThinNeigh(
       PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
       PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE);

  // see iMenthalpy.cc
  virtual PetscErrorCode compute_enthalpy_cold(IceModelVec3 &temperature, IceModelVec3 &result);
  virtual PetscErrorCode compute_enthalpy(IceModelVec3 &temperature, IceModelVec3 &liquid_water_fraction,
                                          IceModelVec3 &result);
  virtual PetscErrorCode compute_liquid_water_fraction(IceModelVec3 &enthalpy, IceModelVec3 &result);

  virtual PetscErrorCode setCTSFromEnthalpy(IceModelVec3 &useForCTS);

  virtual PetscErrorCode getEnthalpyCTSColumn(PetscScalar p_air, //!< atmospheric pressure
					      PetscScalar thk,	 //!< ice thickness
					      PetscInt ks,	 //!< index of the level just below the surface
					      PetscScalar **Enth_s //!< enthalpy of pressure-melting temperature cold ice
					      );

  virtual PetscErrorCode getlambdaColumn(PetscInt ks,	       //!< index of the level just below the surface
					 PetscScalar ice_rho_c,//!< default value only
                                         PetscScalar ice_k,    //!< default value only
					 const PetscScalar *Enth,   //!< enthalpy in the column
					 const PetscScalar *Enth_s, //!< enthalpy of pressure-melting temperature cold ice
					 const PetscScalar *w, //!< vert. velocity
					 PetscScalar *lambda //!< constant controlling choice of implicit method
					 );

  virtual PetscErrorCode enthalpyAndDrainageStep(
                PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                PetscScalar* bulgeCount);

  // see iMgeometry.cc
  virtual PetscErrorCode updateSurfaceElevationAndMask();
  virtual PetscErrorCode update_mask();
  virtual PetscErrorCode update_surface_elevation();
  virtual PetscErrorCode cell_interface_diffusive_flux(IceModelVec2Stag &Qstag, int i, int j,
                                                       planeStar<PetscScalar> &Q_output);
  virtual PetscErrorCode massContExplicitStep();

  // see iMhydrology.cc
  virtual PetscErrorCode diffuse_bwat();

  // see iMicebergs.cc
  virtual PetscErrorCode killIceBergs();           // call this one to do proper sequence
  virtual PetscErrorCode findIceBergCandidates();
  virtual PetscErrorCode identifyNotAnIceBerg();
  virtual PetscErrorCode killIdentifiedIceBergs();
  virtual PetscErrorCode killEasyIceBergs();       // FIXME: do we want this one to happen even if eigencalving does not happen?  should we be calling this one before any time that principal values need to be computed?

  // see iMIO.cc
  virtual PetscErrorCode dumpToFile(string filename);
  virtual PetscErrorCode regrid(int dimensions);
  virtual PetscErrorCode regrid_variables(string filename, set<string> regrid_vars, int ndims);

  // see iMpartgrid.cc
  virtual PetscErrorCode cell_interface_velocities(bool do_part_grid,
                                                   int i, int j,
                                                   planeStar<PetscScalar> &vel_output);
  PetscReal get_average_thickness(bool do_redist, planeStar<int> M,
                                  planeStar<PetscScalar> H);
  virtual PetscErrorCode redistResiduals();
  virtual PetscErrorCode calculateRedistResiduals();

  // see iMreport.cc
  virtual PetscErrorCode volumeArea(
                       PetscScalar& gvolume,PetscScalar& garea);
  virtual PetscErrorCode energyStats(
                       PetscScalar iarea,PetscScalar &gmeltfrac);
  virtual PetscErrorCode ageStats(PetscScalar ivol, PetscScalar &gorigfrac);
  virtual PetscErrorCode summary(bool tempAndAge);
  virtual PetscErrorCode summaryPrintLine(
              PetscBool printPrototype, bool tempAndAge,
              string date, PetscScalar delta_t, 
              PetscScalar volume, PetscScalar area,
              PetscScalar meltfrac, PetscScalar max_diffusivity);

  // see iMreport.cc;  methods for computing diagnostic quantities:
  // scalar:
  virtual PetscErrorCode compute_ice_volume(PetscScalar &result);
  virtual PetscErrorCode compute_sealevel_volume(PetscScalar &result);
  virtual PetscErrorCode compute_ice_volume_temperate(PetscScalar &result);
  virtual PetscErrorCode compute_ice_volume_cold(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area_temperate(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area_cold(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area_grounded(PetscScalar &result);
  virtual PetscErrorCode compute_ice_area_floating(PetscScalar &result);
  virtual PetscErrorCode compute_ice_enthalpy(PetscScalar &result);

  // see iMtemp.cc
  virtual PetscErrorCode excessToFromBasalMeltLayer(
                      const PetscScalar rho, const PetscScalar c, const PetscScalar L,
                      const PetscScalar z, const PetscScalar dz,
                      PetscScalar *Texcess, PetscScalar *bwat);
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);

  // see iMutil.cc
  virtual int            endOfTimeStepHook();
  virtual PetscErrorCode stampHistoryCommand();
  virtual PetscErrorCode stampHistoryEnd();
  virtual PetscErrorCode stampHistory(string);
  virtual PetscErrorCode check_maximum_thickness();
  virtual PetscErrorCode check_maximum_thickness_hook(const int old_Mz);
  virtual bool           issounding(const PetscInt i, const PetscInt j);

protected:
  // working space (a convenience)
  static const PetscInt nWork2d=2;
  IceModelVec2S vWork2d[nWork2d];
  IceModelVec2V vWork2dV;

  // 3D working space
  IceModelVec3 vWork3d;

  PISMStressBalance *stress_balance;

  map<string,PISMDiagnostic*> diagnostics;
  map<string,PISMTSDiagnostic*> ts_diagnostics;

  // Set of variables to put in the output file:
  set<string> output_vars;

  // This is related to the snapshot saving feature
  string snapshots_filename;
  bool save_snapshots, snapshots_file_is_ready, split_snapshots;
  vector<double> snapshot_times;
  set<string> snapshot_vars;
  unsigned int current_snapshot;
  PetscErrorCode init_snapshots();
  PetscErrorCode write_snapshot();

  // scalar time-series
  bool save_ts;			//! true if the user requested time-series output
  string ts_filename;		//! file to write time-series to
  vector<double> ts_times;	//! times requested
  unsigned int current_ts;	//! index of the current time
  set<string> ts_vars;		//! variables requested
  PetscErrorCode init_timeseries();
  PetscErrorCode flush_timeseries();
  PetscErrorCode write_timeseries();
  PetscErrorCode ts_max_timestep(double my_t, double& my_dt, bool &restrict);

  // spatially-varying time-series
  bool save_extra, extra_file_is_ready, split_extra;
  string extra_filename;
  vector<double> extra_times;
  unsigned int next_extra;
  double last_extra;
  set<string> extra_vars;
  NCTimeBounds extra_bounds;
  NCTimeseries timestamp;
  PetscErrorCode init_extras();
  PetscErrorCode write_extras();
  PetscErrorCode extras_max_timestep(double my_t, double& my_dt, bool &restrict);

  // automatic backups
  double backup_interval;
  string backup_filename;
  PetscReal last_backup_time;
  set<string> backup_vars;
  PetscErrorCode init_backups();
  PetscErrorCode write_backup();

  // diagnostic viewers; see iMviewers.cc
  virtual PetscErrorCode init_viewers();
  virtual PetscErrorCode update_viewers();
  set<string> map_viewers, slice_viewers, sounding_viewers;
  PetscInt     id, jd;	     // sounding indices
  map<string,PetscViewer> viewers;

private:
  PetscLogDouble start_time;    // this is used in the wall-clock-time backup code

  int event_step,		//!< total time spent doing time-stepping
    event_velocity,		//!< total velocity computation
    event_energy,		//!< energy balance computation
    event_mass,			//!< mass continuity computation
    event_age,			//!< age computation
    event_beddef,		//!< bed deformation step
    event_output,		//!< time spent writing the output file
    event_output_define,        //!< time spent defining variables
    event_snapshots,            //!< time spent writing snapshots
    event_backups;              //!< time spent writing backups files
};

#endif /* __iceModel_hh */

