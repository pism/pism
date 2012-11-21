// Copyright (C) 2007-2012 Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <petsc.h>
#include "icePSTexModel.hh"
#include "Timeseries.hh"
#include "SSA.hh"
#include "PISMStressBalance.hh"
#include "PIO.hh"
#include "pism_options.hh"

const PetscScalar 
  DEFAULT_PHI_STRONG = 15.0,  // till friction angle outside of stream
  stream_offset = 100.0e3,    // distance from grid center to start of stream
  stream_length = 650.0e3,    // length used in computing slope
  stream_change = 400.0e3,    // distance along stream at which  change from
                              //     'upstream' to 'downstream' till friction
                              //     angle occurs
  xi_slop = 0.15,             // fractions by which reduced till friction
  eta_slop = 0.5;             //     angle extends outside of stream; see
                              //     setTillPhi() and inStreamNbhd()

const int Nexpers = 6,
          NAME_LENGTH = 10;

struct ExperDescription {
  char         name[NAME_LENGTH];
  PetscScalar  bed_end_depth[4];   // drop for stream (see stream_length)
  PetscScalar  stream_width[4];    // width in km
  PetscScalar  upstream_phi[4];    // till friction angle in degrees
  PetscScalar  downstream_phi[4];
};

ExperDescription e[Nexpers] = {

  { "P0A",
    {0.0, 0.0, 0.0, 0.0}, 
    {70.0, 70.0, 70.0, 70.0}, // INACTIVE
    {5.0, 5.0, 5.0, 5.0},     // INACTIVE
    {5.0, 5.0, 5.0, 5.0}      // INACTIVE
  },

  { "P0I",
    {0.0, 500.0, 1000.0, 2000.0},
    {70.0, 70.0, 70.0, 70.0}, // INACTIVE
    {5.0, 5.0, 5.0, 5.0},     // INACTIVE
    {5.0, 5.0, 5.0, 5.0}      // INACTIVE
  },

  { "P1",
    {0.0, 0.0, 0.0, 0.0},
    {70.0, 30.0, 100.0, 50.0},
    {5.0, 5.0, 5.0, 5.0},
    {5.0, 5.0, 5.0, 5.0}
  },

  { "P2",                     // only three streams, variable grid orient
    {0.0, 0.0, 0.0, 0.0},     // last one INACTIVE
    {70.0, 70.0, 70.0, 70.0}, // last *three* INACTIVE; constant width
    {5.0, 5.0, 5.0, 5.0},     // last one INACTIVE
    {5.0, 5.0, 5.0, 5.0}      // last one INACTIVE
  },

  { "P3",
    {0.0, 500.0, 1000.0, 2000.0},
    {70.0, 70.0, 70.0, 70.0},
    {5.0, 5.0, 5.0, 5.0},
    {5.0, 5.0, 5.0, 5.0}
  },

  { "P4",
    {0.0, 0.0, 0.0, 0.0},
    {70.0, 70.0, 70.0, 70.0},
    {5.0, 5.0, 5.0, 5.0},
    {2.0, 3.0, 8.0, 10.0}
  }

};

PetscScalar stream_angle_P2[3] = {0.0, 100.0, 225.0};  // degrees

//! Say whether we are in the stream (strictly or not), and give local coords.
static bool inStreamNbhd(bool strictly_in_stream,
                         const PetscScalar angle, const PetscScalar width,
                         const PetscScalar x, const PetscScalar y,
                         PetscScalar &x_loc, PetscScalar &y_loc) {
  const PetscScalar sinrot = sin(angle),
                    cosrot = cos(angle);
  x_loc =  cosrot * x + sinrot * y - stream_offset;
  y_loc = -sinrot * x + cosrot * y;
  if (strictly_in_stream) {
    //if ( (x_loc > 0.0) && (fabs(y_loc) < width / 2.0) )
    if ( (x_loc > -1.0) && (fabs(y_loc) < width / 2.0) )
      return true;
    else
      return false;
  } else {
    //if ( (x_loc > - xi_slop * stream_change)
    if ( (x_loc > - xi_slop * stream_change - 1.0)
         && (fabs(y_loc) < (1.0 + eta_slop) * (width / 2.0)) )
      return true;
    else
      return false;
  }
}


//! Say whether we are strictly in the stream, and give local coords.
static bool inStream(const PetscScalar angle, const PetscScalar width,
                            const PetscScalar x, const PetscScalar y,
                            PetscScalar &x_loc, PetscScalar &y_loc) {
  return inStreamNbhd(true, angle,width, x,y, x_loc,y_loc);
}

IcePSTexModel::IcePSTexModel(IceGrid &g, NCConfigVariable &conf, NCConfigVariable &conf_overrides)
  : IceEISModel(g, conf, conf_overrides) {  // do almost nothing; derived need constructors
  expername = 'A';      // PST expers are closest to EISMINT II exper A
  ivol = PETSC_NULL;    // if we write series then this is not NULL
}


IcePSTexModel::~IcePSTexModel() {

  if (ivol != PETSC_NULL) {
    delete dt_ser;
    delete ivol;
    delete iarea;
    delete maxcbar;

    delete avup0;  delete avup1;  delete avup2;
    delete avdwn0; delete avdwn1; delete avdwn2;
    if (exper_chosen != 3) {
      delete avup3;
      delete avdwn3;
    }
  }
}


PetscErrorCode IcePSTexModel::prepare_series() {
  PetscErrorCode      ierr;

  // set up the file with name seriesname
  char outname[PETSC_MAX_PATH_LEN];
  PetscBool o_set;
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", outname, PETSC_MAX_PATH_LEN, &o_set); CHKERRQ(ierr);
  if (!o_set)
    strncpy(outname, "unnamedpst.nc", PETSC_MAX_PATH_LEN);
  strncpy(seriesname,"ser_pst_", PETSC_MAX_PATH_LEN);
  strcat(seriesname,outname);

  string time_dimension_name = config.get_string("time_dimension_name");

  ierr = verbPrintf(2,grid.com, 
    "  will write time series with special PST information to %s ...\n",
    seriesname); CHKERRQ(ierr);
  PIO nc(grid, grid.config.get_string("output_format"));
  ierr = nc.open(seriesname, PISM_WRITE); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  // set-up each scalar time series
  dt_ser = new DiagnosticTimeseries(&grid, "dt", time_dimension_name);
  dt_ser->set_units("s", "years");
  dt_ser->set_dimension_units("years", "");
  dt_ser->output_filename = seriesname;
  dt_ser->set_attr("long_name", "PST result: mass continuity time-step");
  dt_ser->set_attr("valid_min", 0.0);
  
  ivol = new DiagnosticTimeseries(&grid, "ivol", time_dimension_name);
  ivol->set_units("m3", "");
  ivol->set_dimension_units("years", "");
  ivol->output_filename = seriesname;
  ivol->set_attr("long_name", "PST result: total ice volume");
  
  iarea = new DiagnosticTimeseries(&grid, "iarea", time_dimension_name);
  iarea->set_units("m2", "");
  iarea->set_dimension_units("years", "");
  iarea->output_filename = seriesname;
  iarea->set_attr("long_name", "PST result: total ice area");
  
  maxcbar = new DiagnosticTimeseries(&grid, "maxcbar", time_dimension_name);
  maxcbar->set_units("m s-1", "m a-1");
  maxcbar->set_dimension_units("years", "");
  maxcbar->output_filename = seriesname;
  maxcbar->set_attr("long_name", 
    "PST result: maximum vertically-averaged ice speed anywhere on ice sheet");
  
  avup0 = new DiagnosticTimeseries(&grid, "avup0", time_dimension_name);
  avup0->set_units("m s-1", "m a-1");
  avup0->set_dimension_units("years", "");
  avup0->output_filename = seriesname;
  avup0->set_attr("long_name",
    "PST result: average speed in upper part of stream 0");
  
  avup1 = new DiagnosticTimeseries(&grid, "avup1", time_dimension_name);
  avup1->set_units("m s-1", "m a-1");
  avup1->set_dimension_units("years", "");
  avup1->output_filename = seriesname;
  avup1->set_attr("long_name",
    "PST result: average speed in upper part of stream 1");
  
  avup2 = new DiagnosticTimeseries(&grid, "avup2", time_dimension_name);
  avup2->set_units("m s-1", "m a-1");
  avup2->set_dimension_units("years", "");
  avup2->output_filename = seriesname;
  avup2->set_attr("long_name",
    "PST result: average speed in upper part of stream 2");
  
  avdwn0 = new DiagnosticTimeseries(&grid, "avdwn0", time_dimension_name);
  avdwn0->set_units("m s-1", "m a-1");
  avdwn0->set_dimension_units("years", "");
  avdwn0->output_filename = seriesname;
  avdwn0->set_attr("long_name",
    "PST result: average speed in downstream part of stream 0");
  
  avdwn1 = new DiagnosticTimeseries(&grid, "avdwn1", time_dimension_name);
  avdwn1->set_units("m s-1", "m a-1");
  avdwn1->set_dimension_units("years", "");
  avdwn1->output_filename = seriesname;
  avdwn1->set_attr("long_name",
    "PST result: average speed in downstream part of stream 1");
  
  avdwn2 = new DiagnosticTimeseries(&grid, "avdwn2", time_dimension_name);
  avdwn2->set_units("m s-1", "m a-1");
  avdwn2->set_dimension_units("years", "");
  avdwn2->output_filename = seriesname;
  avdwn2->set_attr("long_name",
    "PST result: average speed in downstream part of stream 2");

  if (exper_chosen != 3) {
    avup3 = new DiagnosticTimeseries(&grid, "avup3", time_dimension_name);
    avup3->set_units("m s-1", "m a-1");
    avup3->set_dimension_units("years", "");
    avup3->output_filename = seriesname;
    avup3->set_attr("long_name",
      "PST result: average speed in upper part of stream 3");

    avdwn3 = new DiagnosticTimeseries(&grid, "avdwn3", time_dimension_name);
    avdwn3->set_units("m s-1", "m a-1");
    avdwn3->set_dimension_units("years", "");
    avdwn3->output_filename = seriesname;
    avdwn3->set_attr("long_name",
      "PST result: average speed in downstream part of stream 3");
  }
  
  return 0;
}


PetscErrorCode IcePSTexModel::setFromOptions() {
  PetscErrorCode      ierr;

  ierr = IceEISModel::setFromOptions();  CHKERRQ(ierr);

  exper_chosen = -1;
  for (int j=0; j<Nexpers; j++) {
    bool  optionset;
    char optionname[20] = "-";
    strcat(optionname,e[j].name);
    ierr = PISMOptionsIsSet(optionname, optionset); CHKERRQ(ierr);
    if (optionset == PETSC_TRUE) {
      if (exper_chosen >= 0) {
        ierr = PetscPrintf(grid.com,
          "IcePSTexModel ERROR:  Only one experiment name option allowed.\n");
          CHKERRQ(ierr);
	PISMEnd();
      } else {
        exper_chosen = j;
        strcpy(exper_chosen_name,e[j].name);
      }
    }
  }
  if (exper_chosen < 0) {
    ierr = PetscPrintf(grid.com,
      "IcePSTexModel ERROR:  Unrecognized experiment name.\n"
      "  An experiment name option like '-P2' must be chosen.\n");
      CHKERRQ(ierr);
    PISMEnd();
  }
  
  ierr = verbPrintf(2,grid.com, 
    "setting up PST (Plastic till Stream w Thermocoupling) experiment %s ...\n",
    exper_chosen_name); CHKERRQ(ierr);

  bool optionset;
  ierr = PISMOptionsIsSet("-no_ser", optionset); CHKERRQ(ierr);
  if (optionset == PETSC_FALSE) {
    ierr = prepare_series(); CHKERRQ(ierr);
  }

  config.set("default_till_phi", DEFAULT_PHI_STRONG); 
  config.set_flag("include_bmr_in_continuity", true);
  config.set_flag("use_eta_transformation", true);

  if (exper_chosen <= 1) { // P0A and P0I are nonsliding SIA
    config.set_flag("use_ssa_velocity", false);
    config.set_flag("use_ssa_when_grounded", false);
  } else {
    // these options equiv to "-ssa -super -plastic"
    config.set_flag("use_ssa_velocity", true);
    config.set_flag("use_ssa_when_grounded", true);
  }  

  config.set_flag("do_diffuse_bwat", true);

  return 0;
}

PetscErrorCode IcePSTexModel::allocate_basal_yield_stress() {

  if (basal_yield_stress != NULL)
    return 0;

  basal_yield_stress = new PSTYieldStress(grid, config, exper_chosen, exper_chosen_name);
  
  return 0;
}

PetscErrorCode IcePSTexModel::allocate_stressbalance() {
  PetscErrorCode ierr;

  ierr = IceModel::allocate_stressbalance(); CHKERRQ(ierr);

  // typical strain rate is 100 m/yr per 100km in an ice shelf or fast ice stream
  const PetscScalar TYPICAL_STRAIN_RATE = (100.0 / secpera) / (100.0 * 1.0e3);
  const PetscScalar H_SSA_EXTENSION = 50.0; // m; thickness of ice shelf extension
  const PetscScalar constantHardnessForSSA = 1.9e8;  // Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
  const PetscScalar PSTconstantNuHForSSA = H_SSA_EXTENSION * constantHardnessForSSA
                     / (2.0 * pow(TYPICAL_STRAIN_RATE,2./3.)); // Pa s m

  // ssa == NULL means that the user chose non-sliding SIA
  SSA *ssa = dynamic_cast<SSA*>(stress_balance->get_stressbalance());
  if (ssa != NULL)
    ssa->strength_extension->set_notional_strength(PSTconstantNuHForSSA);
  
  return 0;
}

PetscErrorCode IcePSTexModel::initFromFile(string fname) {
  PetscErrorCode      ierr;

  ierr = IceEISModel::initFromFile(fname); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, 
                    "starting PST (Plastic till Stream w Thermocoupling) experiment %s from file  %s ...\n",
                    exper_chosen_name, fname.c_str()); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "  values of mask and phi = (till friction angle) in file will be ignored ...\n");
  CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "  bed topography from file is kept ...\n"); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IcePSTexModel::set_vars_from_options() {
  PetscErrorCode      ierr;

  ierr = IceEISModel::set_vars_from_options(); CHKERRQ(ierr);  

  ierr = verbPrintf(2,grid.com, 
    "setting variables for PST experiment %s ...\n", exper_chosen_name); CHKERRQ(ierr);

  ierr = setBedElev(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
    "  bed topography set for PST exper '%s' ...\n", exper_chosen_name); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IcePSTexModel::setBedElev() {
  PetscErrorCode ierr;
  
  const PetscScalar    width = 200.0e3,  // trough width = 200km; not the same
                                         //   as stream width
                       plateau = 2000.0;
  PetscScalar x_loc, y_loc;

  ierr = vbed.set(plateau); CHKERRQ(ierr);

  ierr = vbed.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar x = grid.x[i], y = grid.y[j];
      // note we treat exper P2 like others; it is flat anyway (slope=0)
      for (PetscInt m=0; m<4; m++) {
        PetscScalar drop = e[exper_chosen].bed_end_depth[m],
                    slope = drop / stream_length;
        if (inStream((pi/2.0)*m,width,x,y,x_loc,y_loc))
          vbed(i,j) = plateau - slope * x_loc * cos(pi * y_loc / width);
      }
    }
  }

  ierr = vbed.end_access(); CHKERRQ(ierr);

  // communicate b because it will be horizontally differentiated
  ierr = vbed.beginGhostComm(); CHKERRQ(ierr);
  ierr = vbed.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IcePSTexModel::additionalAtEndTimestep() {
  PetscErrorCode ierr;

  if (ivol == PETSC_NULL)  return 0;

  // all are in MKS units
  // variables starting with "g" are global across all processors

  // not all of these are used ...
  PetscScalar     gvolume, garea;
  ierr = volumeArea(gvolume, garea); CHKERRQ(ierr);

  PetscScalar     maxcbarALL = 0.0, 
                  areaup[4] = {0.0, 0.0, 0.0, 0.0},
                  avcup[4] = {0.0, 0.0, 0.0, 0.0},
                  areadown[4] = {0.0, 0.0, 0.0, 0.0},
                  avcdown[4] = {0.0, 0.0, 0.0, 0.0};
  PetscScalar     gmaxcbarALL, gareaup[4], gavcup[4], gareadown[4], gavcdown[4];
  PetscScalar     x_loc, y_loc;
  const PetscScalar darea = grid.dx * grid.dy;
  
  IceModelVec *tmp;
  ierr = diagnostics["velbar"]->compute(tmp); CHKERRQ(ierr);
  IceModelVec2V *vel_bar = dynamic_cast<IceModelVec2V*>(tmp);
  if (!vel_bar) SETERRQ(grid.com, 1, "dynamic_cast failure");

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vel_bar->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
        const PetscScalar cbar = (*vel_bar)(i,j).magnitude();
        const PetscScalar x = -grid.Ly + grid.dy * j, // note reversal (FIXME!
						      // do we need this now?)
                          y = -grid.Lx + grid.dx * i,
                          r = sqrt(PetscSqr(x) + PetscSqr(y));
        if (cbar > maxcbarALL)  maxcbarALL = cbar;
        if (exper_chosen == 3) { // exper P2
          const PetscScalar width = e[exper_chosen].stream_width[0] * 1000.0;
          for (int m=0; m<3; m++) {
            if (inStream((pi/180.0) * stream_angle_P2[m],width,x,y,x_loc,y_loc)) {
              if (r > stream_change) {
                areadown[m] += darea;
                avcdown[m] += cbar * darea;
              } else {
                areaup[m] += darea;
                avcup[m] += cbar * darea;
              }
            }
          }
        } else {
          for (int m=0; m<4; m++) {
            const PetscScalar width = e[exper_chosen].stream_width[m] * 1000.0;
            if (inStream((pi/2.0)*m,width,x,y,x_loc,y_loc)) {
              if (x_loc > stream_change - 1.0) {
                areadown[m] += darea;
                avcdown[m] += cbar * darea;
              } else {
                areaup[m] += darea;
                avcup[m] += cbar * darea;
              }
            }
          }
        }
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vel_bar->end_access(); CHKERRQ(ierr);

  delete tmp;

  // globalize and actually compute averages
  ierr = PISMGlobalMax(&maxcbarALL, &gmaxcbarALL, grid.com); CHKERRQ(ierr);
  for (PetscInt m=0; m<4; m++) {
    ierr = PISMGlobalSum(&areaup[m], &gareaup[m], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&avcup[m], &gavcup[m], grid.com); CHKERRQ(ierr);
    if (gareaup[m] > 0.0) {
      gavcup[m] = gavcup[m] / gareaup[m];
    } else {
      gavcup[m] = 0.0;
    }
    ierr = PISMGlobalSum(&areadown[m], &gareadown[m], grid.com); CHKERRQ(ierr);
    ierr = PISMGlobalSum(&avcdown[m], &gavcdown[m], grid.com); CHKERRQ(ierr);
    if (gareadown[m] > 0.0) {
      gavcdown[m] = gavcdown[m] / gareadown[m];
    } else {
      gavcdown[m] = 0.0;
    }
  }

  double dt_years = convert(dt, "seconds", "years"),
    a = grid.time->seconds_to_years(grid.time->current() - dt),
    b = grid.time->seconds_to_years(grid.time->current());

  dt_ser->append(dt_years, a, b);
  dt_ser->interp(a, b);

  ivol->append(gvolume, a, b);
  ivol->interp(a, b);
  iarea->append(garea, a, b);
  iarea->interp(a, b);
  maxcbar->append(gmaxcbarALL, a, b);
  maxcbar->interp(a, b);

  avup0->append(gavcup[0], a, b);
  avup0->interp(a, b);
  avup1->append(gavcup[1], a, b);
  avup1->interp(a, b);
  avup2->append(gavcup[2], a, b);
  avup2->interp(a, b);

  avdwn0->append(gavcdown[0], a, b);
  avdwn0->interp(a, b);
  avdwn1->append(gavcdown[1], a, b);
  avdwn1->interp(a, b);
  avdwn2->append(gavcdown[2], a, b);
  avdwn2->interp(a, b);

  if (exper_chosen != 3) {
    avup3->append(gavcup[3], a, b);
    avup3->interp(a, b);
    avdwn3->append(gavcdown[3], a, b);
    avdwn3->interp(a, b);
  }

  return 0;
}

/*! */
PetscErrorCode PSTYieldStress::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = PISMMohrCoulombYieldStress::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
                    "  setting phi = (till friction angle) for PST exper '%s' ...\n", experiment_name.c_str());
  CHKERRQ(ierr);

  ierr = init_till_phi(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSTYieldStress::init_till_phi() {
  PetscErrorCode ierr;

  const PetscScalar    dx = grid.dx, dy = grid.dy;
  PetscScalar          x_loc, y_loc;

  if (experiment <= 1)
    return 0;  // nothing further for P0A and P0I

  ierr = till_phi.begin_access(); CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      const PetscScalar x = -grid.Ly + dy * j, // note reversal
                        y = -grid.Lx + dx * i;
      if (experiment == 3) { // experiment P2
        const PetscScalar width = e[experiment].stream_width[0] * 1000.0;
        for (PetscInt m = 0; m < 3; m++) {
          if (inStreamNbhd(false, (pi / 180.0) * stream_angle_P2[m], width, x, y, x_loc, y_loc))
            till_phi(i, j) = phiLocal(width, x_loc, y_loc, DEFAULT_PHI_STRONG,
                                     e[experiment].upstream_phi[m],
                                     e[experiment].downstream_phi[m]);
        }
      } else {
        for (PetscInt m = 0; m < 4; m++) { // go through four sectors
          const PetscScalar width = e[experiment].stream_width[m] * 1000.0;
          if (inStreamNbhd(false, (pi / 2.0)*m, width, x, y, x_loc, y_loc)) {
            till_phi(i, j) = phiLocal(width, x_loc, y_loc, DEFAULT_PHI_STRONG,
                                     e[experiment].upstream_phi[m],
                                     e[experiment].downstream_phi[m]);
          }
        }
      }
    }
  }

  ierr = till_phi.end_access(); CHKERRQ(ierr);

  // communicate ghosts so that the tauc computation can be performed locally
  // (including ghosts of tauc, that is)
  ierr = till_phi.beginGhostComm(); CHKERRQ(ierr);
  ierr = till_phi.endGhostComm(); CHKERRQ(ierr);

  return 0;
}

int PSTYieldStress::sectorNumberP2(const PetscScalar x, const PetscScalar y) {
  if (x > 0.0) {
    if (y < x)
      return 0;
    else
      return 1;
  } else {
    if (y < 0.0)
      return 2;
    else 
      return 1;
  }
}


//! Compute the till friction angle in local strip coordinates.
PetscScalar PSTYieldStress::phiLocal(const PetscScalar width,
              const PetscScalar x, const PetscScalar y,
              const PetscScalar STRONG, 
              const PetscScalar UP, const PetscScalar DOWN) {

  const PetscScalar eta   = y / (width/2.0),   // normalized local y
                    xi    = x / stream_change; // normalized local x

  // compute lambda(eta) which is even and in [0,1]
  PetscScalar lambda = 0.0; // for big eta
  if (PetscAbs(eta) <= 1.0 - eta_slop)
    lambda = 1.0;
  else if (PetscAbs(eta) < 1.0 + eta_slop)
    lambda = 0.5 - 0.5 * sin((pi/2.0) * (PetscAbs(eta) - 1.0) / eta_slop);

  if (x > stream_change)
    return DOWN * lambda + STRONG * (1.0 - lambda); // downstream value
  else { // f(xi) is for upstream part only
    PetscScalar f = STRONG;
    if (xi >= xi_slop)
      f = UP;
    else if (xi > - xi_slop) {
      const PetscScalar fav = 0.5 * (STRONG + UP);
      f = fav - 0.5 * (STRONG - UP) * sin((pi/2.0) * (xi / xi_slop));
    }
    return f * lambda + STRONG * (1.0 - lambda); // upstream value
  }
}
