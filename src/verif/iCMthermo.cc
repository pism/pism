// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include <petscdmda.h>
#include "tests/exactTestsFG.h"
#include "tests/exactTestK.h"
#include "tests/exactTestO.h"
#include "iceCompModel.hh"
#include "PISMStressBalance.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "error_handling.hh"

namespace pism {

// boundary conditions for tests F, G (same as EISMINT II Experiment F)
const double IceCompModel::Ggeo = 0.042;
const double IceCompModel::ST = 1.67e-5;
const double IceCompModel::Tmin = 223.15;  // K
const double IceCompModel::LforFG = 750000; // m
const double IceCompModel::ApforG = 200; // m


/*! Re-implemented so that we can add compensatory strain_heating in Tests F and G. */
PetscErrorCode IceCompModel::temperatureStep(double* vertSacrCount, double* bulgeCount) {
  PetscErrorCode  ierr;

  if ((testname == 'F') || (testname == 'G')) {
    IceModelVec3 *strain_heating3;
    ierr = stress_balance->get_volumetric_strain_heating(strain_heating3); CHKERRQ(ierr);

    ierr = strain_heating3->add(1.0, strain_heating3_comp); CHKERRQ(ierr);      // strain_heating = strain_heating + strain_heating_c
    ierr = IceModel::temperatureStep(vertSacrCount,bulgeCount); CHKERRQ(ierr);
    ierr = strain_heating3->add(-1.0, strain_heating3_comp); CHKERRQ(ierr); // strain_heating = strain_heating - strain_heating_c
  } else {
    ierr = IceModel::temperatureStep(vertSacrCount,bulgeCount); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceCompModel::initTestFG() {
  PetscErrorCode  ierr;
  int        Mz = grid.Mz;
  double     H, accum;
  double     *dummy1, *dummy2, *dummy3, *dummy4;

  dummy1 = new double[Mz];  dummy2 = new double[Mz];
  dummy3 = new double[Mz];  dummy4 = new double[Mz];

  ierr = bed_topography.set(0); CHKERRQ(ierr);
  ierr = geothermal_flux.set(Ggeo); CHKERRQ(ierr);

  double *T = new double[grid.Mz];

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(T3);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = grid.radius(i, j);

    if (r > LforFG - 1.0) { // if (essentially) outside of sheet
      ice_thickness(i, j) = 0.0;
      for (int k = 0; k < Mz; k++) {
        T[k] = Tmin + ST * r;
      }
    } else {
      r = PetscMax(r, 1.0); // avoid singularity at origin
      if (testname == 'F') {
        bothexact(0.0, r, &grid.zlevels[0], Mz, 0.0,
                  &H, &accum, T, dummy1, dummy2, dummy3, dummy4);
        ice_thickness(i, j) = H;

      } else {
        bothexact(grid.time->current(), r, &grid.zlevels[0], Mz, ApforG,
                  &H, &accum, T, dummy1, dummy2, dummy3, dummy4);
        ice_thickness(i, j) = H;

      }
    }
    ierr = T3.setInternalColumn(i, j, T); CHKERRQ(ierr);
  }

  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);

  ierr = T3.update_ghosts(); CHKERRQ(ierr);

  ierr = ice_thickness.copy_to(ice_surface_elevation); CHKERRQ(ierr);

  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  delete [] T;

  return 0;
}


PetscErrorCode IceCompModel::getCompSourcesTestFG() {
  PetscErrorCode  ierr;
  double     accum, dummy0, *dummy1, *dummy2, *dummy3, *dummy4;

  dummy1=new double[grid.Mz];
  dummy2=new double[grid.Mz];
  dummy3=new double[grid.Mz];
  dummy4=new double[grid.Mz];

  double *strain_heating_C;
  strain_heating_C = new double[grid.Mz];

  const double
    ice_rho   = config.get("ice_density"),
    ice_c     = config.get("ice_specific_heat_capacity");

  // before temperature and flow step, set strain_heating_c from exact values

  IceModelVec::AccessList list;
  list.add(strain_heating3_comp);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = grid.radius(i, j);
    if (r > LforFG - 1.0) {  // outside of sheet
      ierr = strain_heating3_comp.setColumn(i, j, 0.0); CHKERRQ(ierr);
    } else {
      r = PetscMax(r, 1.0); // avoid singularity at origin
      if (testname == 'F') {
        bothexact(0.0, r, &grid.zlevels[0], grid.Mz, 0.0,
                  &dummy0, &accum, dummy1, dummy2, dummy3, dummy4, strain_heating_C);
      } else {
        bothexact(grid.time->current(), r, &grid.zlevels[0], grid.Mz, ApforG,
                  &dummy0, &accum, dummy1, dummy2, dummy3, dummy4, strain_heating_C);
      }
      for (unsigned int k=0;  k<grid.Mz;  k++) {
        // scale strain_heating to J/(s m^3)
        strain_heating_C[k] = strain_heating_C[k] * ice_rho * ice_c;
      }
      ierr = strain_heating3_comp.setInternalColumn(i, j, strain_heating_C); CHKERRQ(ierr);
    }
  }

  delete [] dummy1;
  delete [] dummy2;
  delete [] dummy3;
  delete [] dummy4;
  delete [] strain_heating_C;

  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestFG() {
  // fills Vecs ice_thickness, ice_surface_elevation, vAccum, T3, u3, v3, w3, strain_heating3, v_strain_heating_Comp
  PetscErrorCode  ierr;
  double     H, accum;
  double     Ts, *Uradial;

  IceModelVec3 *u3, *v3, *w3, *strain_heating3;
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);
  ierr = stress_balance->get_volumetric_strain_heating(strain_heating3); CHKERRQ(ierr);

  Uradial = new double[grid.Mz];

  double *T, *u, *v, *w, *strain_heating, *strain_heating_C;
  T = new double[grid.Mz];
  u = new double[grid.Mz];
  v = new double[grid.Mz];
  w = new double[grid.Mz];
  strain_heating = new double[grid.Mz];
  strain_heating_C = new double[grid.Mz];

  const double
    ice_rho   = config.get("ice_density"),
    ice_c     = config.get("ice_specific_heat_capacity");

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(T3);
  list.add(*u3);
  list.add(*v3);
  list.add(*w3);
  list.add(*strain_heating3);
  list.add(strain_heating3_comp);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double xx = grid.x[i], yy = grid.y[j], r = grid.radius(i, j);
    if (r > LforFG - 1.0) {  // outside of sheet

      ice_thickness(i, j) = 0.0;
      Ts = Tmin + ST * r;
      ierr = T3.setColumn(i, j, Ts); CHKERRQ(ierr);
      ierr = u3->setColumn(i, j, 0.0); CHKERRQ(ierr);
      ierr = v3->setColumn(i, j, 0.0); CHKERRQ(ierr);
      ierr = w3->setColumn(i, j, 0.0); CHKERRQ(ierr);
      ierr = strain_heating3->setColumn(i, j, 0.0); CHKERRQ(ierr);
      ierr = strain_heating3_comp.setColumn(i, j, 0.0); CHKERRQ(ierr);
    } else {  // inside the sheet
      r = PetscMax(r, 1.0); // avoid singularity at origin
      if (testname == 'F') {
        bothexact(0.0, r, &grid.zlevels[0], grid.Mz, 0.0,
                  &H, &accum, T, Uradial, w, strain_heating, strain_heating_C);
        ice_thickness(i,j)   = H;

      } else {
        bothexact(grid.time->current(), r, &grid.zlevels[0], grid.Mz, ApforG,
                  &H, &accum, T, Uradial, w, strain_heating, strain_heating_C);
        ice_thickness(i,j)   = H;

      }
      for (unsigned int k = 0; k < grid.Mz; k++) {
        u[k] = Uradial[k]*(xx/r);
        v[k] = Uradial[k]*(yy/r);
        strain_heating[k] = strain_heating[k] * ice_rho * ice_c; // scale strain_heating to J/(s m^3)
        strain_heating_C[k] = strain_heating_C[k] * ice_rho * ice_c; // scale strain_heating_C to J/(s m^3)
      }
      ierr = T3.setInternalColumn(i, j, T); CHKERRQ(ierr);
      ierr = u3->setInternalColumn(i, j, u); CHKERRQ(ierr);
      ierr = v3->setInternalColumn(i, j, v); CHKERRQ(ierr);
      ierr = w3->setInternalColumn(i, j, w); CHKERRQ(ierr);
      ierr = strain_heating3->setInternalColumn(i, j, strain_heating); CHKERRQ(ierr);
      ierr = strain_heating3_comp.setInternalColumn(i, j, strain_heating_C); CHKERRQ(ierr);
    }
  }

  delete [] Uradial;
  delete [] T;
  delete [] u;
  delete [] v;
  delete [] w;
  delete [] strain_heating;
  delete [] strain_heating_C;

  ierr = ice_thickness.update_ghosts(); CHKERRQ(ierr);
  ierr = ice_thickness.copy_to(ice_surface_elevation); CHKERRQ(ierr);

  ierr = T3.update_ghosts(); CHKERRQ(ierr);

  ierr = u3->update_ghosts(); CHKERRQ(ierr);

  ierr = v3->update_ghosts(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceCompModel::computeTemperatureErrors(double &gmaxTerr,
                                                      double &gavTerr) {

  PetscErrorCode ierr;
  double    maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;

  double   *dummy1, *dummy2, *dummy3, *dummy4, *Tex;
  double   junk0, junk1;

  Tex = new double[grid.Mz];
  dummy1 = new double[grid.Mz];  dummy2 = new double[grid.Mz];
  dummy3 = new double[grid.Mz];  dummy4 = new double[grid.Mz];

  double *T;

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(T3);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = grid.radius(i, j);
    ierr = T3.getInternalColumn(i, j, &T); CHKERRQ(ierr);
    if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet
      // and not at central singularity
      switch (testname) {
      case 'F':
        bothexact(0.0, r, &grid.zlevels[0], grid.Mz, 0.0,
                  &junk0, &junk1, Tex, dummy1, dummy2, dummy3, dummy4);
        break;
      case 'G':
        bothexact(grid.time->current(), r, &grid.zlevels[0], grid.Mz, ApforG,
                  &junk0, &junk1, Tex, dummy1, dummy2, dummy3, dummy4);
        break;
      default:
        throw RuntimeError("temperature errors only computable for tests F and G");
      }
      const int ks = grid.kBelowHeight(ice_thickness(i,j));
      for (int k = 0; k < ks; k++) {  // only eval error if below num surface
        const double Terr = PetscAbs(T[k] - Tex[k]);
        maxTerr = PetscMax(maxTerr, Terr);
        avcount += 1.0;
        avTerr += Terr;
      }
    }
  }

  delete [] Tex;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;

  ierr = GlobalMax(grid.com, &maxTerr,  &gmaxTerr); CHKERRQ(ierr);
  ierr = GlobalSum(grid.com, &avTerr,  &gavTerr); CHKERRQ(ierr);
  double  gavcount;
  ierr = GlobalSum(grid.com, &avcount,  &gavcount); CHKERRQ(ierr);
  gavTerr = gavTerr/PetscMax(gavcount, 1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeIceBedrockTemperatureErrors(
                                double &gmaxTerr, double &gavTerr,
                                double &gmaxTberr, double &gavTberr) {
  PetscErrorCode ierr;

  if ((testname != 'K') && (testname != 'O')) {
    throw RuntimeError("ice and bedrock temperature errors only computable for tests K and O");
  }

  double    maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;
  double    maxTberr = 0.0, avTberr = 0.0, avbcount = 0.0;

  double    *Tex, *Tbex, *T, *Tb;
  double    FF;
  Tex = new double[grid.Mz];

  IceModelVec3BTU *bedrock_temp;

  BTU_Verification *my_btu = dynamic_cast<BTU_Verification*>(btu);
  if (my_btu == NULL) {
    throw RuntimeError("my_btu == NULL");
  }
  ierr = my_btu->get_temp(bedrock_temp); CHKERRQ(ierr);

  std::vector<double> zblevels = bedrock_temp->get_levels();
  unsigned int Mbz = (unsigned int)zblevels.size();
  Tbex = new double[Mbz];

  switch (testname) {
    case 'K':
      for (unsigned int k = 0; k < grid.Mz; k++) {
        ierr = exactK(grid.time->current(), grid.zlevels[k], &Tex[k], &FF,
                      (bedrock_is_ice_forK==PETSC_TRUE)); CHKERRQ(ierr);
      }
      for (unsigned int k = 0; k < Mbz; k++) {
        ierr = exactK(grid.time->current(), zblevels[k], &Tbex[k], &FF,
                      (bedrock_is_ice_forK==PETSC_TRUE)); CHKERRQ(ierr);
      }
      break;
    case 'O':
      double dum1, dum2, dum3, dum4;
      for (unsigned int k = 0; k < grid.Mz; k++) {
        ierr = exactO(grid.zlevels[k], &Tex[k], &dum1, &dum2, &dum3, &dum4);
             CHKERRQ(ierr);
      }
      for (unsigned int k = 0; k < Mbz; k++) {
        ierr = exactO(zblevels[k], &Tbex[k], &dum1, &dum2, &dum3, &dum4);
             CHKERRQ(ierr);
      }
      break;
    default:
      throw RuntimeError("ice and bedrock temperature errors only for tests K and O");
  }

  IceModelVec::AccessList list;
  list.add(T3);
  list.add(*bedrock_temp);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    ierr = bedrock_temp->getInternalColumn(i,j,&Tb); CHKERRQ(ierr);
    for (unsigned int kb = 0; kb < Mbz; kb++) {
      const double Tberr = PetscAbs(Tb[kb] - Tbex[kb]);
      maxTberr = PetscMax(maxTberr,Tberr);
      avbcount += 1.0;
      avTberr += Tberr;
    }
    ierr = T3.getInternalColumn(i,j,&T); CHKERRQ(ierr);
    for (unsigned int k = 0; k < grid.Mz; k++) {
      const double Terr = PetscAbs(T[k] - Tex[k]);
      maxTerr = PetscMax(maxTerr,Terr);
      avcount += 1.0;
      avTerr += Terr;
    }
  }

  delete [] Tex;  delete [] Tbex;

  ierr = GlobalMax(grid.com, &maxTerr,  &gmaxTerr); CHKERRQ(ierr);
  ierr = GlobalSum(grid.com, &avTerr,  &gavTerr); CHKERRQ(ierr);
  double  gavcount;
  ierr = GlobalSum(grid.com, &avcount,  &gavcount); CHKERRQ(ierr);
  gavTerr = gavTerr/PetscMax(gavcount,1.0);  // avoid div by zero

  ierr = GlobalMax(grid.com, &maxTberr,  &gmaxTberr); CHKERRQ(ierr);
  ierr = GlobalSum(grid.com, &avTberr,  &gavTberr); CHKERRQ(ierr);
  double  gavbcount;
  ierr = GlobalSum(grid.com, &avbcount,  &gavbcount); CHKERRQ(ierr);
  gavTberr = gavTberr/PetscMax(gavbcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeBasalTemperatureErrors(
      double &gmaxTerr, double &gavTerr, double &centerTerr) {

  PetscErrorCode  ierr;
  double     domeT, domeTexact, Terr, avTerr;

  double     dummy, z, Texact, dummy1, dummy2, dummy3, dummy4, dummy5;

  IceModelVec::AccessList list(T3);

  domeT=0; domeTexact = 0; Terr=0; avTerr=0;

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = grid.radius(i,j);
    switch (testname) {
    case 'F':
      if (r > LforFG - 1.0) {  // outside of sheet
        Texact=Tmin + ST * r;  // = Ts
      } else {
        r=PetscMax(r,1.0);
        z=0.0;
        bothexact(0.0,r,&z,1,0.0,
                  &dummy5,&dummy,&Texact,&dummy1,&dummy2,&dummy3,&dummy4);
      }
      break;
    case 'G':
      if (r > LforFG -1.0) {  // outside of sheet
        Texact=Tmin + ST * r;  // = Ts
      } else {
        r=PetscMax(r,1.0);
        z=0.0;
        bothexact(grid.time->current(),r,&z,1,ApforG,
                  &dummy5,&dummy,&Texact,&dummy1,&dummy2,&dummy3,&dummy4);
      }
      break;
    default:
      throw RuntimeError("temperature errors only computable for tests F and G");
    }

    const double Tbase = T3.getValZ(i,j,0.0);
    if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
      domeT = Tbase;
      domeTexact = Texact;
    }
    // compute maximum errors
    Terr = PetscMax(Terr,PetscAbsReal(Tbase - Texact));
    // add to sums for average errors
    avTerr += PetscAbs(Tbase - Texact);
  }

  double gdomeT, gdomeTexact;

  ierr = GlobalMax(grid.com, &Terr,  &gmaxTerr); CHKERRQ(ierr);
  ierr = GlobalSum(grid.com, &avTerr,  &gavTerr); CHKERRQ(ierr);
  gavTerr = gavTerr/(grid.Mx*grid.My);
  ierr = GlobalMax(grid.com, &domeT,  &gdomeT); CHKERRQ(ierr);
  ierr = GlobalMax(grid.com, &domeTexact,  &gdomeTexact); CHKERRQ(ierr);
  centerTerr = PetscAbsReal(gdomeT - gdomeTexact);

  return 0;
}


PetscErrorCode IceCompModel::compute_strain_heating_errors(
                      double &gmax_strain_heating_err, double &gav_strain_heating_err) {

  PetscErrorCode ierr;
  double    max_strain_heating_err = 0.0, av_strain_heating_err = 0.0, avcount = 0.0;

  double   *dummy1, *dummy2, *dummy3, *dummy4, *strain_heating_exact;
  double   junk0, junk1;

  strain_heating_exact = new double[grid.Mz];
  dummy1 = new double[grid.Mz];  dummy2 = new double[grid.Mz];
  dummy3 = new double[grid.Mz];  dummy4 = new double[grid.Mz];

  const double
    ice_rho   = config.get("ice_density"),
    ice_c     = config.get("ice_specific_heat_capacity");

  double *strain_heating;
  IceModelVec3 *strain_heating3;
  ierr = stress_balance->get_volumetric_strain_heating(strain_heating3); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(*strain_heating3);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = grid.radius(i,j);
    if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet
      // and not at central singularity
      switch (testname) {
      case 'F':
        bothexact(0.0,r,&grid.zlevels[0],grid.Mz,0.0,
                  &junk0,&junk1,dummy1,dummy2,dummy3,strain_heating_exact,dummy4);
        break;
      case 'G':
        bothexact(grid.time->current(),r,&grid.zlevels[0],grid.Mz,ApforG,
                  &junk0,&junk1,dummy1,dummy2,dummy3,strain_heating_exact,dummy4);
        break;
      default:
        throw RuntimeError("strain-heating (strain_heating) errors only computable for tests F and G");
      }
      for (unsigned int k = 0; k < grid.Mz; k++) {
        // scale exact strain_heating to J/(s m^3)
        strain_heating_exact[k] *= ice_rho * ice_c;
      }
      const unsigned int ks = grid.kBelowHeight(ice_thickness(i,j));
      ierr = strain_heating3->getInternalColumn(i,j,&strain_heating); CHKERRQ(ierr);
      for (unsigned int k = 0; k < ks; k++) {  // only eval error if below num surface
        const double strain_heating_err = PetscAbs(strain_heating[k] - strain_heating_exact[k]);
        max_strain_heating_err = PetscMax(max_strain_heating_err,strain_heating_err);
        avcount += 1.0;
        av_strain_heating_err += strain_heating_err;
      }
    }
  }

  delete [] strain_heating_exact;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;

  ierr = GlobalMax(grid.com, &max_strain_heating_err,  &gmax_strain_heating_err); CHKERRQ(ierr);
  ierr = GlobalSum(grid.com, &av_strain_heating_err,  &gav_strain_heating_err); CHKERRQ(ierr);
  double  gavcount;
  ierr = GlobalSum(grid.com, &avcount,  &gavcount); CHKERRQ(ierr);
  gav_strain_heating_err = gav_strain_heating_err/PetscMax(gavcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeSurfaceVelocityErrors(double &gmaxUerr, double &gavUerr,
                                                          double &gmaxWerr, double &gavWerr) {

  PetscErrorCode ierr;
  double    maxUerr = 0.0, maxWerr = 0.0, avUerr = 0.0, avWerr = 0.0;

  IceModelVec3 *u3, *v3, *w3;
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(*u3);
  list.add(*v3);
  list.add(*w3);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double xx = grid.x[i], yy = grid.y[j], r = grid.radius(i, j);
    if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet
      // and not at central singularity
      double radialUex, wex;
      double dummy0, dummy1, dummy2, dummy3, dummy4;
      switch (testname) {
      case 'F':
        bothexact(0.0, r, &ice_thickness(i,j), 1, 0.0,
                  &dummy0, &dummy1, &dummy2, &radialUex, &wex, &dummy3, &dummy4);
        break;
      case 'G':
        bothexact(grid.time->current(), r, &ice_thickness(i,j), 1, ApforG,
                  &dummy0, &dummy1, &dummy2, &radialUex, &wex, &dummy3, &dummy4);
        break;
      default:
        throw RuntimeError("surface velocity errors only computed for tests F and G");
      }
      const double uex = (xx/r) * radialUex;
      const double vex = (yy/r) * radialUex;
      // note that because getValZ does linear interpolation and H[i][j] is not exactly at
      // a grid point, this causes nonzero errors even with option -eo
      const double Uerr = sqrt(PetscSqr(u3->getValZ(i, j, ice_thickness(i,j)) - uex)
                               + PetscSqr(v3->getValZ(i, j, ice_thickness(i,j)) - vex));
      maxUerr = PetscMax(maxUerr, Uerr);
      avUerr += Uerr;
      const double Werr = PetscAbs(w3->getValZ(i, j, ice_thickness(i,j)) - wex);
      maxWerr = PetscMax(maxWerr, Werr);
      avWerr += Werr;
    }
  }

  ierr = GlobalMax(grid.com, &maxUerr,  &gmaxUerr); CHKERRQ(ierr);
  ierr = GlobalMax(grid.com, &maxWerr,  &gmaxWerr); CHKERRQ(ierr);
  ierr = GlobalSum(grid.com, &avUerr,  &gavUerr); CHKERRQ(ierr);
  gavUerr = gavUerr/(grid.Mx*grid.My);
  ierr = GlobalSum(grid.com, &avWerr,  &gavWerr); CHKERRQ(ierr);
  gavWerr = gavWerr/(grid.Mx*grid.My);
  return 0;
}


PetscErrorCode IceCompModel::computeBasalMeltRateErrors(
                   double &gmaxbmelterr, double &gminbmelterr) {
  PetscErrorCode ierr;
  double    maxbmelterr = -9.99e40, minbmelterr = 9.99e40, err;
  double    bmelt,dum1,dum2,dum3,dum4;

  if (testname != 'O') {
    throw RuntimeError("basal melt rate errors are only computable for test O");
  }

  // we just need one constant from exact solution:
  ierr = exactO(0.0, &dum1, &dum2, &dum3, &dum4, &bmelt); CHKERRQ(ierr);

  IceModelVec::AccessList list(basal_melt_rate);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    err = PetscAbs(basal_melt_rate(i,j) - bmelt);
    maxbmelterr = PetscMax(maxbmelterr, err);
    minbmelterr = PetscMin(minbmelterr, err);
  }

  ierr = GlobalMax(grid.com, &maxbmelterr,  &gmaxbmelterr); CHKERRQ(ierr);
  ierr = GlobalMin(grid.com, &minbmelterr,  &gminbmelterr); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::fillTemperatureSolnTestsKO() {
  PetscErrorCode    ierr;

  double       *Tcol;
  double       dum1, dum2, dum3, dum4;
  double    FF;
  Tcol = new double[grid.Mz];

  // evaluate exact solution in a column; all columns are the same
  switch (testname) {
    case 'K':
      for (unsigned int k=0; k<grid.Mz; k++) {
        ierr = exactK(grid.time->current(), grid.zlevels[k], &Tcol[k], &FF,
                      (bedrock_is_ice_forK==PETSC_TRUE)); CHKERRQ(ierr);
      }
      break;
    case 'O':
      for (unsigned int k=0; k<grid.Mz; k++) {
        ierr = exactO(grid.zlevels[k], &Tcol[k], &dum1, &dum2, &dum3, &dum4); CHKERRQ(ierr);
      }
      break;
    default:
      throw RuntimeError("only fills temperature solutions for tests K and O");
  }

  // copy column values into 3D arrays
  IceModelVec::AccessList list(T3);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    ierr = T3.setInternalColumn(i,j,Tcol); CHKERRQ(ierr);
  }

  delete [] Tcol;

  // communicate T
  ierr = T3.update_ghosts(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::fillBasalMeltRateSolnTestO() {
  PetscErrorCode    ierr;
  double       bmelt, dum1, dum2, dum3, dum4;
  if (testname != 'O') {
    throw RuntimeError("only fills basal melt rate soln for test O");
  }

  // we just need one constant from exact solution:
  ierr = exactO(0.0, &dum1, &dum2, &dum3, &dum4, &bmelt); CHKERRQ(ierr);

  ierr = basal_melt_rate.set(bmelt); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::initTestsKO() {
  PetscErrorCode    ierr;

  if (testname == 'K') {
    bool Mbz_set;
    int Mbz;
    ierr = OptionsInt("-Mbz", "Number of levels in the bedrock thermal model",
                          Mbz, Mbz_set); CHKERRQ(ierr);
    if (Mbz_set && Mbz < 2) {
      throw RuntimeError("pismv test K requires a bedrock thermal layer 1000m deep");
    }
  }

  ierr = bed_topography.set(0.0); CHKERRQ(ierr);
  ierr = geothermal_flux.set(0.042); CHKERRQ(ierr);
  ierr = ice_thickness.set(3000.0); CHKERRQ(ierr);
  ierr = ice_thickness.copy_to(ice_surface_elevation); CHKERRQ(ierr);

  ierr = fillTemperatureSolnTestsKO(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode BTU_Verification::get_temp(IceModelVec3BTU* &result) {
  result = &temp;
  return 0;
}

PetscErrorCode BTU_Verification::bootstrap() {
  PetscErrorCode ierr;

  if (Mbz < 2) {
    return 0;
  }

  std::vector<double> Tbcol(Mbz),
    zlevels = temp.get_levels();
  double dum1, dum2, dum3, dum4;
  double FF;

  // evaluate exact solution in a column; all columns are the same
  switch (testname) {
    case 'K':
      for (unsigned int k=0; k<Mbz; k++) {
        if (exactK(grid.time->current(), zlevels[k], &Tbcol[k], &FF,
                   (bedrock_is_ice==PETSC_TRUE))) {
          throw RuntimeError::formatted("exactK() reports that level %9.7f is below B0 = -1000.0 m",
                                        zlevels[k]);
        }
      }
      break;
    case 'O':
      for (unsigned int k=0; k<Mbz; k++) {
        ierr = exactO(zlevels[k], &Tbcol[k], &dum1, &dum2, &dum3, &dum4); CHKERRQ(ierr);
      }
      break;
    default:
      {
        ierr = BedThermalUnit::bootstrap(); CHKERRQ(ierr);
      }
  }

  // copy column values into 3D arrays
  IceModelVec::AccessList list(temp);

  for (Points p(grid); p; p.next()) {
    ierr = temp.setInternalColumn(p.i(), p.j(), &Tbcol[0]); CHKERRQ(ierr);
  }

  return 0;
}

} // end of namespace pism
