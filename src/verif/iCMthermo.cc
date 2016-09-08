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

#include <cmath>

#include "tests/exactTestsFG.hh"
#include "tests/exactTestK.h"
#include "tests/exactTestO.h"
#include "iceCompModel.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/PISMTime.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_options.hh"
#include "base/util/error_handling.hh"
#include "earth/PISMBedDef.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/pism_utilities.hh"
#include "BTU_Verification.hh"

namespace pism {

// boundary conditions for tests F, G (same as EISMINT II Experiment F)
const double IceCompModel::ST = 1.67e-5;
const double IceCompModel::Tmin = 223.15;  // K
const double IceCompModel::LforFG = 750000; // m
const double IceCompModel::ApforG = 200; // m


/*! Re-implemented so that we can add compensatory strain_heating in Tests F and G. */
void IceCompModel::temperatureStep(unsigned int *vertSacrCount, unsigned int *bulgeCount) {

  if ((testname == 'F') || (testname == 'G')) {
    // FIXME: This code messes with the strain heating field owned by
    // stress_balance. This is BAD.
    IceModelVec3 &strain_heating3 = const_cast<IceModelVec3&>(m_stress_balance->volumetric_strain_heating());

    strain_heating3.add(1.0, strain_heating3_comp);      // strain_heating = strain_heating + strain_heating_c
    IceModel::temperatureStep(vertSacrCount, bulgeCount);
    strain_heating3.add(-1.0, strain_heating3_comp); // strain_heating = strain_heating - strain_heating_c
  } else {
    IceModel::temperatureStep(vertSacrCount, bulgeCount);
  }
}


void IceCompModel::initTestFG() {

  IceModelVec2S bed_topography;
  bed_topography.create(m_grid, "topg", WITHOUT_GHOSTS);
  bed_topography.set(0.0);
  m_beddef->set_elevation(bed_topography);

  IceModelVec::AccessList list;
  list.add(m_ice_thickness);
  list.add(m_ice_temperature);

  const double time = testname == 'F' ? 0.0 : m_time->current();
  const double A    = testname == 'F' ? 0.0 : ApforG;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // avoid singularity at origin
    const double r = std::max(radius(*m_grid, i, j), 1.0);

    if (r > LforFG - 1.0) { // if (essentially) outside of sheet
      m_ice_thickness(i, j) = 0.0;
      m_ice_temperature.set_column(i, j, Tmin + ST * r);
    } else {
      TestFGParameters P = exactFG(time, r, m_grid->z(), A);
      m_ice_thickness(i, j) = P.H;
      m_ice_temperature.set_column(i, j, &P.T[0]);
    }
  }

  m_ice_thickness.update_ghosts();

  m_ice_temperature.update_ghosts();

  m_ice_surface_elevation.copy_from(m_ice_thickness);
}


void IceCompModel::getCompSourcesTestFG() {

  const double
    ice_rho   = m_config->get_double("constants.ice.density"),
    ice_c     = m_config->get_double("constants.ice.specific_heat_capacity");

  // before temperature and flow step, set strain_heating_c from exact values

  IceModelVec::AccessList list;
  list.add(strain_heating3_comp);

  const double time = testname == 'F' ? 0.0 : m_time->current();
  const double A    = testname == 'F' ? 0.0 : ApforG;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = std::max(radius(*m_grid, i, j), 1.0); // avoid singularity at origin

    if (r > LforFG - 1.0) {  // outside of sheet
      strain_heating3_comp.set_column(i, j, 0.0);
    } else {
      TestFGParameters P = exactFG(time, r, m_grid->z(), A);

      strain_heating3_comp.set_column(i, j, &P.Sigc[0]);
    }
  }

  // scale strain_heating to J/(s m^3)
  strain_heating3_comp.scale(ice_rho * ice_c);
}


void IceCompModel::fillSolnTestFG() {
  // fills Vecs ice_thickness, ice_surface_elevation, vAccum, m_ice_temperature, u3, v3, w3, strain_heating3, v_strain_heating_Comp

  // FIXME: This code messes with the fields owned by stress_balance.
  // This is BAD.
  IceModelVec3
    &strain_heating3 = const_cast<IceModelVec3&>(m_stress_balance->volumetric_strain_heating()),
    &u3 = const_cast<IceModelVec3&>(m_stress_balance->velocity_u()),
    &v3 = const_cast<IceModelVec3&>(m_stress_balance->velocity_v()),
    &w3 = const_cast<IceModelVec3&>(m_stress_balance->velocity_w());

  std::vector<double> u(m_grid->Mz());
  std::vector<double> v(m_grid->Mz());

  const double
    ice_rho   = m_config->get_double("constants.ice.density"),
    ice_c     = m_config->get_double("constants.ice.specific_heat_capacity");

  const double time = testname == 'F' ? 0.0 : m_time->current();
  const double A    = testname == 'F' ? 0.0 : ApforG;

  IceModelVec::AccessList list;
  list.add(m_ice_thickness);
  list.add(m_ice_temperature);
  list.add(u3);
  list.add(v3);
  list.add(w3);
  list.add(strain_heating3);
  list.add(strain_heating3_comp);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      x = m_grid->x(i),
      y = m_grid->y(j),
      r = std::max(radius(*m_grid, i, j), 1.0); // avoid singularity at origin

    if (r > LforFG - 1.0) {  // outside of sheet
      m_ice_thickness(i, j) = 0.0;

      m_ice_temperature.set_column(i, j, Tmin + ST * r);

      strain_heating3.set_column(i, j, 0.0);
      strain_heating3_comp.set_column(i, j, 0.0);

      u3.set_column(i, j, 0.0);
      v3.set_column(i, j, 0.0);
      w3.set_column(i, j, 0.0);
    } else {  // inside the sheet
      TestFGParameters P = exactFG(time, r, m_grid->z(), A);

      m_ice_thickness(i, j) = P.H;

      m_ice_temperature.set_column(i, j, &P.T[0]);

      strain_heating3.set_column(i, j, &P.Sig[0]);
      strain_heating3_comp.set_column(i, j, &P.Sigc[0]);

      for (unsigned int k = 0; k < m_grid->Mz(); k++) {
        u[k] = P.U[k]*(x/r);
        v[k] = P.U[k]*(y/r);
      }

      u3.set_column(i, j, &u[0]);
      v3.set_column(i, j, &v[0]);
      w3.set_column(i, j, &P.w[0]);
    }
  }

  // scale strain_heating to J/(s m^3)
  strain_heating3.scale(ice_rho * ice_c);
  strain_heating3_comp.scale(ice_rho * ice_c);

  m_ice_thickness.update_ghosts();
  m_ice_surface_elevation.copy_from(m_ice_thickness);

  m_ice_temperature.update_ghosts();

  u3.update_ghosts();

  v3.update_ghosts();
}

void IceCompModel::computeTemperatureErrors(double &gmaxTerr,
                                            double &gavTerr) {
  double maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;

  if (testname != 'F' and testname != 'G') {
    throw RuntimeError("temperature errors only computable for tests F and G");
  }

  const double time = testname == 'F' ? 0.0 : m_time->current();
  const double A    = testname == 'F' ? 0.0 : ApforG;

  IceModelVec::AccessList list;
  list.add(m_ice_thickness);
  list.add(m_ice_temperature);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double r = radius(*m_grid, i, j);
      const double *T = m_ice_temperature.get_column(i, j);

      // only evaluate error if inside sheet and not at central
      // singularity
      if ((r >= 1.0) and (r <= LforFG - 1.0)) {
        TestFGParameters P = exactFG(time, r, m_grid->z(), A);

        // only evaluate error if below ice surface
        const int ks = m_grid->kBelowHeight(m_ice_thickness(i, j));
        for (int k = 0; k < ks; k++) {
          const double Terr = fabs(T[k] - P.T[k]);
          maxTerr = std::max(maxTerr, Terr);
          avcount += 1.0;
          avTerr += Terr;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  gmaxTerr = GlobalMax(m_grid->com, maxTerr);
  gavTerr = GlobalSum(m_grid->com, avTerr);
  double gavcount = GlobalSum(m_grid->com, avcount);
  gavTerr = gavTerr / std::max(gavcount, 1.0);  // avoid division by zero
}


void IceCompModel::computeIceBedrockTemperatureErrors(double &gmaxTerr, double &gavTerr,
                                                      double &gmaxTberr, double &gavTberr) {

  if (testname != 'K' and testname != 'O') {
    throw RuntimeError("ice and bedrock temperature errors only computable for tests K and O");
  }

  double    maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;
  double    maxTberr = 0.0, avTberr = 0.0, avbcount = 0.0;

  const double *Tb, *T;
  std::vector<double> Tex(m_grid->Mz());

  energy::BTU_Verification *my_btu = dynamic_cast<energy::BTU_Verification*>(m_btu);
  if (my_btu == NULL) {
    throw RuntimeError("BTU_Verification is required");
  }
  const IceModelVec3Custom &bedrock_temp = my_btu->temperature();

  std::vector<double> zblevels = bedrock_temp.get_levels();
  unsigned int Mbz = (unsigned int)zblevels.size();
  std::vector<double> Tbex(Mbz);

  switch (testname) {
    case 'K':
      for (unsigned int k = 0; k < m_grid->Mz(); k++) {
        TestKParameters K = exactK(m_time->current(), m_grid->z(k), bedrock_is_ice_forK == true);
        Tex[k] = K.T;
      }
      for (unsigned int k = 0; k < Mbz; k++) {
        TestKParameters K = exactK(m_time->current(), zblevels[k], bedrock_is_ice_forK == true);
        Tbex[k] = K.T;
      }
      break;
    case 'O':
      double dum1, dum2, dum3, dum4;
      for (unsigned int k = 0; k < m_grid->Mz(); k++) {
        exactO(m_grid->z(k), &Tex[k], &dum1, &dum2, &dum3, &dum4);
      }
      for (unsigned int k = 0; k < Mbz; k++) {
        exactO(zblevels[k], &Tbex[k], &dum1, &dum2, &dum3, &dum4);
      }
      break;
    default:
      throw RuntimeError("ice and bedrock temperature errors only for tests K and O");
  }

  IceModelVec::AccessList list;
  list.add(m_ice_temperature);
  list.add(bedrock_temp);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    Tb = bedrock_temp.get_column(i, j);
    for (unsigned int kb = 0; kb < Mbz; kb++) {
      const double Tberr = fabs(Tb[kb] - Tbex[kb]);
      maxTberr = std::max(maxTberr, Tberr);
      avbcount += 1.0;
      avTberr += Tberr;
    }
    T = m_ice_temperature.get_column(i, j);
    for (unsigned int k = 0; k < m_grid->Mz(); k++) {
      const double Terr = fabs(T[k] - Tex[k]);
      maxTerr = std::max(maxTerr, Terr);
      avcount += 1.0;
      avTerr += Terr;
    }
  }

  gmaxTerr = GlobalMax(m_grid->com, maxTerr);
  gavTerr = GlobalSum(m_grid->com, avTerr);
  double gavcount = GlobalSum(m_grid->com, avcount);
  gavTerr = gavTerr/std::max(gavcount, 1.0);  // avoid division by zero

  gmaxTberr = GlobalMax(m_grid->com, maxTberr);
  gavTberr = GlobalSum(m_grid->com, avTberr);
  double gavbcount = GlobalSum(m_grid->com, avbcount);
  gavTberr = gavTberr/std::max(gavbcount, 1.0);  // avoid division by zero
}


void IceCompModel::computeBasalTemperatureErrors(double &gmaxTerr, double &gavTerr, double &centerTerr) {

  if (testname != 'F' and testname != 'G') {
    throw RuntimeError("temperature errors only computable for tests F and G");
  }

  double
    Texact     = 0.0,
    domeT      = 0.0,
    domeTexact = 0.0,
    Terr       = 0.0,
    avTerr     = 0.0;

  const double time = testname == 'F' ? 0.0 : m_time->current();
  const double A    = testname == 'F' ? 0.0 : ApforG;
  std::vector<double> z(1, 0.0);

  IceModelVec::AccessList list(m_ice_temperature);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double r = std::max(radius(*m_grid, i, j), 1.0);

      if (r > LforFG - 1.0) { // outside of sheet
        Texact = Tmin + ST * r; // = Ts
      } else {
        Texact = exactFG(time, r, z, A).T[0];
      }

      const double Tbase = m_ice_temperature.get_column(i,j)[0];
      if (i == ((int)m_grid->Mx() - 1) / 2 and
          j == ((int)m_grid->My() - 1) / 2) {
        domeT      = Tbase;
        domeTexact = Texact;
      }
      // compute maximum errors
      Terr = std::max(Terr, fabs(Tbase - Texact));
      // add to sums for average errors
      avTerr += fabs(Tbase - Texact);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  double gdomeT, gdomeTexact;

  gmaxTerr    = GlobalMax(m_grid->com, Terr);
  gavTerr     = GlobalSum(m_grid->com, avTerr);
  gavTerr     = gavTerr/(m_grid->Mx()*m_grid->My());
  gdomeT      = GlobalMax(m_grid->com, domeT);
  gdomeTexact = GlobalMax(m_grid->com, domeTexact);
  centerTerr  = fabs(gdomeT - gdomeTexact);
}


void IceCompModel::compute_strain_heating_errors(double &gmax_strain_heating_err, double &gav_strain_heating_err) {
  double max_strain_heating_err = 0.0, av_strain_heating_err = 0.0, avcount = 0.0;

  if (testname != 'F' and testname != 'G') {
    throw RuntimeError("strain-heating (strain_heating) errors only computable for tests F and G");
  }

  const double
    ice_rho   = m_config->get_double("constants.ice.density"),
    ice_c     = m_config->get_double("constants.ice.specific_heat_capacity");

  const IceModelVec3 &strain_heating3 = m_stress_balance->volumetric_strain_heating();

  IceModelVec::AccessList list;
  list.add(m_ice_thickness);
  list.add(strain_heating3);

  const double time = testname == 'F' ? 0.0 : m_time->current();
  const double A    = testname == 'F' ? 0.0 : ApforG;

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double r = radius(*m_grid, i, j);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {
        // only evaluate error if inside sheet and not at central singularity

        TestFGParameters P = exactFG(time, r, m_grid->z(), A);

        for (unsigned int k = 0; k < m_grid->Mz(); k++) {
          // scale exact strain_heating to J/(s m^3)
          P.Sig[k] *= ice_rho * ice_c;
        }

        const unsigned int ks = m_grid->kBelowHeight(m_ice_thickness(i, j));
        const double *strain_heating = strain_heating3.get_column(i, j);

        for (unsigned int k = 0; k < ks; k++) {  // only evaluate error if below ice surface
          const double strain_heating_err = fabs(strain_heating[k] - P.Sig[k]);
          max_strain_heating_err = std::max(max_strain_heating_err, strain_heating_err);
          avcount += 1.0;
          av_strain_heating_err += strain_heating_err;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  gmax_strain_heating_err = GlobalMax(m_grid->com, max_strain_heating_err);
  gav_strain_heating_err = GlobalSum(m_grid->com, av_strain_heating_err);
  double gavcount = GlobalSum(m_grid->com, avcount);
  gav_strain_heating_err = gav_strain_heating_err/std::max(gavcount, 1.0);  // avoid div by zero
}


void IceCompModel::computeSurfaceVelocityErrors(double &gmaxUerr, double &gavUerr,
                                                double &gmaxWerr, double &gavWerr) {
  double
    maxUerr = 0.0,
    maxWerr = 0.0,
    avUerr  = 0.0,
    avWerr  = 0.0;

  const IceModelVec3
    &u3 = m_stress_balance->velocity_u(),
    &v3 = m_stress_balance->velocity_v(),
    &w3 = m_stress_balance->velocity_w();

  IceModelVec::AccessList list;
  list.add(m_ice_thickness);
  list.add(u3);
  list.add(v3);
  list.add(w3);

  const double time = testname == 'F' ? 0.0 : m_time->current();
  const double A    = testname == 'F' ? 0.0 : ApforG;

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double H = m_ice_thickness(i, j);
      std::vector<double> z(1, H);

      const double
        x = m_grid->x(i),
        y = m_grid->y(j),
        r = radius(*m_grid, i, j);

      if ((r >= 1.0) and (r <= LforFG - 1.0)) {
        // only evaluate error if inside sheet and not at central singularity

        TestFGParameters P = exactFG(time, r, z, A);

        const double
          uex = (x/r) * P.U[0],
          vex = (y/r) * P.U[0];
        // note that because getValZ does linear interpolation and H(i, j) is not exactly at
        // a grid point, this causes nonzero errors even with option -eo
        const double Uerr = sqrt(PetscSqr(u3.getValZ(i, j, H) - uex) +
                                 PetscSqr(v3.getValZ(i, j, H) - vex));
        maxUerr = std::max(maxUerr, Uerr);
        avUerr += Uerr;
        const double Werr = fabs(w3.getValZ(i, j, H) - P.w[0]);
        maxWerr = std::max(maxWerr, Werr);
        avWerr += Werr;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  gmaxUerr = GlobalMax(m_grid->com, maxUerr);
  gmaxWerr = GlobalMax(m_grid->com, maxWerr);
  gavUerr  = GlobalSum(m_grid->com, avUerr);
  gavUerr  = gavUerr/(m_grid->Mx()*m_grid->My());
  gavWerr  = GlobalSum(m_grid->com, avWerr);
  gavWerr  = gavWerr/(m_grid->Mx()*m_grid->My());
}


void IceCompModel::computeBasalMeltRateErrors(double &gmaxbmelterr, double &gminbmelterr) {
  double    maxbmelterr = -9.99e40, minbmelterr = 9.99e40, err;
  double    bmelt, dum1, dum2, dum3, dum4;

  if (testname != 'O') {
    throw RuntimeError("basal melt rate errors are only computable for test O");
  }

  // we just need one constant from exact solution:
  exactO(0.0, &dum1, &dum2, &dum3, &dum4, &bmelt);

  IceModelVec::AccessList list(m_basal_melt_rate);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    err = fabs(m_basal_melt_rate(i, j) - bmelt);
    maxbmelterr = std::max(maxbmelterr, err);
    minbmelterr = std::min(minbmelterr, err);
  }

  gmaxbmelterr = GlobalMax(m_grid->com, maxbmelterr);
  gminbmelterr = GlobalMin(m_grid->com, minbmelterr);
}


void IceCompModel::fillTemperatureSolnTestsKO() {

  double       dum1, dum2, dum3, dum4;
  std::vector<double> Tcol(m_grid->Mz());

  // evaluate exact solution in a column; all columns are the same
  switch (testname) {
    case 'K':
      for (unsigned int k=0; k<m_grid->Mz(); k++) {
        TestKParameters K = exactK(m_time->current(), m_grid->z(k), bedrock_is_ice_forK == true);
        Tcol[k] = K.T;
      }
      break;
    case 'O':
      for (unsigned int k=0; k<m_grid->Mz(); k++) {
        exactO(m_grid->z(k), &Tcol[k], &dum1, &dum2, &dum3, &dum4);
      }
      break;
    default:
      throw RuntimeError("only fills temperature solutions for tests K and O");
  }

  // copy column values into 3D arrays
  IceModelVec::AccessList list(m_ice_temperature);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      m_ice_temperature.set_column(i, j, &Tcol[0]);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  // communicate T
  m_ice_temperature.update_ghosts();
}


void IceCompModel::fillBasalMeltRateSolnTestO() {
  double       bmelt, dum1, dum2, dum3, dum4;
  if (testname != 'O') {
    throw RuntimeError("only fills basal melt rate soln for test O");
  }

  // we just need one constant from exact solution:
  exactO(0.0, &dum1, &dum2, &dum3, &dum4, &bmelt);

  m_basal_melt_rate.set(bmelt);
}


void IceCompModel::initTestsKO() {

  IceModelVec2S bed_topography;
  bed_topography.create(m_grid, "topg", WITHOUT_GHOSTS);
  bed_topography.set(0);
  m_beddef->set_elevation(bed_topography);

  m_ice_thickness.set(3000.0);
  m_ice_surface_elevation.copy_from(m_ice_thickness);

  fillTemperatureSolnTestsKO();
}

} // end of namespace pism
