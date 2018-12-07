// Copyright (C) 2004-2018 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/Time.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/error_handling.hh"
#include "pism/earth/BedDef.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/pism_utilities.hh"
#include "BTU_Verification.hh"
#include "pism/energy/TemperatureModel.hh"
#include "pism/coupler/SurfaceModel.hh"
#include "pism/coupler/OceanModel.hh"
#include "pism/hydrology/Hydrology.hh"

namespace pism {

// boundary conditions for tests F, G (same as EISMINT II Experiment F)
const double IceCompModel::m_ST     = 1.67e-5;
const double IceCompModel::m_Tmin   = 223.15; // K
const double IceCompModel::m_LforFG = 750000; // m
const double IceCompModel::m_ApforG = 200; // m

void IceCompModel::energy_step() {

  energy::EnergyModelStats stats;

  IceModelVec2S &bedtoptemp              = m_work2d[1];
  IceModelVec2S &basal_enthalpy          = m_work2d[2];
  m_energy_model->enthalpy().getHorSlice(basal_enthalpy, 0.0);

  bedrock_surface_temperature(m_geometry.sea_level_elevation,
                              m_geometry.cell_type,
                              m_geometry.bed_elevation,
                              m_geometry.ice_thickness,
                              basal_enthalpy,
                              m_surface->temperature(),
                              bedtoptemp);

  m_btu->update(bedtoptemp, t_TempAge, dt_TempAge);

  energy::Inputs inputs;
  {
    inputs.basal_frictional_heating = &m_stress_balance->basal_frictional_heating();
    inputs.basal_heat_flux          = &m_btu->flux_through_top_surface(); // bedrock thermal layer
    inputs.cell_type                = &m_geometry.cell_type;              // geometry
    inputs.ice_thickness            = &m_geometry.ice_thickness;          // geometry
    inputs.shelf_base_temp          = &m_ocean->shelf_base_temperature(); // ocean model
    inputs.surface_liquid_fraction  = &m_surface->liquid_water_fraction(); // surface model
    inputs.surface_temp             = &m_surface->temperature(); // surface model
    inputs.till_water_thickness     = &m_subglacial_hydrology->till_water_thickness();

    inputs.volumetric_heating_rate  = &m_stress_balance->volumetric_strain_heating();
    inputs.u3                       = &m_stress_balance->velocity_u();
    inputs.v3                       = &m_stress_balance->velocity_v();
    inputs.w3                       = &m_stress_balance->velocity_w();

    inputs.check();             // make sure all data members were set
  }

  if ((m_testname == 'F') || (m_testname == 'G')) {
    // Compute compensatory strain heating (fills strain_heating3_comp).
    getCompSourcesTestFG();

    // Add computed strain heating to the compensatory part.
    m_strain_heating3_comp.add(1.0, *inputs.volumetric_heating_rate);

    // Use the result.
    inputs.volumetric_heating_rate = &m_strain_heating3_comp;
  }

  m_energy_model->update(t_TempAge, dt_TempAge, inputs);

  m_stdout_flags = m_energy_model->stdout_flags() + m_stdout_flags;
}

void IceCompModel::initTestFG() {

  IceModelVec::AccessList list{&m_geometry.ice_thickness};

  const double time = m_testname == 'F' ? 0.0 : m_time->current();
  const double A    = m_testname == 'F' ? 0.0 : m_ApforG;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // avoid singularity at origin
    const double r = std::max(radius(*m_grid, i, j), 1.0);

    if (r > m_LforFG - 1.0) { // if (essentially) outside of sheet
      m_geometry.ice_thickness(i, j) = 0.0;
    } else {
      m_geometry.ice_thickness(i, j) = exactFG(time, r, m_grid->z(), A).H;
    }
  }

  m_geometry.ice_thickness.update_ghosts();

  {
    IceModelVec2S bed_topography(m_grid, "topg", WITHOUT_GHOSTS);
    bed_topography.set(0.0);

    IceModelVec2S bed_uplift(m_grid, "uplift", WITHOUT_GHOSTS);
    bed_uplift.set(0.0);

    IceModelVec2S sea_level(m_grid, "sea_level", WITHOUT_GHOSTS);
    sea_level.set(0.0);

    m_beddef->bootstrap(bed_topography, bed_uplift, m_geometry.ice_thickness,
                        sea_level);
  }
}

void IceCompModel::initTestsKO() {

  IceModelVec2S bed_topography(m_grid, "topg", WITHOUT_GHOSTS);
  bed_topography.set(0.0);

  IceModelVec2S bed_uplift(m_grid, "uplift", WITHOUT_GHOSTS);
  bed_uplift.set(0.0);

  IceModelVec2S sea_level(m_grid, "sea_level", WITHOUT_GHOSTS);
  sea_level.set(0.0);

  m_geometry.ice_thickness.set(3000.0);

  m_beddef->bootstrap(bed_topography, bed_uplift, m_geometry.ice_thickness,
                      sea_level);
}

void IceCompModel::getCompSourcesTestFG() {

  const double
    ice_rho   = m_config->get_double("constants.ice.density"),
    ice_c     = m_config->get_double("constants.ice.specific_heat_capacity");

  // before temperature and flow step, set strain_heating_c from exact values

  IceModelVec::AccessList list{&m_strain_heating3_comp};

  const double time = m_testname == 'F' ? 0.0 : m_time->current();
  const double A    = m_testname == 'F' ? 0.0 : m_ApforG;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = std::max(radius(*m_grid, i, j), 1.0); // avoid singularity at origin

    if (r > m_LforFG - 1.0) {  // outside of sheet
      m_strain_heating3_comp.set_column(i, j, 0.0);
    } else {
      TestFGParameters P = exactFG(time, r, m_grid->z(), A);

      m_strain_heating3_comp.set_column(i, j, &P.Sigc[0]);
    }
  }

  // scale strain_heating to J/(s m^3)
  m_strain_heating3_comp.scale(ice_rho * ice_c);
}

void IceCompModel::computeTemperatureErrors(double &gmaxTerr,
                                            double &gavTerr) {
  double maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;

  if (m_testname != 'F' and m_testname != 'G') {
    throw RuntimeError(PISM_ERROR_LOCATION, "temperature errors only computable for tests F and G");
  }

  const double time = m_testname == 'F' ? 0.0 : m_time->current();
  const double A    = m_testname == 'F' ? 0.0 : m_ApforG;

  energy::TemperatureModel *m = dynamic_cast<energy::TemperatureModel*>(m_energy_model);
  const IceModelVec3 &ice_temperature = m->temperature();

  IceModelVec::AccessList list{&m_geometry.ice_thickness, &ice_temperature};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double r = radius(*m_grid, i, j);
      const double *T = ice_temperature.get_column(i, j);

      // only evaluate error if inside sheet and not at central
      // singularity
      if ((r >= 1.0) and (r <= m_LforFG - 1.0)) {
        TestFGParameters P = exactFG(time, r, m_grid->z(), A);

        // only evaluate error if below ice surface
        const int ks = m_grid->kBelowHeight(m_geometry.ice_thickness(i, j));
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

  if (m_testname != 'K' and m_testname != 'O') {
    throw RuntimeError(PISM_ERROR_LOCATION, "ice and bedrock temperature errors only computable for tests K and O");
  }

  double    maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;
  double    maxTberr = 0.0, avTberr = 0.0, avbcount = 0.0;

  const double *Tb, *T;
  std::vector<double> Tex(m_grid->Mz());

  energy::BTU_Verification *my_btu = dynamic_cast<energy::BTU_Verification*>(m_btu);
  if (my_btu == NULL) {
    throw RuntimeError(PISM_ERROR_LOCATION, "BTU_Verification is required");
  }
  const IceModelVec3Custom &bedrock_temp = my_btu->temperature();

  std::vector<double> zblevels = bedrock_temp.levels();
  unsigned int Mbz = (unsigned int)zblevels.size();
  std::vector<double> Tbex(Mbz);

  switch (m_testname) {
    case 'K':
      for (unsigned int k = 0; k < m_grid->Mz(); k++) {
        TestKParameters K = exactK(m_time->current(), m_grid->z(k), m_bedrock_is_ice_forK);
        Tex[k] = K.T;
      }
      for (unsigned int k = 0; k < Mbz; k++) {
        TestKParameters K = exactK(m_time->current(), zblevels[k], m_bedrock_is_ice_forK);
        Tbex[k] = K.T;
      }
      break;
    case 'O':
      for (unsigned int k = 0; k < m_grid->Mz(); k++) {
        Tex[k] = exactO(m_grid->z(k)).TT;
      }
      for (unsigned int k = 0; k < Mbz; k++) {
        Tbex[k] = exactO(zblevels[k]).TT;
      }
      break;
    default:
      throw RuntimeError(PISM_ERROR_LOCATION, "ice and bedrock temperature errors only for tests K and O");
  }

  energy::TemperatureModel *m = dynamic_cast<energy::TemperatureModel*>(m_energy_model);
  const IceModelVec3 &ice_temperature = m->temperature();

  IceModelVec::AccessList list{&ice_temperature, &bedrock_temp};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    Tb = bedrock_temp.get_column(i, j);
    for (unsigned int kb = 0; kb < Mbz; kb++) {
      const double Tberr = fabs(Tb[kb] - Tbex[kb]);
      maxTberr = std::max(maxTberr, Tberr);
      avbcount += 1.0;
      avTberr += Tberr;
    }
    T = ice_temperature.get_column(i, j);
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

  if (m_testname != 'F' and m_testname != 'G') {
    throw RuntimeError(PISM_ERROR_LOCATION, "temperature errors only computable for tests F and G");
  }

  double
    Texact     = 0.0,
    domeT      = 0.0,
    domeTexact = 0.0,
    Terr       = 0.0,
    avTerr     = 0.0;

  const double time = m_testname == 'F' ? 0.0 : m_time->current();
  const double A    = m_testname == 'F' ? 0.0 : m_ApforG;
  std::vector<double> z(1, 0.0);

  energy::TemperatureModel *m = dynamic_cast<energy::TemperatureModel*>(m_energy_model);
  const IceModelVec3 &ice_temperature = m->temperature();

  IceModelVec::AccessList list(ice_temperature);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double r = std::max(radius(*m_grid, i, j), 1.0);

      if (r > m_LforFG - 1.0) { // outside of sheet
        Texact = m_Tmin + m_ST * r; // = Ts
      } else {
        Texact = exactFG(time, r, z, A).T[0];
      }

      const double Tbase = ice_temperature.get_column(i,j)[0];
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

  if (m_testname != 'F' and m_testname != 'G') {
    throw RuntimeError(PISM_ERROR_LOCATION, "strain-heating (strain_heating) errors only computable for tests F and G");
  }

  const double
    ice_rho   = m_config->get_double("constants.ice.density"),
    ice_c     = m_config->get_double("constants.ice.specific_heat_capacity");

  const IceModelVec3 &strain_heating3 = m_stress_balance->volumetric_strain_heating();

  IceModelVec::AccessList list{&m_geometry.ice_thickness, &strain_heating3};

  const double time = m_testname == 'F' ? 0.0 : m_time->current();
  const double A    = m_testname == 'F' ? 0.0 : m_ApforG;

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double r = radius(*m_grid, i, j);
      if ((r >= 1.0) && (r <= m_LforFG - 1.0)) {
        // only evaluate error if inside sheet and not at central singularity

        TestFGParameters P = exactFG(time, r, m_grid->z(), A);

        for (unsigned int k = 0; k < m_grid->Mz(); k++) {
          // scale exact strain_heating to J/(s m^3)
          P.Sig[k] *= ice_rho * ice_c;
        }

        const unsigned int ks = m_grid->kBelowHeight(m_geometry.ice_thickness(i, j));
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

  IceModelVec::AccessList list{&m_geometry.ice_thickness, &u3, &v3, &w3};

  const double time = m_testname == 'F' ? 0.0 : m_time->current();
  const double A    = m_testname == 'F' ? 0.0 : m_ApforG;

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double H = m_geometry.ice_thickness(i, j);
      std::vector<double> z(1, H);

      const double
        x = m_grid->x(i),
        y = m_grid->y(j),
        r = radius(*m_grid, i, j);

      if ((r >= 1.0) and (r <= m_LforFG - 1.0)) {
        // only evaluate error if inside sheet and not at central singularity

        TestFGParameters P = exactFG(time, r, z, A);

        const double
          uex = (x/r) * P.U[0],
          vex = (y/r) * P.U[0];
        // note that because getValZ does linear interpolation and H(i, j) is not exactly at
        // a grid point, this causes nonzero errors
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
  double
    maxbmelterr = -9.99e40,
    minbmelterr = 9.99e40;

  if (m_testname != 'O') {
    throw RuntimeError(PISM_ERROR_LOCATION, "basal melt rate errors are only computable for test O");
  }

  // we just need one constant from exact solution:
  double bmelt = exactO(0.0).bmelt;

  const IceModelVec2S &basal_melt_rate = m_energy_model->basal_melt_rate();

  IceModelVec::AccessList list(basal_melt_rate);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double err = fabs(basal_melt_rate(i, j) - bmelt);
    maxbmelterr = std::max(maxbmelterr, err);
    minbmelterr = std::min(minbmelterr, err);
  }

  gmaxbmelterr = GlobalMax(m_grid->com, maxbmelterr);
  gminbmelterr = GlobalMin(m_grid->com, minbmelterr);
}

} // end of namespace pism
