// Copyright (C) 2011-2018 Torsten Albrecht and Constantine Khroulev
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


#include <petscsys.h>           // PetscSqr
#include <cmath>                // atan, sqrt
#include <gsl/gsl_math.h>       // M_PI

#include "IceModel.hh"

#include "pism/energy/EnergyModel.hh"

#include "pism/rheology/FlowLaw.hh"

#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_options.hh"

namespace pism {

//! \file fracture_density.cc implementing calculation of fracture density with PIK options -fractures.

void IceModel::update_fracture_density() {
  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();
  const int
    Mx = m_grid->Mx(),
    My = m_grid->My();

  IceModelVec2S
    &D        = m_fracture->density,
    &A        = m_fracture->age,
    &D_new    = m_work2d[0],
    &A_new    = m_work2d[1],
    &hardness = m_work2d[2];

  // get SSA velocities and related strain rates and stresses
  const IceModelVec2V &ssa_velocity        = m_stress_balance->advective_velocity();
  IceModelVec2        &strain_rates        = m_fracture->strain_rates;
  IceModelVec2        &deviatoric_stresses = m_fracture->deviatoric_stresses;

  stressbalance::compute_2D_principal_strain_rates(ssa_velocity,
                                                   m_geometry.cell_type,
                                                   strain_rates);

  averaged_hardness_vec(*m_stress_balance->shallow()->flow_law(),
                        m_geometry.ice_thickness,
                        m_energy_model->enthalpy(),
                        hardness);

  m_stress_balance->compute_2D_stresses(ssa_velocity,
                                        hardness,
                                        m_geometry.cell_type,
                                        deviatoric_stresses);

  IceModelVec::AccessList list{&ssa_velocity, &strain_rates, &deviatoric_stresses,
      &m_geometry.ice_thickness, &D, &D_new, &m_geometry.cell_type};

  D_new.copy_from(D);

  const bool dirichlet_bc = m_config->get_boolean("stress_balance.ssa.dirichlet_bc");
  if (dirichlet_bc) {
    list.add(m_ssa_dirichlet_bc_mask);
    list.add(m_ssa_dirichlet_bc_values);
  }

  const bool write_fd = m_config->get_boolean("fracture_density.write_fields");
  if (write_fd) {
    list.add({&m_fracture->growth_rate, &m_fracture->healing_rate,
          &m_fracture->flow_enhancement, &m_fracture->toughness, &A, &A_new});
    A_new.copy_from(A);
  }

  //options
  /////////////////////////////////////////////////////////
  double soft_residual = options::Real("-fracture_softening", "soft_residual", 1.0);
  // assume linear response function: E_fr = (1-(1-soft_residual)*phi) -> 1-phi
  //
  // more: T. Albrecht, A. Levermann; Fracture-induced softening for
  // large-scale ice dynamics; (2013), The Cryosphere Discussions 7;
  // 4501-4544; DOI:10.5194/tcd-7-4501-2013

  // get four options for calculation of fracture density.
  // 1st: fracture growth constant gamma
  // 2nd: fracture initiation stress threshold sigma_cr
  // 3rd: healing rate constant gamma_h
  // 4th: healing strain rate threshold
  // more: T. Albrecht, A. Levermann; Fracture field for large-scale
  // ice dynamics; (2012), Journal of Glaciology, Vol. 58, No. 207,
  // 165-176, DOI: 10.3189/2012JoG11J191.

  double gamma = 1.0, initThreshold = 7.0e4, gammaheal = 0.0, healThreshold = 2.0e-10;
  {
    options::RealList fractures("-fracture_parameters",
                                "gamma, initThreshold, gammaheal, healThreshold",
                                {gamma, initThreshold, gammaheal, healThreshold});
    if (fractures->size() != 4) {
      throw RuntimeError(PISM_ERROR_LOCATION, "option -fracture_parameters requires exactly 4 arguments");
    }
    gamma         = fractures[0];
    initThreshold = fractures[1];
    gammaheal     = fractures[2];
    healThreshold = fractures[3];
  }

  m_log->message(3, "PISM-PIK INFO: fracture density is found with parameters:\n"
                    " gamma=%.2f, sigma_cr=%.2f, gammah=%.2f, healing_cr=%.1e and soft_res=%f \n",
                 gamma, initThreshold, gammaheal, healThreshold, soft_residual);

  bool do_fracground = m_config->get_boolean("fracture_density.include_grounded_ice");

  double fdBoundaryValue = m_config->get_double("fracture_density.phi0");

  bool constant_healing = m_config->get_boolean("fracture_density.constant_healing");

  bool fracture_weighted_healing = m_config->get_boolean("fracture_density.fracture_weighted_healing");

  bool max_shear_stress = m_config->get_boolean("fracture_density.max_shear_stress");

  bool lefm = m_config->get_boolean("fracture_density.lefm");

  bool constant_fd = m_config->get_boolean("fracture_density.constant_fd");

  bool fd2d_scheme = m_config->get_boolean("fracture_density.fd2d_scheme");

  const double one_year = units::convert(m_sys, 1.0, "year", "seconds");

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double tempFD = 0.0;
    //SSA: v . grad memField

    double uvel = ssa_velocity(i, j).u;
    double vvel = ssa_velocity(i, j).v;

    if (fd2d_scheme) {
      if (uvel >= dx * vvel / dy && vvel >= 0.0) { //1
        tempFD = uvel * (D(i, j) - D(i - 1, j)) / dx + vvel * (D(i - 1, j) - D(i - 1, j - 1)) / dy;
      } else if (uvel <= dx * vvel / dy && uvel >= 0.0) { //2
        tempFD = uvel * (D(i, j - 1) - D(i - 1, j - 1)) / dx + vvel * (D(i, j) - D(i, j - 1)) / dy;
      } else if (uvel >= -dx * vvel / dy && uvel <= 0.0) { //3
        tempFD = -uvel * (D(i, j - 1) - D(i + 1, j - 1)) / dx + vvel * (D(i, j) - D(i, j - 1)) / dy;
      } else if (uvel <= -dx * vvel / dy && vvel >= 0.0) { //4
        tempFD = -uvel * (D(i, j) - D(i + 1, j)) / dx + vvel * (D(i + 1, j) - D(i + 1, j - 1)) / dy;
      } else if (uvel <= dx * vvel / dy && vvel <= 0.0) { //5
        tempFD = -uvel * (D(i, j) - D(i + 1, j)) / dx - vvel * (D(i + 1, j) - D(i + 1, j + 1)) / dy;
      } else if (uvel >= dx * vvel / dy && uvel <= 0.0) { //6
        tempFD = -uvel * (D(i, j + 1) - D(i + 1, j + 1)) / dx - vvel * (D(i, j) - D(i, j + 1)) / dy;
      } else if (uvel <= -dx * vvel / dy && uvel >= 0.0) { //7
        tempFD = uvel * (D(i, j + 1) - D(i - 1, j + 1)) / dx - vvel * (D(i, j) - D(i, j + 1)) / dy;
      } else if (uvel >= -dx * vvel / dy && vvel <= 0.0) { //8
        tempFD = uvel * (D(i, j) - D(i - 1, j)) / dx - vvel * (D(i - 1, j) - D(i - 1, j + 1)) / dy;
      } else {
        m_log->message(3, "######### missing case of angle %f of %f and %f at %d, %d \n",
                       atan(vvel / uvel) / M_PI * 180., uvel * 3e7, vvel * 3e7, i, j);
      }
    } else {
      tempFD += uvel * (uvel < 0 ? D(i + 1, j) - D(i, j) : D(i, j) - D(i - 1, j)) / dx;
      tempFD += vvel * (vvel < 0 ? D(i, j + 1) - D(i, j) : D(i, j) - D(i, j - 1)) / dy;
    }

    D_new(i, j) -= tempFD * m_dt;

    //sources /////////////////////////////////////////////////////////////////
    ///von mises criterion

    double
      txx    = deviatoric_stresses(i, j, 0),
      tyy    = deviatoric_stresses(i, j, 1),
      txy    = deviatoric_stresses(i, j, 2),
      T1     = 0.5 * (txx + tyy) + sqrt(0.25 * PetscSqr(txx - tyy) + PetscSqr(txy)), //Pa
      T2     = 0.5 * (txx + tyy) - sqrt(0.25 * PetscSqr(txx - tyy) + PetscSqr(txy)), //Pa
      sigmat = sqrt(PetscSqr(T1) + PetscSqr(T2) - T1 * T2);


    ///max shear stress criterion (more stringent than von mises)
    if (max_shear_stress) {
      double maxshear = fabs(T1);
      maxshear        = std::max(maxshear, fabs(T2));
      maxshear        = std::max(maxshear, fabs(T1 - T2));

      sigmat = maxshear;
    }

    ///lefm mixed-mode criterion
    if (lefm) {
      double sigmamu = 0.1; //friction coefficient between crack faces

      double sigmac = 0.64 / M_PI; //initial crack depth 20cm

      double sigmabetatest, sigmanor, sigmatau, Kone, Ktwo, KSI, KSImax = 0.0, sigmatetanull;

      for (int l = 46; l <= 90; ++l) { //optimize for various precursor angles beta
        sigmabetatest = l * M_PI / 180.0;

        //rist_sammonds99
        sigmanor = 0.5 * (T1 + T2) - (T1 - T2) * cos(2 * sigmabetatest);
        sigmatau = 0.5 * (T1 - T2) * sin(2 * sigmabetatest);
        //shayam_wu90
        if (sigmamu * sigmanor < 0.0) { //compressive case
          if (fabs(sigmatau) <= fabs(sigmamu * sigmanor)) {
            sigmatau = 0.0;
          } else {
            if (sigmatau > 0) { //coulomb friction opposing sliding
              sigmatau += (sigmamu * sigmanor);
            } else {
              sigmatau -= (sigmamu * sigmanor);
            }
          }
        }

        //stress intensity factors
        Kone = sigmanor * sqrt(M_PI * sigmac); //normal
        Ktwo = sigmatau * sqrt(M_PI * sigmac); //shear

        if (Ktwo == 0.0) {
          sigmatetanull = 0.0;
        } else { //eq15 in hulbe_ledoux10 or eq15 shayam_wu90
          sigmatetanull = -2.0 * atan((sqrt(PetscSqr(Kone) + 8.0 * PetscSqr(Ktwo)) - Kone) / (4.0 * Ktwo));
        }

        KSI = cos(0.5 * sigmatetanull) *
              (Kone * cos(0.5 * sigmatetanull) * cos(0.5 * sigmatetanull) - 0.5 * 3.0 * Ktwo * sin(sigmatetanull));
        // mode I stress intensity

        KSImax = std::max(KSI, KSImax);
      }
      sigmat = KSImax;
    }

    //////////////////////////////////////////////////////////////////////////////

    //fracture density
    double fdnew = gamma * (strain_rates(i, j, 0) - 0.0) * (1 - D_new(i, j));
    if (sigmat > initThreshold) {
      D_new(i, j) += fdnew * m_dt;
    }

    //healing
    double fdheal = gammaheal * (strain_rates(i, j, 0) - healThreshold);
    if (m_geometry.ice_thickness(i, j) > 0.0) {
      if (constant_healing) {
        fdheal = gammaheal * (-healThreshold);
        if (fracture_weighted_healing) {
          D_new(i, j) += fdheal * m_dt * (1 - D(i, j));
        } else {
          D_new(i, j) += fdheal * m_dt;
        }
      } else if (strain_rates(i, j, 0) < healThreshold) {
        if (fracture_weighted_healing) {
          D_new(i, j) += fdheal * m_dt * (1 - D(i, j));
        } else {
          D_new(i, j) += fdheal * m_dt;
        }
      }
    }

    //bounding
    if (D_new(i, j) < 0.0) {
      D_new(i, j) = 0.0;
    }

    if (D_new(i, j) > 1.0) {
      D_new(i, j) = 1.0;
    }

    //################################################################################
    // write related fracture quantities to nc-file
    // if option -write_fd_fields is set
    if (write_fd && m_geometry.ice_thickness(i, j) > 0.0) {
      //fracture toughness
      m_fracture->toughness(i, j) = sigmat;

      // fracture growth rate
      if (sigmat > initThreshold) {
        m_fracture->growth_rate(i, j) = fdnew;
        //m_fracture->growth_rate(i,j)=gamma*(vPrinStrain1(i,j)-0.0)*(1-D_new(i,j));
      } else {
        m_fracture->growth_rate(i, j) = 0.0;
      }

      // fracture healing rate
      if (m_geometry.ice_thickness(i, j) > 0.0) {
        if (constant_healing || (strain_rates(i, j, 0) < healThreshold)) {
          if (fracture_weighted_healing) {
            m_fracture->healing_rate(i, j) = fdheal * (1 - D(i, j));
          } else {
            m_fracture->healing_rate(i, j) = fdheal;
          }
        } else {
          m_fracture->healing_rate(i, j) = 0.0;
        }
      }

      //fracture age since fracturing occured
      A_new(i, j) -= m_dt * uvel * (uvel < 0 ? A(i + 1, j) - A(i, j) : A(i, j) - A(i - 1, j)) / dx;
      A_new(i, j) -= m_dt * vvel * (vvel < 0 ? A(i, j + 1) - A(i, j) : A(i, j) - A(i, j - 1)) / dy;
      A_new(i, j) += m_dt / one_year;
      if (sigmat > initThreshold) {
        A_new(i, j) = 0.0;
      }

      // additional flow enhancement due to fracture softening
      double phi_exp   = 3.0; //flow_law->exponent();
      double softening = pow((1.0 - (1.0 - soft_residual) * D_new(i, j)), -phi_exp);
      if (m_geometry.ice_thickness(i, j) > 0.0) {
        m_fracture->flow_enhancement(i, j) = 1.0 / pow(softening, 1 / 3.0);
      } else {
        m_fracture->flow_enhancement(i, j) = 1.0;
      }
    }

    //boundary condition
    if (dirichlet_bc && !do_fracground) {
      if (m_ssa_dirichlet_bc_mask.as_int(i, j) == 1) {
        if (m_ssa_dirichlet_bc_values(i, j).u != 0.0 || m_ssa_dirichlet_bc_values(i, j).v != 0.0) {
          D_new(i, j) = fdBoundaryValue;
        }

        if (write_fd) {
          A_new(i, j)                       = 0.0;
          m_fracture->growth_rate(i, j)      = 0.0;
          m_fracture->healing_rate(i, j)     = 0.0;
          m_fracture->flow_enhancement(i, j) = 1.0;
          m_fracture->toughness(i, j)        = 0.0;
        }
      }
    }
    // ice free regions and boundary of computational domain
    if (m_geometry.ice_thickness(i, j) == 0.0 || i == 0 || j == 0 || i == Mx - 1 || j == My - 1) {
      D_new(i, j) = 0.0;
      if (write_fd) {
        A_new(i, j)                       = 0.0;
        m_fracture->growth_rate(i, j)      = 0.0;
        m_fracture->healing_rate(i, j)     = 0.0;
        m_fracture->flow_enhancement(i, j) = 1.0;
        m_fracture->toughness(i, j)        = 0.0;
      }
    }

    if (constant_fd) { //no fd evolution
      D_new(i, j) = D(i, j);
    }
  }

  if (write_fd) {
    A_new.update_ghosts(A);
  }

  D_new.update_ghosts(D);
}

} // end of namespace pism
