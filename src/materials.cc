// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "materials.hh"

PetscScalar MaterialType::gasConst_R = 8.31441; // J/(mol K)    Gas Constant
PetscScalar MaterialType::grav  = 9.81;         // m/s^2        acceleration of gravity


PetscScalar BedrockType::k      = 3.0;          // J/(m K s) = W/(m K)    thermal conductivity
PetscScalar BedrockType::c_p    = 1000;         // J/(kg K)     specific heat capacity
// for following, reference Lingle & Clark (1985),  Bueler, Lingle, Kallen-Brown (2006)
// D = E T^3/(12 (1-nu^2)) for Young's modulus E = 6.6e10 N/m^2, lithosphere thickness T = 88 km,
//    and Poisson's ratio nu = 0.5
PetscScalar BedrockType::rho    = 3300;         // kg/(m^3)     density
PetscScalar BedrockType::D      = 5.0e24;       // N m          lithosphere flexural rigidity
PetscScalar BedrockType::eta    = 1.0e21;       // Pa s         half-space (mantle) viscosity

// DumbOceanType has nothing but density
PetscScalar DumbOceanType::rho      = 1027;     // kg/m         density


PetscScalar IceType::beta_CC_grad = 8.66e-4;// K/m          Clausius-Clapeyron gradient
PetscScalar IceType::rho    = 910;          // kg/m^3       density
PetscScalar IceType::k      = 2.10;         // J/(m K s) = W/(m K)    thermal conductivity
PetscScalar IceType::c_p    = 2009;         // J/(kg K)     specific heat capacity
PetscScalar IceType::latentHeat = 3.35e5;   // J/kg         latent heat capacity
PetscScalar IceType::meltingTemp = 273.15;   // K

PetscScalar IceType::flow(const PetscScalar stress, const PetscScalar temp,
                          const PetscScalar pressure) const {
  return 0;
}

PetscScalar IceType::flow(const PetscScalar stress, const PetscScalar temp,
                          const PetscScalar pressure, const PetscScalar gs) const {
  return flow(stress, temp, pressure);
}

PetscScalar IceType::effectiveViscosity(const PetscScalar regularization,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar temp, const PetscScalar pressure) const {
  const PetscScalar alpha = 0.5 * PetscSqr(u_x) + 0.5 * PetscSqr(v_y)
                             + 0.5 * PetscSqr(u_x + v_y)
                             + 0.25 * PetscSqr(u_y + v_x);
  const PetscScalar B = 135720960;
  return B / 2 * pow(regularization + alpha, -1.0 / 3.0);
}

PetscScalar IceType::effectiveViscosityColumn(const PetscScalar regularization,
                           const PetscScalar H, const PetscScalar dz,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar *T1, const PetscScalar *T2) const  {
  const PetscScalar alpha = 0.5 * PetscSqr(u_x) + 0.5 * PetscSqr(v_y)
                             + 0.5 * PetscSqr(u_x + v_y)
                             + 0.25 * PetscSqr(u_y + v_x);
  const PetscScalar B = 135720960;
  return H * B / 2 * pow(regularization + alpha, -1.0 / 3.0);
}

PetscScalar GlenIce::hardness_a = 135720960.;
PetscScalar GlenIce::softness_A = 4.0e-25;
PetscInt    GlenIce::n = 3;

PetscScalar GlenIce::flow(const PetscScalar stress, const PetscScalar temp,
                          const PetscScalar pressure) const {
  // ignor temp and pressure
  return softness_A * pow(stress, n-1);
}

PetscScalar GlenIce::exponent() const {
  return n;
}

// ThermoGlenIce is Paterson-Budd
PetscScalar ThermoGlenIce::A_cold = 3.61e-13;   // Pa^-3 / s
PetscScalar ThermoGlenIce::A_warm = 1.73e3; // Pa^-3 / s
PetscScalar ThermoGlenIce::Q_cold = 6.0e4;      // J / mol
PetscScalar ThermoGlenIce::Q_warm = 13.9e4; // J / mol
PetscScalar ThermoGlenIce::crit_temp = 263.15;  // K

PetscScalar ThermoGlenIce::flow(const PetscScalar stress, const PetscScalar temp,
                                const PetscScalar pressure) const {
  const PetscScalar T = temp + (beta_CC_grad / (rho * grav)) * pressure; // homologous temp
  return softnessParameter(T) * pow(stress,n-1);
}

PetscScalar ThermoGlenIce::effectiveViscosity(const PetscScalar regularization,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar temp, const PetscScalar pressure) const  {
  const PetscScalar T = temp + (beta_CC_grad / (rho * grav)) * pressure; // homologous temp
  const PetscScalar B = hardnessParameter(T);
  const PetscScalar alpha = 0.5 * PetscSqr(u_x) + 0.5 * PetscSqr(v_y)
                             + 0.5 * PetscSqr(u_x + v_y) + 0.25 * PetscSqr(u_y + v_x);
  return 0.5 * B * pow(regularization + alpha, -(n-1.0)/(2.0*n)); // regularization fixed at 1e-25 formerly
}

PetscScalar ThermoGlenIce::effectiveViscosityColumn(const PetscScalar regularization,
                           const PetscScalar H, const PetscScalar dz,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar *T1, const PetscScalar *T2) const  {
  // DESPITE NAME, does *not* return effective viscosity.
  // The result is \nu_e H, i.e. viscosity times thickness.
  // B is really hardness times thickness.
  const PetscInt  ks = static_cast<PetscInt>(floor(H/dz));
  // Integrate the hardness parameter using the trapezoid rule.
  PetscScalar B = 0;
  for (PetscInt m=1; m<ks; m++) {
    B += dz * hardnessParameter(0.5 * (T1[m] + T2[m]) + beta_CC_grad * (H - k*dz));
  }
  if (ks > 0) {
    B += 0.5 * dz * hardnessParameter(0.5 * (T1[0] + T2[0]) + beta_CC_grad * H);
    B += 0.5 * dz * hardnessParameter(0.5 * (T1[ks] + T2[ks]) + beta_CC_grad * (H - ks*dz));
  }
  const PetscScalar alpha = 0.5 * PetscSqr(u_x) + 0.5 * PetscSqr(v_y)
                             + 0.5 * PetscSqr(u_x + v_y) + 0.25 * PetscSqr(u_y + v_x);
  return 0.5 * B * pow(regularization + alpha, -(n-1.0)/(2.0*n));
}

PetscScalar ThermoGlenIce::softnessParameter(PetscScalar T) const {
  if (T < crit_temp) {
    return A_cold * exp(-Q_cold/(gasConst_R * T));
  }
  return A_warm * exp(-Q_warm/(gasConst_R * T));
}

PetscScalar ThermoGlenIce::hardnessParameter(PetscScalar T) const {
  return pow(softnessParameter(T), -1.0/n);
}


// ThermoGlenIceHooke: only change A(T) factor from ThermoGlenIce, which is
// Paterson-Budd
PetscScalar ThermoGlenIceHooke::Q_Hooke = 7.88e4;       // J / mol
// A_Hooke = (1/B_0)^n where n=3 and B_0 = 1.928 a^(1/3) Pa
PetscScalar ThermoGlenIceHooke::A_Hooke = 4.42165e-9;    // s^-1 Pa^-3
PetscScalar ThermoGlenIceHooke::C_Hooke = 0.16612;       // Kelvin^K_Hooke
PetscScalar ThermoGlenIceHooke::K_Hooke = 1.17;          // [unitless]
PetscScalar ThermoGlenIceHooke::Tr_Hooke = 273.39;       // Kelvin
PetscScalar ThermoGlenIceHooke::R_Hooke = 8.321;         // J mol^-1 K^-1

PetscScalar ThermoGlenIceHooke::softnessParameter(PetscScalar T) const {

  return A_Hooke * exp( -Q_Hooke/(R_Hooke * T)
                       + 3.0 * C_Hooke * pow(Tr_Hooke - T,-K_Hooke));
}


// ThermoGlenArrIce and ThermoGlenArrIceWarm: cold and warm parts, respectively,
// (i.e. simple Arrhenius) of Paterson-Budd
PetscScalar ThermoGlenArrIce::softnessParameter(PetscScalar T) const {
  return A() * exp(-Q()/(gasConst_R * T));
}

PetscScalar ThermoGlenArrIce::flow(const PetscScalar stress, const PetscScalar temp,
                                   const PetscScalar pressure) const {
  // ignors pressure
  return softnessParameter(temp) * pow(stress,n-1);  // uses NON-homologous temp
}

PetscScalar ThermoGlenArrIce::A() const {
  return A_cold;
}

PetscScalar ThermoGlenArrIce::Q() const {
  return Q_cold;
}

PetscScalar ThermoGlenArrIceWarm::A() const {
  return A_warm;
}

PetscScalar ThermoGlenArrIceWarm::Q() const {
  return Q_warm;
}


// HybridIce is Goldsby-Kohlstedt in ice sheets, Glen-Paterson-Budd in MacAyeal regions
PetscScalar HybridIce::V_act_vol    = -13.e-6;  // m^3/mol
PetscScalar HybridIce::d_grain_size = 1.0e-3;   // m  (see p. ?? of G&K paper)
//--- dislocation creep ---
PetscScalar
HybridIce::disl_crit_temp=258.0,    // Kelvin
  //disl_A_cold=4.0e5,                  // MPa^{-4.0} s^{-1}
  //disl_A_warm=6.0e28,                 // MPa^{-4.0} s^{-1}
  HybridIce::disl_A_cold=4.0e-19,     // Pa^{-4.0} s^{-1}
  HybridIce::disl_A_warm=6.0e4,       // Pa^{-4.0} s^{-1} (GK)
  HybridIce::disl_n=4.0,              // stress exponent
  HybridIce::disl_Q_cold=60.e3,       // J/mol Activation energy
  HybridIce::disl_Q_warm=180.e3;      // J/mol Activation energy (GK)
//--- grain boundary sliding ---
PetscScalar
HybridIce::gbs_crit_temp=255.0,     // Kelvin
  //gbs_A_cold=3.9e-3,                  // MPa^{-1.8} m^{1.4} s^{-1}
  //gbs_A_warm=3.e26,                   // MPa^{-1.8} m^{1.4} s^{-1}
  HybridIce::gbs_A_cold=6.1811e-14,   // Pa^{-1.8} m^{1.4} s^{-1}
  HybridIce::gbs_A_warm=4.7547e15,    // Pa^{-1.8} m^{1.4} s^{-1}
  HybridIce::gbs_n=1.8,               // stress exponent
  HybridIce::gbs_Q_cold=49.e3,        // J/mol Activation energy
  HybridIce::gbs_Q_warm=192.e3,       // J/mol Activation energy
  HybridIce::p_grain_sz_exp=1.4;      // from Peltier
//--- easy slip (basal) ---
PetscScalar
//basal_A=5.5e7,                      // MPa^{-2.4} s^{-1}
HybridIce::basal_A=2.1896e-7,       // Pa^{-2.4} s^{-1}
  HybridIce::basal_n=2.4,             // stress exponent
  HybridIce::basal_Q=60.e3;           // J/mol Activation energy
//--- diffusional flow ---
PetscScalar
HybridIce::diff_crit_temp=258.0,    // when to use enhancement factor
  HybridIce::diff_V_m=1.97e-5,        // Molar volume (m^3/mol)
  HybridIce::diff_D_0v=9.10e-4,       // Preexponential volume diffusion (m^2/s)
  HybridIce::diff_Q_v=59.4e3,         // activation energy, vol. diff. (J/mol)
  HybridIce::diff_D_0b=5.8e-4,        // preexponential grain boundary coeff.
  HybridIce::diff_Q_b=49.e3,          // activation energy, g.b. (J/mol)
  HybridIce::diff_delta=9.04e-10;     // grain boundary width (m)

PetscScalar HybridIce::flow(const PetscScalar stress, const PetscScalar temp,
                            const PetscScalar pressure) const {
  return flow(stress, temp, pressure, d_grain_size);
}

PetscScalar HybridIce::flow(const PetscScalar stress, const PetscScalar temp,
                            const PetscScalar pressure, const PetscScalar gs) const {
  /*
  This is the (forward) Goldsby-Kohlstedt flow law.  See:
  D. L. Goldsby & D. L. Kohlstedt (2001), "Superplastic deformation
  of ice: experimental observations", J. Geophys. Res. 106(M6), 11017-11030.
  */
  PetscScalar eps_diff, eps_disl, eps_basal, eps_gbs, diff_D_b;

  if (PetscAbs(stress) < 1e-10) return 0;
  const PetscScalar T = temp + (beta_CC_grad / (rho * grav)) * pressure;
  const PetscScalar pV = pressure * V_act_vol;
  const PetscScalar RT = gasConst_R * T;
  // Diffusional Flow
  const PetscScalar diff_D_v = diff_D_0v * exp(-diff_Q_v/RT);
  diff_D_b = diff_D_0b * exp(-diff_Q_b/RT);
  if (T > diff_crit_temp) diff_D_b *= 1000; // Coble creep scaling
  eps_diff = 14 * diff_V_m *
    (diff_D_v + M_PI * diff_delta * diff_D_b / gs) / (RT*PetscSqr(gs));
  // Dislocation Creep
  if (T > disl_crit_temp)
    eps_disl = disl_A_warm * pow(stress, disl_n-1) * exp(-(disl_Q_warm + pV)/RT);
  else
    eps_disl = disl_A_cold * pow(stress, disl_n-1) * exp(-(disl_Q_cold + pV)/RT);
  // Basal Slip
  eps_basal = basal_A * pow(stress, basal_n-1) * exp(-(basal_Q + pV)/RT);
  // Grain Boundary Sliding
  if (T > gbs_crit_temp)
    eps_gbs = gbs_A_warm * (pow(stress, gbs_n-1) / pow(gs, p_grain_sz_exp)) *
      exp(-(gbs_Q_warm + pV)/RT);
  else
    eps_gbs = gbs_A_cold * (pow(stress, gbs_n-1) / pow(gs, p_grain_sz_exp)) *
      exp(-(gbs_Q_cold + pV)/RT);

  return eps_diff + eps_disl + (eps_basal * eps_gbs) / (eps_basal + eps_gbs);
}

/*****************
THE NEXT PROCEDURE REPEATS CODE; INTENDED ONLY FOR DEBUGGING
*****************/
GKparts HybridIce::flowParts(const PetscScalar stress, const PetscScalar temp,
                             const PetscScalar pressure) const {
  PetscScalar gs, eps_diff, eps_disl, eps_basal, eps_gbs, diff_D_b;
  GKparts p;

  if (PetscAbs(stress) < 1e-10) {
    p.eps_total=0.0;
    p.eps_diff=0.0; p.eps_disl=0.0; p.eps_gbs=0.0; p.eps_basal=0.0;
    return p;
  }
  const PetscScalar T = temp + (beta_CC_grad / (rho * grav)) * pressure;
  const PetscScalar pV = pressure * V_act_vol;
  const PetscScalar RT = gasConst_R * T;
  // Diffusional Flow
  const PetscScalar diff_D_v = diff_D_0v * exp(-diff_Q_v/RT);
  diff_D_b = diff_D_0b * exp(-diff_Q_b/RT);
  if (T > diff_crit_temp) diff_D_b *= 1000; // Coble creep scaling
  gs = d_grain_size;
  eps_diff = 14 * diff_V_m *
    (diff_D_v + M_PI * diff_delta * diff_D_b / gs) / (RT*PetscSqr(gs));
  // Dislocation Creep
  if (T > disl_crit_temp)
    eps_disl = disl_A_warm * pow(stress, disl_n-1) * exp(-(disl_Q_warm + pV)/RT);
  else
    eps_disl = disl_A_cold * pow(stress, disl_n-1) * exp(-(disl_Q_cold + pV)/RT);
  // Basal Slip
  eps_basal = basal_A * pow(stress, basal_n-1) * exp(-(basal_Q + pV)/RT);
  // Grain Boundary Sliding
  if (T > gbs_crit_temp)
    eps_gbs = gbs_A_warm * (pow(stress, gbs_n-1) / pow(gs, p_grain_sz_exp)) *
      exp(-(gbs_Q_warm + pV)/RT);
  else
    eps_gbs = gbs_A_cold * (pow(stress, gbs_n-1) / pow(gs, p_grain_sz_exp)) *
      exp(-(gbs_Q_cold + pV)/RT);

  p.eps_diff=eps_diff;
  p.eps_disl=eps_disl;
  p.eps_basal=eps_basal;
  p.eps_gbs=eps_gbs;
  p.eps_total=eps_diff + eps_disl + (eps_basal * eps_gbs) / (eps_basal + eps_gbs);
  return p;
}
/*****************/


// HybridIceStripped is a simplification of Goldsby-Kohlstedt; compare that
// used in Peltier et al 2000, which is even simpler
PetscScalar HybridIceStripped::d_grain_size_stripped = 3.0e-3;
                                    // m; = 3mm  (see Peltier et al 2000 paper)

PetscScalar HybridIceStripped::flow(const PetscScalar stress, const PetscScalar temp,
                            const PetscScalar pressure) const {
  return flow(stress, temp, pressure, d_grain_size_stripped);
}

PetscScalar HybridIceStripped::flow(const PetscScalar stress, const PetscScalar temp,
                            const PetscScalar pressure, const PetscScalar gs) const {
  // note value of gs is ignored
  // note pressure only effects the temperature; the "P V" term is dropped
  // note no diffusional flow
  PetscScalar eps_disl, eps_basal, eps_gbs;

  if (PetscAbs(stress) < 1e-10) return 0;
  const PetscScalar T = temp + (beta_CC_grad / (rho * grav)) * pressure;
  const PetscScalar RT = gasConst_R * T;
  // NO Diffusional Flow
  // Dislocation Creep
  if (T > disl_crit_temp)
    eps_disl = disl_A_warm * pow(stress, disl_n-1) * exp(-disl_Q_warm/RT);
  else
    eps_disl = disl_A_cold * pow(stress, disl_n-1) * exp(-disl_Q_cold/RT);
  // Basal Slip
  eps_basal = basal_A * pow(stress, basal_n-1) * exp(-basal_Q/RT);
  // Grain Boundary Sliding
  if (T > gbs_crit_temp)
    eps_gbs = gbs_A_warm *
              (pow(stress, gbs_n-1) / pow(d_grain_size_stripped, p_grain_sz_exp)) *
              exp(-gbs_Q_warm/RT);
  else
    eps_gbs = gbs_A_cold *
              (pow(stress, gbs_n-1) / pow(d_grain_size_stripped, p_grain_sz_exp)) *
              exp(-gbs_Q_cold/RT);

  return eps_disl + (eps_basal * eps_gbs) / (eps_basal + eps_gbs);
}
