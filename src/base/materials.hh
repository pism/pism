// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

#ifndef __materials_hh
#define __materials_hh

#include <cmath>
#include <petsc.h>
#include "pism_const.hh"

/*******************
REGARDING IceType and HybridIce:  The hierarchy is:
  IceType <- GlenIce <- ThermoGlenIce <- HybridIce <- HybridIceStripped,
where "<-" means "derived class".  IceType is a virtual
class; it should never be used "as is".
Here:
GlenIce           means *isothermal* Glen flow law ice.
ThermoGlenIce     means *Paterson-Budd* version of Arhennius relation
HybridIce         means *Goldsby-Kohlstedt* flow law where vMask=SHEET,
                  and otherwise Paterson-Budd
HybridIceStripped means, where SHEET, G-K without the pressure dependence and without the
                  diffusional part; also grain size fixed at 3mm
Note each IceType has both a forward flow law ("flow") and an
inverted-and-vertically-integrated flow law ("effectiveViscosityColumn").  Only the
former form of the flow law is known for Goldsby-Kohlstedt.  If one can
invert-and-vertically-integrate the G-K law then one can build a "trueGKIce"
derived class.
*******************/

class IceType {
public:
  static PetscScalar beta_CC_grad;
  static PetscScalar rho;
  static PetscScalar k;
  static PetscScalar c_p;
  static PetscScalar latentHeat;
  static PetscScalar meltingTemp;

  virtual ~IceType() {}
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure, const PetscScalar gs) const;
  // this one returns nu; applies to ice shelf/stream approximation
  virtual PetscScalar effectiveViscosity(const PetscScalar regularization,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar temp, const PetscScalar pressure) const;
  // this one returns nu * H; it is adapted to a staggered grid so T1,T2 get averaged
  virtual PetscScalar effectiveViscosityColumn(const PetscScalar regularization,
                           const PetscScalar H, const PetscScalar dz,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar *T1, const PetscScalar *T2) const;
};


class GlenIce : public IceType {
public:
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar exponent() const;
protected:
  static PetscInt n;
  static PetscScalar hardness_a;
  static PetscScalar softness_A;
};


class ThermoGlenIce : public GlenIce {
public:
  static PetscScalar A_cold;  // these four constants from Paterson & Budd (1982)
  static PetscScalar A_warm;
  static PetscScalar Q_cold;
  static PetscScalar Q_warm;
  static PetscScalar crit_temp;
public:
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar effectiveViscosity(const PetscScalar regularization,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar temp, const PetscScalar pressure) const;
  virtual PetscScalar effectiveViscosityColumn(const PetscScalar regularization,
                           const PetscScalar H, const PetscScalar dz,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar *T1, const PetscScalar *T2) const;
protected:
  virtual PetscScalar softnessParameter(PetscScalar T) const;
  virtual PetscScalar hardnessParameter(PetscScalar T) const;
};

class ThermoGlenIceHooke : public ThermoGlenIce {
public:
  static PetscScalar A_Hooke;  // these six constants from Hooke (1981)
  static PetscScalar Q_Hooke;
  static PetscScalar C_Hooke;
  static PetscScalar K_Hooke;
  static PetscScalar Tr_Hooke;
  static PetscScalar R_Hooke;
protected:
  virtual PetscScalar softnessParameter(PetscScalar T) const;
};


class ThermoGlenArrIce : public ThermoGlenIce {
public:
  virtual PetscScalar softnessParameter(PetscScalar T) const;
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar A() const;  // returns A_cold for Paterson-Budd
  virtual PetscScalar Q() const;  // returns Q_cold for Paterson-Budd
};


class ThermoGlenArrIceWarm : public ThermoGlenArrIce {
public:
  virtual PetscScalar A() const;  // returns A_warm for Paterson-Budd
  virtual PetscScalar Q() const;  // returns Q_warm for Paterson-Budd
};



struct GKparts {
  PetscScalar eps_total, eps_diff, eps_disl, eps_basal, eps_gbs;
};

class HybridIce : public ThermoGlenIce {
public:
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure, const PetscScalar gs) const;
  GKparts flowParts(const PetscScalar stress, const PetscScalar temp,
                    const PetscScalar pressure) const;

protected:
  static PetscScalar  V_act_vol;
  static PetscScalar  d_grain_size;
  //--- diffusional flow ---
  static PetscScalar
  diff_crit_temp, diff_V_m, diff_D_0v, diff_Q_v, diff_D_0b, diff_Q_b, diff_delta;
  //--- dislocation creep ---
  static PetscScalar
  disl_crit_temp, disl_A_cold, disl_A_warm, disl_n, disl_Q_cold, disl_Q_warm;
  //--- easy slip (basal) ---
  static PetscScalar
  basal_A, basal_n, basal_Q;
  //--- grain boundary sliding ---
  static PetscScalar
  gbs_crit_temp, gbs_A_cold, gbs_A_warm, gbs_n, gbs_Q_cold, p_grain_sz_exp, gbs_Q_warm;
};


class HybridIceStripped : public HybridIce {
public:
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure, const PetscScalar gs) const;

protected:
  static PetscScalar d_grain_size_stripped;
};


class BedrockType {
public:
  static PetscScalar rho;
  static PetscScalar k;
  static PetscScalar c_p;
  static PetscScalar D;
  static PetscScalar eta;
};


class SeaWaterType {
public:
  static PetscScalar rho;
};


class FreshWaterType {
public:
  static PetscScalar rho;
};


class BasalType {
public:
  virtual PetscErrorCode printInfo(const int thresh, MPI_Comm com) = 0;
  virtual PetscScalar velocity(PetscScalar sliding_coefficient,
                               PetscScalar stress) = 0;
  virtual PetscScalar drag(PetscScalar beta, PetscScalar tauc,
                           PetscScalar vx, PetscScalar vy) = 0;
  virtual ~BasalType() {};
};

class ViscousBasalType : public BasalType {
public:
  virtual PetscErrorCode printInfo(const int thresh, MPI_Comm com);
  virtual PetscScalar velocity(PetscScalar sliding_coefficient,
                               PetscScalar stress);
  virtual PetscScalar drag(PetscScalar coeff, PetscScalar tauc,
                           PetscScalar vx, PetscScalar vy);
};


class PlasticBasalType : public ViscousBasalType {
public:
  PlasticBasalType(const PetscScalar pReg);
  virtual PetscErrorCode printInfo(const int thresh, MPI_Comm com);
  virtual PetscScalar drag(PetscScalar coeff, PetscScalar tauc,
                           PetscScalar vx, PetscScalar vy);
private:
  PetscScalar plastic_regularize;
};


#endif /* __materials_hh */
