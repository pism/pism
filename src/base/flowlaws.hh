// Copyright (C) 2004-2010 Jed Brown, Ed Bueler, and Constantine Khroulev
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

#ifndef __flowlaws_hh
#define __flowlaws_hh

#include <petsc.h>
#include "enthalpyConverter.hh"
#include "NCVariable.hh"


//! Abstract class containing the constitutive relation for the flow of ice.
/*!
This is the interface which most of PISM uses for rheology.
 */
class IceFlowLaw {
public:
  // ideally these would be protected or private:
  PetscScalar rho,          //!< density
              beta_CC_grad, //!< Clausius-Clapeyron gradient
              k,            //!< thermal conductivity
              c_p,          //!< specific heat capacity
              latentHeat,   //!< latent heat capacity
              meltingTemp,  //!< melting temperature
              standard_gravity,
    ideal_gas_constant,
    exponent_n;

  IceFlowLaw(MPI_Comm c,const char pre[], const NCConfigVariable &config);
  virtual ~IceFlowLaw();
  virtual PetscErrorCode setFromOptions() {return 0;}
  virtual PetscErrorCode printInfo(PetscInt) const;
  virtual PetscErrorCode view(PetscViewer) const {return 0;}

  // virtual PetscScalar flow_from_temp(PetscScalar stress, PetscScalar temp,
  //                          PetscScalar pressure, PetscScalar gs) const = 0;

  virtual PetscScalar flow_from_enth(PetscScalar stress, PetscScalar temp,
                           PetscScalar pressure, PetscScalar gs) const = 0;

  virtual PetscScalar effectiveViscosity(PetscScalar hardness,
                                         PetscScalar u_x, PetscScalar u_y,
                                         PetscScalar v_x, PetscScalar v_y) const = 0;

  virtual PetscScalar hardnessParameter(PetscScalar T) const = 0;
  virtual PetscScalar hardnessParameter_from_enth(PetscScalar E, PetscScalar p) const = 0;

  virtual PetscScalar averagedHardness_from_enth(PetscScalar thickness,
                                                 PetscInt kbelowH,
                                                 const PetscScalar *zlevels,
                                                 const PetscScalar *enthalpy) const;
  virtual PetscScalar exponent() const {return exponent_n;}
protected:
  MPI_Comm comm;
  char prefix[256];
  EnthalpyConverter *EC;
};

//! Glen (1955) and Paterson-Budd (1982) flow law with additional water fraction factor from Lliboutry & Duval (1985).
/*!
See [\ref AschwandenBlatter].  The basic references are [\ref Glen] and [\ref PatersonBudd] 
and [\ref LliboutryDuval1985].
 */
class GPBLDIce : public IceFlowLaw {
public:
  GPBLDIce(MPI_Comm c, const char pre[], const NCConfigVariable &config);
  virtual ~GPBLDIce() {}

  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode view(PetscViewer viewer) const;

  virtual PetscScalar softnessParameter_from_enth(
              PetscScalar enthalpy, PetscScalar pressure) const;
  virtual PetscScalar hardnessParameter_from_enth(
              PetscScalar enthalpy, PetscScalar pressure) const;

  virtual PetscScalar flow_from_enth(
              PetscScalar stress, PetscScalar enthalpy, PetscScalar pressure,
              PetscScalar gs) const; // grainsize arg gs not used

  virtual PetscScalar effectiveViscosity(PetscScalar hardness,
                                         PetscScalar u_x, PetscScalar u_y,
                                         PetscScalar v_x, PetscScalar v_y) const;

protected:
  PetscReal T_0, water_frac_coeff;
  PetscReal schoofLen,schoofVel,schoofReg,
    A_cold, A_warm, Q_cold, Q_warm,  // see Paterson & Budd (1982)
    crit_temp, n;

  virtual PetscScalar softnessParameter_pb(PetscScalar T_pa) const;
  virtual PetscScalar hardnessParameter_pb(PetscScalar T_pa) const;
};

// This uses the definition of second invariant from Hutter and several others, namely
// \f$ \frac 1 2 D_{ij} D_{ij} \f$ where incompressibility is used to compute \f$ D_{zz} \f$
static inline PetscScalar secondInvariant(PetscScalar u_x, PetscScalar u_y,
                                          PetscScalar v_x, PetscScalar v_y)
{ return 0.5 * (PetscSqr(u_x) + PetscSqr(v_y) + PetscSqr(u_x + v_y) + 0.5*PetscSqr(u_y + v_x)); }

// The second invariant of a symmetric strain rate tensor in compressed form [u_x, v_y, 0.5(u_y+v_x)]
static inline PetscScalar secondInvariantDu(const PetscScalar Du[])
{ return 0.5 * (PetscSqr(Du[0]) + PetscSqr(Du[1]) + PetscSqr(Du[0]+Du[1]) + 2*PetscSqr(Du[2])); }

#endif // __flowlaws_hh
