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

#include "flowlaws.hh"
#include "pism_const.hh"
#include "enthalpyConverter.hh"

IceFlowLaw::IceFlowLaw(MPI_Comm c,const char pre[], const NCConfigVariable &config) : comm(c) {
  PetscMemzero(prefix,sizeof(prefix));
  if (pre) PetscStrncpy(prefix,pre,sizeof(prefix));

  standard_gravity = config.get("standard_gravity");
  ideal_gas_constant = config.get("ideal_gas_constant");
  rho          = config.get("ice_density");
  beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * standard_gravity;
  k            = config.get("ice_thermal_conductivity");
  c_p          = config.get("ice_specific_heat_capacity");
  latentHeat   = config.get("water_latent_heat_fusion");
  meltingTemp  = config.get("water_melting_temperature");

  if (config.get_flag("verification_mode")) {
    EC = new ICMEnthalpyConverter(config);
  } else {
    EC = new EnthalpyConverter(config);
  }
}

IceFlowLaw::~IceFlowLaw() {
  delete EC;
}

PetscErrorCode IceFlowLaw::printInfo(PetscInt) const {
  PetscPrintf(comm,"WARNING:  IceFlowLaw::printInfo() called but base class IceFlowLaw is virtual!!\n");
  return 0;
}

//! Computes vertical average of B(E,pressure) ice hardness, namely \f$\bar
//! B(E,p)\f$. See comment for hardnessParameter().
/*! Note E[0],...,E[kbelowH] must be valid.  */
PetscScalar IceFlowLaw::averagedHardness_from_enth(
                PetscScalar thickness, PetscInt kbelowH, const PetscScalar *zlevels,
                const PetscScalar *enthalpy) const {

  PetscScalar B = 0;
  if (kbelowH > 0) {
    PetscScalar dz = zlevels[1] - zlevels[0];
    B += 0.5 * dz * hardnessParameter_from_enth(enthalpy[0],
                                                EC->getPressureFromDepth(thickness) );
    for (PetscInt m=1; m < kbelowH; m++) {
      const PetscScalar dzNEXT = zlevels[m+1] - zlevels[m],
                        depth  = thickness - 0.5 * (zlevels[m+1] + zlevels[m]);
      B += 0.5 * (dz + dzNEXT) * hardnessParameter_from_enth(enthalpy[m],
                                                             EC->getPressureFromDepth(depth) );
      dz = dzNEXT;
    }
    // use last dz from for loop
    const PetscScalar depth  = 0.5 * (thickness - zlevels[kbelowH]);
    B += 0.5 * dz * hardnessParameter_from_enth(enthalpy[kbelowH],
                                                EC->getPressureFromDepth(depth) );
  }

  // so far B is an integral of ice hardness; compute the average now:
  if (thickness > 0)
    B = B / thickness;
  else
    B = 0;

  return B;
}

PetscTruth IceFlowLawUsesGrainSize(IceFlowLaw *ice) {
  static const PetscReal gs[] = {1e-4,1e-3,1e-2,1},s=1e4,E=500000,p=1e6;
  PetscReal ref = ice->flow_from_enth(s,E,p,gs[0]);
  for (int i=1; i<4; i++) {
    if (ice->flow_from_enth(s,E,p,gs[i]) != ref) return PETSC_TRUE;
  }
  return PETSC_FALSE;
}


/*!
This constructor just sets flow law factor for nonzero water content, from
\ref AschwandenBlatter and \ref LliboutryDuval1985.
 */
GPBLDIce::GPBLDIce(MPI_Comm c,const char pre[],
					 const NCConfigVariable &config) : IceFlowLaw(c,pre,config) {
  T_0  = config.get("water_melting_temperature");    // K
  water_frac_coeff = config.get("gpbld_water_frac_coeff");                
}

PetscErrorCode GPBLDIce::setFromOptions() {
  PetscErrorCode ierr;
  PetscReal slen=schoofLen/1e3,	// convert to km
    svel=schoofVel*secpera;	// convert to m/year

  ierr = PetscOptionsBegin(comm,prefix,"GPBLDIce options",NULL); CHKERRQ(ierr);
  {
    ierr = PetscOptionsReal("-ice_reg_schoof_vel","Regularizing velocity (Schoof definition, m/a)","",svel,&svel,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ice_reg_schoof_length","Regularizing length (Schoof definition, km)","",slen,&slen,NULL);CHKERRQ(ierr);
    schoofVel = svel / secpera;	// convert to m/s
    schoofLen = slen * 1e3;	// convert to meters
    schoofReg = PetscSqr(schoofVel/schoofLen);
    ierr = PetscOptionsReal("-ice_pb_A_cold","Paterson-Budd cold softness parameter (Pa^-3 s^-1)","",A_cold,&A_cold,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ice_pb_A_warm","Paterson-Budd warm softness parameter (Pa^-3 s^-1)","",A_warm,&A_warm,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ice_pb_Q_cold","Paterson-Budd activation energy (J/mol)","",Q_cold,&Q_cold,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ice_pb_Q_warm","Paterson-Budd activation energy (J/mol)","",Q_warm,&Q_warm,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ice_pb_crit_temp","Paterson-Budd critical temperature (K)","",crit_temp,&crit_temp,NULL);CHKERRQ(ierr);

    ierr = PetscOptionsReal("-ice_gpbld_water_frac_coeff",
      "coefficient of softness factor in temperate ice, as function of liquid water fraction (no units)",
      "",water_frac_coeff,&water_frac_coeff,NULL); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode GPBLDIce::view(PetscViewer viewer) const {
  PetscErrorCode ierr;
  PetscTruth iascii;
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(comm,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"ThermoGlenIce object (%s)\n",prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  v_schoof=%4f m/a L_schoof=%4f km\n",schoofVel*secpera,schoofLen/1e3);CHKERRQ(ierr);

    ierr = PetscViewerASCIIPrintf(viewer,"<\nderived GPBLDIce object (%s)\n",prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  T_0             =%10.3f (K)\n",T_0);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  water_frac_coeff=%10.1f\n",water_frac_coeff);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,">\n",water_frac_coeff);CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this object\n");
  }
  return 0;
}

//! Paterson-Budd softness parameter
PetscScalar GPBLDIce::softnessParameter_pb(PetscScalar T_pa) const {
  if (T_pa < crit_temp) {
    return A_cold * exp(-Q_cold/(ideal_gas_constant * T_pa));
  }
  return A_warm * exp(-Q_warm/(ideal_gas_constant * T_pa));
}

//! Paterson-Budd hardness parameter
PetscScalar GPBLDIce::hardnessParameter_pb(PetscScalar T_pa) const {
  return pow(softnessParameter_pb(T_pa), -1.0/n);
}

//! The softness factor in the Glen-Paterson-Budd-Lliboutry-Duval flow law.  For constitutive law form.
/*!
This is a modification of Glen-Paterson-Budd ice, which is ThermoGlenIce.  In particular, if
\f$A()\f$ is the softness factor for ThermoGlenIce, if \f$E\f$ is the enthalpy, and \f$p\f$ is
the pressure then the softness we compute is
   \f[A = A(T_{pa}(E,p))(1+184\omega).\f]
The pressure-melting temperature \f$T_{pa}(E,p)\f$ is computed by getPATemp().
 */
PetscScalar GPBLDIce::softnessParameter_from_enth(
                PetscScalar enthalpy, PetscScalar pressure) const {
  PetscErrorCode ierr;
#ifdef PISM_DEBUG
  if (enthalpy < 0) {
    SETERRQ(1, "Negative enthalpy in GPBLDIce::softnessParameter_from_enth() ... this should never happen");
  }
#endif

  if (EC == NULL) {
    PetscErrorPrintf("EC is NULL in GPBLDIce::flow_from_enth()\n");
    endPrintRank();
  }
  PetscScalar E_s, E_l;
  EC->getEnthalpyInterval(pressure, E_s, E_l);
  if (enthalpy <= E_s) {       // cold ice
    PetscScalar T_pa;
    ierr = EC->getPATemp(enthalpy,pressure,T_pa);
    if (ierr) {
      PetscErrorPrintf(
        "getPATemp() returned ierr>0 in GPBLDIce::softnessParameter_from_enth()\n");
      endPrintRank();
    }
    return softnessParameter_pb( T_pa );
  } else if (enthalpy < E_l) { // temperate ice
    PetscScalar omega;
    ierr = EC->getWaterFraction(enthalpy,pressure,omega);
    if (ierr) {
      PetscErrorPrintf(
        "getWaterFraction() returned ierr>0 in GPBLDIce::softnessParameter_from_enth()\n");
      endPrintRank();
    }
    // next line implements eqn (23) in \ref AschwandenBlatter2009
    return softnessParameter_pb(T_0) * (1.0 + water_frac_coeff * omega);
  } else { // liquid water not allowed
    PetscErrorPrintf("ERROR in PolyThermalGlenPBLDIce::flow_from_temp(): liquid water not allowed\n\n");
    endPrintRank();
    return 0.0;
  }
}

//! The hardness factor in the Paterson-Budd-Lliboutry-Duval flow law.  For viscosity form.
PetscScalar GPBLDIce::hardnessParameter_from_enth(
                PetscScalar enthalpy, PetscScalar pressure) const {
  return pow(softnessParameter_from_enth(enthalpy,pressure), -1.0/n);
}


//! Glen-Paterson-Budd-Lliboutry-Duval flow law itself.
PetscScalar GPBLDIce::flow_from_enth(PetscScalar stress, PetscScalar enthalpy,
                                     PetscScalar pressure, PetscScalar /* gs */) const {
  return softnessParameter_from_enth(enthalpy,pressure) * pow(stress,n-1);
}

//! Returns viscosity and \b not the nu * H product.
PetscScalar GPBLDIce::effectiveViscosity(PetscScalar hardness,
                                         PetscScalar u_x, PetscScalar u_y,
                                         PetscScalar v_x, PetscScalar v_y) const {
  const PetscScalar alpha = secondInvariant(u_x, u_y, v_x, v_y);
  return 0.5 * hardness * pow(schoofReg + alpha, (1-n)/(2*n));
}


// Rather than make this part of the base class, we just check at some reference values.
// PetscTruth IceFlowLawIsPatersonBuddCold(IceFlowLaw *ice, const NCConfigVariable &config) {
//   static const struct {PetscReal s,T,p,gs;} v[] = {
//     {1e3,223,1e6,1e-3},{2e4,254,3e6,2e-3},{5e4,268,5e6,3e-3},{1e5,273,8e6,5e-3}};
//   ThermoGlenArrIce cpb(PETSC_COMM_SELF,NULL,config); // This is unmodified cold Paterson-Budd
//   for (int i=0; i<4; i++) {
//     const PetscReal left  = ice->flow_from_temp(v[i].s, v[i].T, v[i].p, v[i].gs),
//                     right =  cpb.flow_from_temp(v[i].s, v[i].T, v[i].p, v[i].gs);
//     if (PetscAbs((left - right)/left)>1.0e-15) {
//       return PETSC_FALSE;
//     }
//   }
//   return PETSC_TRUE;
// }
