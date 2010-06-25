// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "materials.hh"
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
}


PetscErrorCode IceFlowLaw::printInfo(PetscInt) const {
  PetscPrintf(comm,"WARNING:  IceFlowLaw::printInfo() called but base class IceFlowLaw is virtual!!\n");
  return 0;
}

// Rather than make this part of the base class, we just check at some reference values.
PetscTruth IceFlowLawIsPatersonBuddCold(IceFlowLaw *ice, const NCConfigVariable &config) {
  static const struct {PetscReal s,T,p,gs;} v[] = {
    {1e3,223,1e6,1e-3},{2e4,254,3e6,2e-3},{5e4,268,5e6,3e-3},{1e5,273,8e6,5e-3}};
  ThermoGlenArrIce cpb(PETSC_COMM_SELF,NULL,config); // This is unmodified cold Paterson-Budd
  for (int i=0; i<4; i++) {
    const PetscReal left  = ice->flow(v[i].s, v[i].T, v[i].p, v[i].gs),
                    right =  cpb.flow(v[i].s, v[i].T, v[i].p, v[i].gs);
    if (PetscAbs((left - right)/left)>1.0e-15) {
      return PETSC_FALSE;
    }
  }
  return PETSC_TRUE;
}

PetscTruth IceFlowLawUsesGrainSize(IceFlowLaw *ice) {
  static const PetscReal gs[] = {1e-4,1e-3,1e-2,1},s=1e4,T=260,p=1e6;
  PetscReal ref = ice->flow(s,T,p,gs[0]);
  for (int i=1; i<4; i++) {
    if (ice->flow(s,T,p,gs[i]) != ref) return PETSC_TRUE;
  }
  return PETSC_FALSE;
}




CustomGlenIce::CustomGlenIce(MPI_Comm c,const char pre[], const NCConfigVariable &config) : IceFlowLaw(c,pre,config)
{
  exponent_n = config.get("Glen_exponent");
  softness_A = config.get("ice_softness");
  hardness_B = pow(softness_A, -1/exponent_n);
  setSchoofRegularization(config.get("Schoof_regularizing_velocity"),
			  config.get("Schoof_regularizing_length"));
}


PetscScalar CustomGlenIce::flow(PetscScalar stress,PetscScalar,PetscScalar,PetscScalar) const
{
  return softness_A * pow(stress,exponent_n-1);
}


PetscScalar CustomGlenIce::effectiveViscosityColumn(PetscScalar H, PetscInt, const PetscScalar *,
                           PetscScalar u_x, PetscScalar u_y, PetscScalar v_x, PetscScalar v_y,
                           const PetscScalar *, const PetscScalar *) const  {
  return H * (hardness_B / 2) * pow(schoofReg + secondInvariant(u_x,u_y,v_x,v_y), (1-exponent_n)/(2*exponent_n));
}


PetscErrorCode CustomGlenIce::setDensity(PetscReal r) {rho = r; return 0;}


PetscErrorCode CustomGlenIce::setExponent(PetscReal n) {exponent_n = n; return 0;}


PetscErrorCode CustomGlenIce::setSchoofRegularization(PetscReal vel,PetscReal len) {
  schoofVel = vel/secpera;  // vel has units m/a
  schoofLen = len*1e3;      // len has units km
  schoofReg = PetscSqr(schoofVel/schoofLen); 
  return 0;
}


PetscErrorCode CustomGlenIce::setSoftness(PetscReal A) {
  softness_A = A;
  hardness_B = pow(A,-1/exponent_n); 
  return 0;
}


PetscErrorCode CustomGlenIce::setHardness(PetscReal B) {
  hardness_B = B;
  softness_A = pow(B,-exponent_n);
  return 0;
}


PetscScalar CustomGlenIce::exponent() const { return exponent_n; }


PetscScalar CustomGlenIce::softnessParameter(PetscScalar /*T*/) const { return softness_A; }


PetscScalar CustomGlenIce::hardnessParameter(PetscScalar /*T*/) const { return hardness_B; }


PetscScalar CustomGlenIce::averagedHardness(
                PetscScalar /* H */, PetscInt /* kbelowH */, const PetscScalar /* zlevels */ [],
                const PetscScalar /* T */[]) const  { return hardness_B; }


PetscErrorCode CustomGlenIce::setFromOptions()
{
  PetscReal n = exponent_n,
    B=hardness_B,
    A=softness_A,
    slen=schoofLen/1e3,		// convert to km
    svel=schoofVel*secpera;	// convert to m/year
  PetscTruth     flg;
  PetscErrorCode ierr;

  ierr = PetscOptionsBegin(comm,prefix,"CustomGlenIce options",NULL);CHKERRQ(ierr);
  {
    ierr = PetscOptionsReal("-ice_custom_n","Power-law exponent","setExponent",n,&n,&flg);CHKERRQ(ierr);
    if (flg) {ierr = setExponent(n);CHKERRQ(ierr);}
    ierr = PetscOptionsReal("-ice_custom_softness","Softness parameter A (Pa^{-n} s^{-1})","setSoftness",A,&A,&flg);CHKERRQ(ierr);
    if (flg) {ierr = setSoftness(A);CHKERRQ(ierr);}
    ierr = PetscOptionsReal("-ice_custom_hardness","Hardness parameter B (Pa s^{1/n})","setHardness",B,&B,&flg);CHKERRQ(ierr);
    if (flg) {ierr = setHardness(B);CHKERRQ(ierr);}
    ierr = PetscOptionsReal("-ice_custom_density","Density rho (km m^{-1})","setDensity",rho,&rho,NULL);CHKERRQ(ierr);
    // use -ice_ instead of -ice_custom_ to be compatible with ThermoGlenIce
    ierr = PetscOptionsReal("-ice_reg_schoof_vel","Regularizing velocity (Schoof definition, m/a)","setSchoofRegularization",svel,&svel,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ice_reg_schoof_length","Regularizing length (Schoof definition, km)","setSchoofRegularization",slen,&slen,NULL);CHKERRQ(ierr);
    ierr = setSchoofRegularization(svel,slen);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  return 0;
}


PetscErrorCode CustomGlenIce::printInfo(PetscInt thresh) const {
  PetscErrorCode ierr;
  if (thresh <= getVerbosityLevel()) {
    ierr = view(NULL);CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode CustomGlenIce::view(PetscViewer viewer) const {
  PetscErrorCode ierr;
  PetscTruth iascii;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(comm,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"CustomGlenIce object (%s)\n",prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,
        "  n=%3g,   A=%8.3e,   B=%8.3e,  rho=%.1f,\n"
        "  v_schoof=%.2f m/a,   L_schoof=%.2f km,  schoofReg = %.2e\n",
        exponent_n,softness_A,hardness_B,rho,schoofVel*secpera,schoofLen/1e3,schoofReg);CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this object\n");
  }
  return 0;
}


ThermoGlenIce::ThermoGlenIce(MPI_Comm c,const char pre[], const NCConfigVariable &config) : IceFlowLaw(c,pre,config) {
  n = 3;			// Paterson-Budd has the fixed Glen exponent, so it's hard-wired.
  A_cold = config.get("Paterson-Budd_A_cold");
  A_warm = config.get("Paterson-Budd_A_warm");
  Q_cold = config.get("Paterson-Budd_Q_cold");
  Q_warm = config.get("Paterson-Budd_Q_warm");
  crit_temp = config.get("Paterson-Budd_critical_temperature");
  schoofLen = config.get("Schoof_regularizing_length") * 1e3; // convert to meters
  schoofVel = config.get("Schoof_regularizing_velocity")/secpera; // convert to m/s
  schoofReg = PetscSqr(schoofVel/schoofLen);
}


PetscErrorCode ThermoGlenIce::setFromOptions() {
  PetscErrorCode ierr;
  PetscReal slen=schoofLen/1e3,	// convert to km
    svel=schoofVel*secpera;	// convert to m/year

  ierr = PetscOptionsBegin(comm,prefix,"ThermoGlenIce options",NULL);CHKERRQ(ierr);
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
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  return 0;
}


PetscErrorCode ThermoGlenIce::printInfo(PetscInt thresh) const {
  PetscErrorCode ierr;
  if (thresh <= getVerbosityLevel()) {
    ierr = view(NULL);CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode ThermoGlenIce::view(PetscViewer viewer) const {
  PetscErrorCode ierr;
  PetscTruth iascii;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(comm,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"ThermoGlenIce object (%s)\n",prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  v_schoof=%4f m/a L_schoof=%4f km\n",schoofVel*secpera,schoofLen/1e3);CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this object\n");
  }
  return 0;
}


PetscScalar ThermoGlenIce::flow(PetscScalar stress,PetscScalar temp,PetscScalar pressure,PetscScalar) const {
  const PetscScalar T = temp + (beta_CC_grad / (rho * standard_gravity)) * pressure; // pressure-adjusted temp
  return softnessParameter(T) * pow(stress,n-1);
}


//! DESPITE NAME, does *not* return effective viscosity. The result is \f$\nu_e H\f$,
//! i.e. viscosity times thickness. \f$B\f$ is really hardness times thickness.
PetscScalar ThermoGlenIce::effectiveViscosityColumn(
                PetscScalar H,  PetscInt kbelowH, const PetscScalar *zlevels,
                PetscScalar u_x,  PetscScalar u_y, PetscScalar v_x,  PetscScalar v_y,
                const PetscScalar *T1, const PetscScalar *T2) const {
//  const PetscInt  ks = static_cast<PetscInt>(floor(H/dz));
  // Integrate the hardness parameter using the trapezoid rule.
  PetscScalar B = 0;
  if (kbelowH > 0) {
    PetscScalar dz = zlevels[1] - zlevels[0];
    B += 0.5 * dz * hardnessParameter(0.5 * (T1[0] + T2[0]) + beta_CC_grad * H);
    for (PetscInt m=1; m < kbelowH; m++) {
      const PetscScalar dzNEXT = zlevels[m+1] - zlevels[m];
      B += 0.5 * (dz + dzNEXT) * hardnessParameter(0.5 * (T1[m] + T2[m])
           + beta_CC_grad * (H - zlevels[m]));
      dz = dzNEXT;
    }
    // use last dz
    B += 0.5 * dz * hardnessParameter(0.5 * (T1[kbelowH] + T2[kbelowH])
                                      + beta_CC_grad * (H - zlevels[kbelowH]));
  }
  const PetscScalar alpha = secondInvariant(u_x, u_y, v_x, v_y);
  return 0.5 * B * pow(schoofReg + alpha, (1-n)/(2*n));
}


/*! This is not a natural part of all IceFlowLaw instances since it doesn't make any sense
    for plenty of rheologies.  Specifically, however, we need some exponent to compute the coordinate
    transformation in IceModel::computeDrivingStress (see iMgeometry.cc).  */
PetscScalar ThermoGlenIce::exponent() const { return n; }


//! Return the softness parameter A(T) for a given temperature T.
/*! This is not a natural part of all IceFlowLaw instances.   */
PetscScalar ThermoGlenIce::softnessParameter(PetscScalar T) const {
  if (T < crit_temp) {
    return A_cold * exp(-Q_cold/(ideal_gas_constant * T));
  }
  return A_warm * exp(-Q_warm/(ideal_gas_constant * T));
}


//! Return the hardness parameter B(T) \f$= A(T)^{-1/n}\f$ for a given temperature T.
/*! This is not a natural part of all IceFlowLaw instances, but it is an important optimization 
    in the SSA stress balance to be able to factor the flow law into a temperature-dependent and 
    a stress dependent part, and to do IceModel::correctSigma().   */
PetscScalar ThermoGlenIce::hardnessParameter(PetscScalar T) const {
  return pow(softnessParameter(T), -1.0/n);
}


//! Computes vertical average of B(T) ice hardness, namely \f$\bar B(T)\f$.  See comment for hardnessParameter().
/*! Note T[0],...,T[kbelowH] must be valid.  In particular, even if kbelowH == 0, we still use
    T[0].  Uses trapezoid rule to do integral.  */
PetscScalar ThermoGlenIce::averagedHardness(PetscScalar H, PetscInt kbelowH, const PetscScalar zlevels[],
                                           const PetscScalar T[]) const {
  PetscScalar B;
  if ((kbelowH > 0) && (H > 1.0)) {
    // use trapezoid rule starting from bottom of ice column
    PetscScalar dz = zlevels[1] - zlevels[0];
    B = 0.5 * dz * hardnessParameter(T[0] + beta_CC_grad * H);
    for (PetscInt m=1; m < kbelowH; m++) {
      dz = zlevels[m+1] - zlevels[m-1];
      B += 0.5 * dz * hardnessParameter(T[m] + beta_CC_grad * (H - zlevels[m]));
    }
    // in top piece, between last zlevel in ice and ice surface, use T at lower end of subinterval
    dz = H - zlevels[kbelowH];
    B += 0.5 * dz * hardnessParameter(T[kbelowH] + beta_CC_grad * (H - zlevels[kbelowH]));
    // B now contains integral, but we want average
    return B / H;
  } else {
    // if no significant ice thickness, use T[0] directly
    return hardnessParameter(T[0]);
  }
}


/*!
This constructor just sets flow law factor for nonzero water content, from
\ref AschwandenBlatter and \ref LliboutryDuval1985.
 */
PolyThermalGPBLDIce::PolyThermalGPBLDIce(MPI_Comm c,const char pre[],
					 const NCConfigVariable &config) : ThermoGlenIce(c,pre,config) {
  EC = new EnthalpyConverter(config);
  T_0  = config.get("water_melting_temperature");    // K
  water_frac_coeff = config.get("gpbld_water_frac_coeff");                
}

PolyThermalGPBLDIce::~PolyThermalGPBLDIce() {
  delete EC;
}


PetscErrorCode PolyThermalGPBLDIce::setFromOptions() {
  PetscErrorCode ierr;

  ierr = ThermoGlenIce::setFromOptions(); CHKERRQ(ierr);
  
  ierr = PetscOptionsBegin(comm,prefix,"PolyThermalGPBLDIce options",NULL); CHKERRQ(ierr);
  {
    ierr = PetscOptionsReal("-ice_gpbld_water_frac_coeff",
      "coefficient of softness factor in temperate ice, as function of liquid water fraction (no units)",
      "",water_frac_coeff,&water_frac_coeff,NULL); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PolyThermalGPBLDIce::view(PetscViewer viewer) const {
  PetscErrorCode ierr;

  ierr = ThermoGlenIce::view(viewer); CHKERRQ(ierr);
  
  PetscTruth iascii;
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(comm,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"<\nderived PolyThermalGPBLDIce object (%s)\n",prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  T_0             =%10.3f (K)\n",T_0);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  water_frac_coeff=%10.1f\n",water_frac_coeff);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,">\n",water_frac_coeff);CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this object\n");
  }
  return 0;
}


//! The softness factor in the Glen-Paterson-Budd-Lliboutry-Duval flow law.  For constitutive law form.
/*!
This is a modification of Glen-Paterson-Budd ice, which is ThermoGlenIce.  In particular, if
\f$A()\f$ is the softness factor for ThermoGlenIce, if \f$E\f$ is the enthalpy, and \f$p\f$ is
the pressure then the softness we compute is
   \f[A = A(T_{pa}(E,p))(1+184\omega).\f]
The pressure-melting temperature \f$T_{pa}(E,p)\f$ is computed by getPATemp().
 */
PetscScalar PolyThermalGPBLDIce::softnessParameterFromEnth(
                PetscScalar enthalpy, PetscScalar pressure) const {
  PetscErrorCode ierr;
#ifdef PISM_DEBUG
  if (enthalpy < 0) {
    SETERRQ(1, "Negative enthalpy in PolyThermalGPBLDIce::softnessParameterFromEnth() ... this should never happen");
  }
#endif

  if (EC == NULL) {
    PetscErrorPrintf("EC is NULL in PolyThermalGPBLDIce::flowFromEnth()\n");
    endPrintRank();
  }
  PetscScalar E_s, E_l;
  EC->getEnthalpyInterval(pressure, E_s, E_l);
  if (enthalpy <= E_s) {       // cold ice
    PetscScalar T_pa;
    ierr = EC->getPATemp(enthalpy,pressure,T_pa);
    if (ierr) {
      PetscErrorPrintf(
        "getPATemp() returned ierr>0 in PolyThermalGPBLDIce::softnessParameterFromEnth()\n");
      endPrintRank();
    }
    return softnessParameter( T_pa ); // uses ThermoGlenIce formula
  } else if (enthalpy < E_l) { // temperate ice
    PetscScalar omega;
    ierr = EC->getWaterFraction(enthalpy,pressure,omega);
    if (ierr) {
      PetscErrorPrintf(
        "getWaterFraction() returned ierr>0 in PolyThermalGPBLDIce::softnessParameterFromEnth()\n");
      endPrintRank();
    }
    // next line implements eqn (23) in \ref AschwandenBlatter2009
    return softnessParameter(T_0) * (1.0 + water_frac_coeff * omega);  // uses ThermoGlenIce formula
  } else { // liquid water not allowed
    PetscErrorPrintf("ERROR in PolyThermalGlenPBLDIce::flow(): liquid water not allowed\n\n");
    endPrintRank();
    return 0.0;
  }
}


//! The hardness factor in the Paterson-Budd-Lliboutry-Duval flow law.  For viscosity form.
PetscScalar PolyThermalGPBLDIce::hardnessParameterFromEnth(
                PetscScalar enthalpy, PetscScalar pressure) const {
  return pow(softnessParameterFromEnth(enthalpy,pressure), -1.0/n);
}


//! Glen-Paterson-Budd-Lliboutry-Duval flow law itself.
PetscScalar PolyThermalGPBLDIce::flowFromEnth(
                PetscScalar stress, PetscScalar enthalpy, PetscScalar pressure, PetscScalar /* gs */) const {
  return softnessParameterFromEnth(enthalpy,pressure) * pow(stress,n-1);
}


PetscScalar PolyThermalGPBLDIce::effectiveViscosityColumnFromEnth(
                PetscScalar thickness,  PetscInt kbelowH, const PetscScalar *zlevels,
                PetscScalar u_x,  PetscScalar u_y, PetscScalar v_x,  PetscScalar v_y,
                const PetscScalar *enthalpy1, const PetscScalar *enthalpy2) const {
  if (EC == NULL) {
    PetscErrorPrintf(
      "EC is NULL in PolyThermalGPBLDIce::effectiveViscosityColumnFromEnth()\n");
    endPrintRank();
  }

  // result is \nu_e H, i.e. viscosity times thickness; B is really hardness times thickness
  // integrates the hardness parameter using the trapezoid rule.
  PetscScalar B = 0;
  if (kbelowH > 0) {
    PetscScalar dz = zlevels[1] - zlevels[0];
    B += 0.5 * dz * hardnessParameterFromEnth( 0.5 * (enthalpy1[0] + enthalpy2[0]),
                                               EC->getPressureFromDepth(thickness) );
    for (PetscInt m=1; m < kbelowH; m++) {
      const PetscScalar dzNEXT = zlevels[m+1] - zlevels[m],
                        depth  = thickness - 0.5 * (zlevels[m+1] + zlevels[m]);
      B += 0.5 * (dz + dzNEXT) * hardnessParameterFromEnth( 0.5 * (enthalpy1[m] + enthalpy2[m]),
                                                            EC->getPressureFromDepth(depth) );
      dz = dzNEXT;
    }
    // use last dz from for loop
    const PetscScalar depth  = 0.5 * (thickness - zlevels[kbelowH]);
    B += 0.5 * dz * hardnessParameterFromEnth( 0.5 * (enthalpy1[kbelowH] + enthalpy2[kbelowH]),
                                               EC->getPressureFromDepth(depth) );
  }
  const PetscScalar alpha = secondInvariant(u_x, u_y, v_x, v_y);
  return 0.5 * B * pow(schoofReg + alpha, (1-n)/(2*n));
}


//! This is the Hooke flow law, see [\ref Hooke].
ThermoGlenIceHooke::ThermoGlenIceHooke(MPI_Comm c,const char pre[],
				       const NCConfigVariable &config) : ThermoGlenIce(c,pre,config) {
  Q_Hooke  = config.get("Hooke_Q");
  A_Hooke  = config.get("Hooke_A");
  C_Hooke  = config.get("Hooke_C");
  K_Hooke  = config.get("Hooke_k");
  Tr_Hooke = config.get("Hooke_Tr");
  R_Hooke  = config.get("ideal_gas_constant");
}

//! This method implements formula (7) in [\ref Hooke].
PetscScalar ThermoGlenIceHooke::softnessParameter(PetscScalar T) const {

  return A_Hooke * exp( -Q_Hooke/(R_Hooke * T)
                       + 3.0 * C_Hooke * pow(Tr_Hooke - T,-K_Hooke));
}


PetscErrorCode ThermoGlenArrIce::view(PetscViewer viewer) const {
  PetscErrorCode ierr;
  PetscTruth iascii;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(comm,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"ThermoGlenArrIce object (%s)\n",prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  No customizable options specific to ThermoGlenArrIce\n");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  Derived from ThermoGlenIce with the following state\n");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = ThermoGlenIce::view(viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this object\n");
  }
  return 0;
}


//! Return the softness parameter A(T) for a given temperature T.
PetscScalar ThermoGlenArrIce::softnessParameter(PetscScalar T) const {
  return A() * exp(-Q()/(ideal_gas_constant * T));
}


//! Return the temperature T corresponding to a given value A=A(T).
PetscScalar ThermoGlenArrIce::tempFromSoftness(PetscScalar myA) const {
  return - Q() / (ideal_gas_constant * (log(myA) - log(A())));
}


PetscScalar ThermoGlenArrIce::flow(PetscScalar stress, PetscScalar temp, PetscScalar,PetscScalar) const {
  // ignores pressure
  return softnessParameter(temp) * pow(stress,n-1);  // uses NON-pressure-adjusted temp
}

PetscScalar ThermoGlenArrIce::A() const {
  return A_cold;
}

PetscScalar ThermoGlenArrIce::Q() const {
  return Q_cold;
}


PetscErrorCode ThermoGlenArrIceWarm::view(PetscViewer viewer) const {
  PetscErrorCode ierr;
  PetscTruth iascii;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(comm,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"ThermoGlenArrIceWarm object (%s)\n",prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  No customizable options specific to ThermoGlenArrIceWarm\n");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  Derived from ThermoGlenArrIce with the following state\n");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = ThermoGlenArrIce::view(viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this object\n");
  }
  return 0;
}

PetscScalar ThermoGlenArrIceWarm::A() const {
  return A_warm;
}

PetscScalar ThermoGlenArrIceWarm::Q() const {
  return Q_warm;
}


HybridIce::HybridIce(MPI_Comm c,const char pre[],
		     const NCConfigVariable &config) : ThermoGlenIce(c,pre,config) {

  V_act_vol    = -13.e-6;  // m^3/mol
  d_grain_size = 1.0e-3;   // m  (see p. ?? of G&K paper)
  //--- dislocation creep ---
  disl_crit_temp=258.0,    // Kelvin
  //disl_A_cold=4.0e5,                  // MPa^{-4.0} s^{-1}
  //disl_A_warm=6.0e28,                 // MPa^{-4.0} s^{-1}
  disl_A_cold=4.0e-19,     // Pa^{-4.0} s^{-1}
  disl_A_warm=6.0e4,       // Pa^{-4.0} s^{-1} (GK)
  disl_n=4.0,              // stress exponent
  disl_Q_cold=60.e3,       // J/mol Activation energy
  disl_Q_warm=180.e3;      // J/mol Activation energy (GK)
  //--- grain boundary sliding ---
  gbs_crit_temp=255.0,     // Kelvin
  //gbs_A_cold=3.9e-3,                  // MPa^{-1.8} m^{1.4} s^{-1}
  //gbs_A_warm=3.e26,                   // MPa^{-1.8} m^{1.4} s^{-1}
  gbs_A_cold=6.1811e-14,   // Pa^{-1.8} m^{1.4} s^{-1}
  gbs_A_warm=4.7547e15,    // Pa^{-1.8} m^{1.4} s^{-1}
  gbs_n=1.8,               // stress exponent
  gbs_Q_cold=49.e3,        // J/mol Activation energy
  gbs_Q_warm=192.e3,       // J/mol Activation energy
  p_grain_sz_exp=1.4;      // from Peltier
  //--- easy slip (basal) ---
  //basal_A=5.5e7,                      // MPa^{-2.4} s^{-1}
  basal_A=2.1896e-7,       // Pa^{-2.4} s^{-1}
  basal_n=2.4,             // stress exponent
  basal_Q=60.e3;           // J/mol Activation energy
  //--- diffusional flow ---
  diff_crit_temp=258.0,    // when to use enhancement factor
  diff_V_m=1.97e-5,        // Molar volume (m^3/mol)
  diff_D_0v=9.10e-4,       // Preexponential volume diffusion (m^2/s)
  diff_Q_v=59.4e3,         // activation energy, vol. diff. (J/mol)
  diff_D_0b=5.8e-4,        // preexponential grain boundary coeff.
  diff_Q_b=49.e3,          // activation energy, g.b. (J/mol)
  diff_delta=9.04e-10;     // grain boundary width (m)
}


/*!
  This is the (forward) Goldsby-Kohlstedt flow law.  See:
  D. L. Goldsby & D. L. Kohlstedt (2001), "Superplastic deformation
  of ice: experimental observations", J. Geophys. Res. 106(M6), 11017-11030.
*/
PetscScalar HybridIce::flow(PetscScalar stress, PetscScalar temp,
                            PetscScalar pressure, PetscScalar gs) const {
  PetscScalar eps_diff, eps_disl, eps_basal, eps_gbs, diff_D_b;

  if (PetscAbs(stress) < 1e-10) return 0;
  const PetscScalar T = temp + (beta_CC_grad / (rho * standard_gravity)) * pressure;
  const PetscScalar pV = pressure * V_act_vol;
  const PetscScalar RT = ideal_gas_constant * T;
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
GKparts HybridIce::flowParts(PetscScalar stress,PetscScalar temp,PetscScalar pressure) const {
  PetscScalar gs, eps_diff, eps_disl, eps_basal, eps_gbs, diff_D_b;
  GKparts p;

  if (PetscAbs(stress) < 1e-10) {
    p.eps_total=0.0;
    p.eps_diff=0.0; p.eps_disl=0.0; p.eps_gbs=0.0; p.eps_basal=0.0;
    return p;
  }
  const PetscScalar T = temp + (beta_CC_grad / (rho * standard_gravity)) * pressure;
  const PetscScalar pV = pressure * V_act_vol;
  const PetscScalar RT = ideal_gas_constant * T;
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

PetscErrorCode HybridIce::printInfo(PetscInt thresh) const {
  PetscErrorCode ierr;
  if (thresh <= getVerbosityLevel()) {
    ierr = view(NULL);CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode HybridIce::view(PetscViewer viewer) const {
  PetscErrorCode ierr;
  PetscTruth iascii;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(comm,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"HybridIce object (%s)\n",prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  No customizable options specific to HybridIce\n");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  Derived from ThermoGlenIce with the following state\n");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = ThermoGlenIce::view(viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this object\n");
  }
  return 0;
}


HybridIceStripped::HybridIceStripped(MPI_Comm c,const char pre[],
				     const NCConfigVariable &config) : HybridIce(c,pre,config) {
  d_grain_size_stripped = 3.0e-3;  // m; = 3mm  (see Peltier et al 2000 paper)
}


PetscScalar HybridIceStripped::flow(PetscScalar stress, PetscScalar temp, PetscScalar pressure, PetscScalar) const {
  // note value of gs is ignored
  // note pressure only effects the temperature; the "P V" term is dropped
  // note no diffusional flow
  PetscScalar eps_disl, eps_basal, eps_gbs;

  if (PetscAbs(stress) < 1e-10) return 0;
  const PetscScalar T = temp + (beta_CC_grad / (rho * standard_gravity)) * pressure;
  const PetscScalar RT = ideal_gas_constant * T;
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


IceBasalResistancePlasticLaw::IceBasalResistancePlasticLaw(
             PetscScalar regularizationConstant, bool pseudoPlastic,
             PetscScalar pseudoExponent, PetscScalar pseudoUThreshold) {
  plastic_regularize = regularizationConstant;
  pseudo_plastic = pseudoPlastic;
  pseudo_q = pseudoExponent;
  pseudo_u_threshold = pseudoUThreshold;
}


PetscErrorCode IceBasalResistancePlasticLaw::printInfo(int verbthresh, MPI_Comm com) {
  PetscErrorCode ierr;
  if (pseudo_plastic == PETSC_TRUE) {
    if (pseudo_q == 1.0) {
      ierr = verbPrintf(verbthresh, com, 
        "Using linearly viscous till with u_threshold = %.2f m/a.\n", 
        pseudo_u_threshold * secpera); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(verbthresh, com, 
        "Using pseudo-plastic till with eps = %10.5e m/a, q = %.4f,"
        " and u_threshold = %.2f m/a.\n", 
        plastic_regularize * secpera, pseudo_q, pseudo_u_threshold * secpera); 
        CHKERRQ(ierr);
    }
  } else {
    ierr = verbPrintf(verbthresh, com, 
      "Using purely plastic till with eps = %10.5e m/a.\n",
      plastic_regularize * secpera); CHKERRQ(ierr);
  }
  return 0;
}


//! Compute the drag coefficient for the basal shear stress.
/*!
The basal shear stress term \f$\tau_b\f$ in the SSA stress balance for ice
is minus the return value here times (vx,vy).

Purely plastic is the pseudo_q = 0.0 case; linear is pseudo_q = 1.0; set 
pseudo_q using IceBasalResistancePlasticLaw constructor.
 */
PetscScalar IceBasalResistancePlasticLaw::drag(PetscScalar tauc,
                                   PetscScalar vx, PetscScalar vy) {
  const PetscScalar magreg2 = PetscSqr(plastic_regularize) + PetscSqr(vx) + PetscSqr(vy);
  if (pseudo_plastic == PETSC_TRUE) {
    return tauc * pow(magreg2, 0.5*(pseudo_q - 1)) * pow(pseudo_u_threshold, -pseudo_q);
  } else { // pure plastic, but regularized
    return tauc / sqrt(magreg2);
  }
}

// Derivative of drag with respect to \f$ \alpha = \frac 1 2 (u_x^2 + u_y^2) \f$
void IceBasalResistancePlasticLaw::dragWithDerivative(PetscReal tauc, PetscScalar vx, PetscScalar vy,
						      PetscScalar *d, PetscScalar *dd) const
{
  const PetscScalar magreg2 = PetscSqr(plastic_regularize) + PetscSqr(vx) + PetscSqr(vy);
  if (pseudo_plastic == PETSC_TRUE) {
    *d = tauc * pow(magreg2, 0.5*(pseudo_q - 1)) * pow(pseudo_u_threshold, -pseudo_q);
    if (dd) *dd = (pseudo_q - 1) * *d / magreg2;
  } else { // pure plastic, but regularized
    *d = tauc / sqrt(magreg2);
    if (dd) *dd = -1 * *d / magreg2;
  }
}


SSAStrengthExtension::SSAStrengthExtension() {
  min_thickness = 50.0;   // m
          // minimum thickness (for SSA velocity computation) at which 
          // NuH switches from vertical integral to constant value
          // this value strongly related to calving front
          // force balance, but the geometry itself is not affected by this value
  const PetscReal
    DEFAULT_CONSTANT_HARDNESS_FOR_SSA = 1.9e8,  // Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
    DEFAULT_TYPICAL_STRAIN_RATE = (100.0 / secpera) / (100.0 * 1.0e3);  // typical strain rate is 100 m/yr per 
  nuH = min_thickness * DEFAULT_CONSTANT_HARDNESS_FOR_SSA
                       / (2.0 * pow(DEFAULT_TYPICAL_STRAIN_RATE,2./3.)); // Pa s m
          // COMPARE: 30.0 * 1e6 * secpera = 9.45e14 is Ritz et al (2001) value of
          //          30 MPa yr for \bar\nu
}

SSAStrengthExtension::~SSAStrengthExtension() {
  // do nothing
}

PetscErrorCode SSAStrengthExtension::set_notional_strength(PetscReal my_nuH) {
  nuH = my_nuH;
  return 0;
}

PetscErrorCode SSAStrengthExtension::set_min_thickness(PetscReal my_min_thickness) {
  min_thickness = my_min_thickness;
  return 0;
}

PetscReal SSAStrengthExtension::notional_strength() const {
  return nuH;
}

PetscReal SSAStrengthExtension::min_thickness_for_extension() const {
  return min_thickness;
}


#undef ALEN
#define ALEN(a) (sizeof(a)/sizeof(a)[0])

IceFlowLawFactory::IceFlowLawFactory(MPI_Comm c,const char pre[], const NCConfigVariable &conf) : config(conf)
{
  comm = c;
  prefix[0] = 0;
  if (pre) {
    PetscStrncpy(prefix,pre,sizeof(prefix));
  }
  if (registerAll()) {
    PetscPrintf(comm,"IceFlowLawFactory::registerAll returned an error but we're in a constructor\n");
    PetscEnd();
  }
  if (setType(ICE_PB)) {       // Set's a default type
    PetscPrintf(comm,"IceFlowLawFactory::setType(\"%s\") returned an error, but we're in a constructor\n",ICE_PB);
    PetscEnd();
  }
}

IceFlowLawFactory::~IceFlowLawFactory()
{
  PetscFListDestroy(&type_list);
}

#undef __FUNCT__
#define __FUNCT__ "IceFlowLawFactory::registerType"
PetscErrorCode IceFlowLawFactory::registerType(const char tname[],
					PetscErrorCode(*icreate)(MPI_Comm,const char[],const NCConfigVariable &, IceFlowLaw**))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFListAdd(&type_list,tname,NULL,(void(*)(void))icreate);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


static PetscErrorCode create_custom(MPI_Comm comm,const char pre[], const NCConfigVariable &config, IceFlowLaw **i) {
  *i = new (CustomGlenIce)(comm, pre, config);  return 0;
}
static PetscErrorCode create_pb(MPI_Comm comm,const char pre[], const NCConfigVariable &config, IceFlowLaw **i) {
  *i = new (ThermoGlenIce)(comm, pre, config);  return 0;
}
static PetscErrorCode create_gpbld(MPI_Comm comm,const char pre[], const NCConfigVariable &config, IceFlowLaw **i) {
  *i = new (PolyThermalGPBLDIce)(comm, pre, config);  return 0;
}
static PetscErrorCode create_hooke(MPI_Comm comm,const char pre[], const NCConfigVariable &config, IceFlowLaw **i) {
  *i = new (ThermoGlenIceHooke)(comm, pre, config);  return 0;
}
static PetscErrorCode create_arr(MPI_Comm comm,const char pre[], const NCConfigVariable &config, IceFlowLaw **i) {
  *i = new (ThermoGlenArrIce)(comm, pre, config);  return 0;
}
static PetscErrorCode create_arrwarm(MPI_Comm comm,const char pre[], const NCConfigVariable &config, IceFlowLaw **i) {
  *i = new (ThermoGlenArrIceWarm)(comm, pre, config);  return 0;
}
static PetscErrorCode create_hybrid(MPI_Comm comm,const char pre[], const NCConfigVariable &config, IceFlowLaw **i) {
  *i = new (HybridIce)(comm, pre, config);  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "IceFlowLawFactory::registerAll"
PetscErrorCode IceFlowLawFactory::registerAll()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMemzero(&type_list,sizeof(type_list));CHKERRQ(ierr);
  ierr = registerType(ICE_CUSTOM, &create_custom); CHKERRQ(ierr);
  ierr = registerType(ICE_PB,     &create_pb);     CHKERRQ(ierr);
  ierr = registerType(ICE_GPBLD,  &create_gpbld);  CHKERRQ(ierr);
  ierr = registerType(ICE_HOOKE,  &create_hooke);  CHKERRQ(ierr);
  ierr = registerType(ICE_ARR,    &create_arr);    CHKERRQ(ierr);
  ierr = registerType(ICE_ARRWARM,&create_arrwarm);CHKERRQ(ierr);
  ierr = registerType(ICE_HYBRID, &create_hybrid); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "IceFlowLawFactory::setType"
PetscErrorCode IceFlowLawFactory::setType(const char type[])
{
  void (*r)(void);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFListFind(type_list,comm,type,(void(**)(void))&r);CHKERRQ(ierr);
  if (!r) {
    ierr = PetscPrintf(comm, "PISM ERROR: Selected ice type \"%s\" is not available.\n",type); CHKERRQ(ierr);
    PetscEnd();
  }
  ierr = PetscStrncpy(type_name,type,sizeof(type_name));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IceFlowLawFactory::setFromOptions"
PetscErrorCode IceFlowLawFactory::setFromOptions()
{
  PetscErrorCode ierr;
  PetscTruth flg;
  char my_type_name[256];

  PetscFunctionBegin;
  // These options will choose Goldsby-Kohlstedt ice by default (see IceModel::setFromOptions()) but if a derived class
  // uses a different initialization procedure, we'll recognize them here as well.  A better long-term solution would be
  // to separate tracking of grain size from a particular flow law (since in principle they are unrelated) but since
  // HYBRID is the only one that currently uses grain size, this solution is acceptable.
  ierr = PetscOptionsHasName(prefix, "-gk_age", &flg); CHKERRQ(ierr);
  if (flg) {
    ierr = setType(ICE_HYBRID);CHKERRQ(ierr);
  }
  // -gk 0 does not make sense, so using PetscOptionsHasName is OK.
  ierr = PetscOptionsHasName(prefix, "-gk", &flg); CHKERRQ(ierr);
  if (flg) {
    ierr = setType(ICE_HYBRID);CHKERRQ(ierr);
  }
  ierr = PetscOptionsBegin(comm,prefix,"IceFlowLawFactory options","IceFlowLaw");CHKERRQ(ierr);
  {
    ierr = PetscOptionsList("-ice_type","Ice type","IceFlowLawFactory::setType",
                            type_list,type_name,my_type_name,sizeof(my_type_name),&flg);CHKERRQ(ierr);
    if (flg) {ierr = setType(my_type_name);CHKERRQ(ierr);}
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  
//  ierr = PetscPrintf(comm,"IceFlowLawFactory::type_name=%s at end of IceFlowLawFactory::setFromOptions()\n",
//                     type_name); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "IceFlowLawFactory::create"
PetscErrorCode IceFlowLawFactory::create(IceFlowLaw **inice)
{
  PetscErrorCode ierr,(*r)(MPI_Comm,const char[],const NCConfigVariable &,IceFlowLaw**);
  IceFlowLaw *ice;

  PetscFunctionBegin;
  PetscValidPointer(inice,3);
  *inice = 0;
  // find the function that can create selected ice type:
  ierr = PetscFListFind(type_list,comm,type_name,(void(**)(void))&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(1,"Selected Ice type %s not available, but we shouldn't be able to get here anyway",type_name);
  // create an IceFlowLaw instance:
  ierr = (*r)(comm,prefix,config,&ice);CHKERRQ(ierr);
  *inice = ice;
  PetscFunctionReturn(0);
}

