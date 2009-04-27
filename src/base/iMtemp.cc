// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petscda.h>
#include "iceModelVec.hh"
#include "columnSystem.hh"
#include "iceModel.hh"


//! Manage the solution of the temperature and age equations, and parallel communication.
PetscErrorCode IceModel::temperatureAgeStep() {
  PetscErrorCode  ierr;

  PetscScalar  myCFLviolcount = 0.0,   // these are counts but they are type "PetscScalar"
               myVertSacrCount = 0.0,  // because that type works with PetscGlobalSum()
               mybulgeCount = 0.0;
  PetscScalar gVertSacrCount, gbulgeCount;

  // values go in vtaunew; note: 3D CFL check happens here:
  ierr = ageStep(&myCFLviolcount); CHKERRQ(ierr);
    
  // new temperature values go in vTnew; also updates Hmelt:
  ierr = temperatureStep(&myVertSacrCount,&mybulgeCount); CHKERRQ(ierr);  

  // start temperature & age communication
  ierr = T3.beginGhostCommTransfer(Tnew3); CHKERRQ(ierr);
  ierr = tau3.beginGhostCommTransfer(taunew3); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&myCFLviolcount, &CFLviolcount, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&myVertSacrCount, &gVertSacrCount, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&mybulgeCount, &gbulgeCount, grid.com); CHKERRQ(ierr);

  if (gVertSacrCount > 0.0) { // count of when BOMBPROOF switches to lower accuracy
    const PetscScalar bfsacrPRCNT = 100.0 * (gVertSacrCount / (grid.Mx * grid.My));
    const PetscScalar BPSACR_REPORT_VERB2_PERCENT = 5.0; // only report if above 5%
    if (bfsacrPRCNT > BPSACR_REPORT_VERB2_PERCENT) {
      ierr = verbPrintf(2,grid.com," [BPsacr=%.4f%%] ", bfsacrPRCNT); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(3,grid.com," [BPsacr=%.4f%%] ", bfsacrPRCNT); CHKERRQ(ierr);
    }
  }

  if (gbulgeCount > 0.0) {   // count of when advection bulges are limited;
                             //    frequently it is identically zero
    const PetscScalar bulgePRCNT
             = 100.0 * (gbulgeCount / (grid.Mx * grid.My * (grid.Mz + grid.Mbz)) );
    const PetscScalar BPBULGE_REPORT_PERCENT = 0.0001; // only report if above 1 / 10^6
    if (bulgePRCNT > BPBULGE_REPORT_PERCENT) {
      ierr = verbPrintf(2,grid.com," [BPbulge=%.5f%%] ", bulgePRCNT); CHKERRQ(ierr);
    }
  }

  // complete temperature & age communication
  ierr = T3.endGhostCommTransfer(Tnew3); CHKERRQ(ierr);
  ierr = tau3.endGhostCommTransfer(taunew3); CHKERRQ(ierr);
  return 0;
}


//! Takes a semi-implicit time-step for the temperature equation.
/*!
In summary, the conservation of energy equation is
    \f[ \rho c_p(T) \frac{dT}{dt} = k \frac{\partial^2 T}{\partial z^2} + \Sigma,\f] 
where \f$T(t,x,y,z)\f$ is the temperature of the ice.  This equation is the shallow approximation
of the full 3D conservation of energy.  Note \f$dT/dt\f$ stands for the material
derivative, so advection is included.  Here \f$\rho\f$ is the density of ice, 
\f$c_p\f$ is its specific heat, and \f$k\f$ is its conductivity.  Also \f$\Sigma\f$ is the volume
strain heating (with SI units of J$/(\text{s} \text{m}^3) = \text{W}\,\text{m}^{-3}$).

\latexonly\index{BOMBPROOF!implementation for temperature equation}\endlatexonly
Note that both the temperature equation and the age equation involve advection.
We handle the horizontal advection explicitly by first-order upwinding.  We handle the
vertical advection implicitly by centered differencing when possible, and a retreat to
implicit first-order upwinding when necessary.  There is a CFL condition
for the horizontal explicit upwinding \lo\cite{MortonMayers}\elo.  We report any CFL violations,
but they are designed to not occur.

The vertical conduction term is handled implicitly (i.e. by backward Euler).

We work from the bottom of the column upward in building the system to solve
(in the semi-implicit time-stepping scheme).  The excess energy above pressure melting
is converted to melt-water, and that a fraction of this melt water is transported to 
the ice base according to the scheme in excessToFromBasalMeltLayer().

The method uses equally-spaced calculation but the methods getValColumn(), 
setValColumn() interpolate back-and-forth from this equally-spaced calculational
grid to the (usually) non-equally space storage grid.

In this procedure four scalar fields are modified: vHmelt, vbasalMeltRate, Tb3, and Tnew3.
But vHmelt, vbasalMeltRate and Tb3 will never need to communicate ghosted values (i.e. horizontal 
stencil neighbors.  The ghosted values for T3 are updated from the values in Tnew3 in the
communication done by temperatureAgeStep().
 
Here is a more complete discussion and derivation.

Consider a column of a slowly flowing and heat conducting material as shown in the left side 
of the next figure.  (The left side shows a general column of flowing and heat conduction 
material showing a small segment \f$V\f$.  The right side shows a more specific column of ice 
flowing and sliding over bedrock.)  This is an \em Eulerian view so the material flows through 
a column which remains fixed (and is notional).  The column is vertical.  We will 
assume when needed that it is rectangular in cross-section with cross-sectional area 
\f$\Delta x\Delta y\f$.

\image latex earlycols.png "Left: a general column of material.  Right: ice over bedrock." width=3in

[FIXME: CONTINUE TO MINE eqns3D.tex FOR MORE]

The application of the geothermal flux at the base of a column is a special case for which 
we give a finite difference argument.  This scheme follows the equation (2.114) in 
\lo\cite{MortonMayers}\elo.  We have the boundary condition
	\f[  -k \frac{\partial T}{\partial z} = G(t,x,y) \f]
where \f$G(t,x,y)\f$ is the applied geothermal flux, and it is applied at level \f$z=-B_0\f$ 
in the bedrock (which is the only case considered here).  We <em> add a virtual lower grid 
point </em> \f$z_{-1} = z_0 - \Delta  z\f$ and we approximate the above boundary condition 
at \f$z_0\f$ by the centered-difference
	\f[  -k \frac{T_{1} - T_{-1}}{2 \Delta z} = G. \f]
Here \f$T_k = T_{ijk}^{l+1}\f$.  We also apply the discretized conduction equation 
at \f$z_0\f$.  These two combined equations yield a simplified form
	\f[(1 + 2 KR) T_0 - 2 K T_1 = T_0 + \frac{2\Delta t}{\rho c_p \Delta z} G \f]
where \f$K = k \Delta t (\rho c \Delta z^2)^{-1}\f$.

The scheme BOMBPROOF is very reliable, but there is still an extreme and rare fjord
situation.  For example, it occurs in one column of ice in one fjord perhaps only once
in a 200ka simulation of the whole sheet, in my (ELB) experience modeling the Greenland 
ice sheet.  It causes the discretized advection bulge to give temperatures below that 
of the coldest ice anywhere, a continuum impossibility.  So as a final protection
there is a "bulge limiter" which sets the temperature to the surface temperature of
the column minus the bulge maximum (15 K) if it is below that level.  The number of 
times this occurs is reported as a "BPbulge" percentage.
 */
PetscErrorCode IceModel::temperatureStep(
     PetscScalar* vertSacrCount, PetscScalar* bulgeCount) {
  PetscErrorCode  ierr;

  PetscInt    Mz, Mbz;
  PetscScalar dzEQ, dzbEQ, *zlevEQ, *zblevEQ;

  ierr = getMzMbzForTempAge(Mz, Mbz); CHKERRQ(ierr);

  zlevEQ = new PetscScalar[Mz];
  zblevEQ = new PetscScalar[Mbz];

  ierr = getVertLevsForTempAge(Mz, Mbz, dzEQ, dzbEQ, zlevEQ, zblevEQ); CHKERRQ(ierr);

  ierr = verbPrintf(5,grid.com,
    "\n  [entering temperatureStep(); Mz = %d, dzEQ = %5.3f, Mbz = %d, dzbEQ = %5.3f]",
    Mz, dzEQ, Mbz, dzbEQ); CHKERRQ(ierr);

  tempSystemCtx system(Mz,Mbz);
  
  ierr = system.tempSetConstants(
           grid.dx, grid.dy, dtTempAge, dzEQ, dzbEQ,
           ice->rho, ice->c_p, ice->k, 
           bed_thermal.rho, bed_thermal.c_p, bed_thermal.k); CHKERRQ(ierr);

  const PetscInt k0 = Mbz - 1;
  PetscScalar *x;  
  x = new PetscScalar[Mz + k0]; // space for solution of system; length = Mz + Mbz - 1 

  // constants needed after solution of system, in insertion phase
  const PetscScalar rho_c_I = ice->rho * ice->c_p,
                    rho_c_br = bed_thermal.rho * bed_thermal.c_p,
                    rho_c_av = (dzEQ * rho_c_I + dzbEQ * rho_c_br) / (dzEQ + dzbEQ);

  PetscScalar *Tb, *Tbnew;
  PetscScalar **Ts, **Tshelfbase, **H, **Ghf, **mask, **Hmelt, **Rb, **basalMeltRate, **bmr_float;

  PetscScalar *u, *v, *w, *Sigma, *T, *Tnew;
  u = new PetscScalar[Mz];
  v = new PetscScalar[Mz];
  w = new PetscScalar[Mz];
  Sigma = new PetscScalar[Mz];
  T = new PetscScalar[Mz];
  Tnew = new PetscScalar[Mz];

  Tb = new PetscScalar[Mbz];
  Tbnew = new PetscScalar[Mbz];

  IceModelVec2    *pccTs, *pccsbt, *pccsbmf;

  if (atmosPCC != PETSC_NULL) {
    // call sets pccTs to point to IceModelVec2 with current surface temps
    ierr = atmosPCC->updateSurfTempAndProvide(
              grid.year, dtTempAge / secpera, (void*)(&info_atmoscoupler), pccTs); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");
  }
  ierr = pccTs->get_array(Ts);  CHKERRQ(ierr);

  if (oceanPCC != PETSC_NULL) {
    ierr = oceanPCC->updateShelfBaseTempAndProvide(
              grid.year, dt / secpera, (void*)(&info_oceancoupler), pccsbt); CHKERRQ(ierr);
    ierr = oceanPCC->updateShelfBaseMassFluxAndProvide(
              grid.year, dt / secpera, (void*)(&info_oceancoupler), pccsbmf); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: oceanPCC == PETSC_NULL");
  }
  ierr = pccsbt->get_array(Tshelfbase);  CHKERRQ(ierr);
  ierr = pccsbmf->get_array(bmr_float);  CHKERRQ(ierr);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vHmelt.get_array(Hmelt); CHKERRQ(ierr);
  ierr = vbasalMeltRate.get_array(basalMeltRate); CHKERRQ(ierr);
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vRb.get_array(Rb); CHKERRQ(ierr);
  ierr = vGhf.get_array(Ghf); CHKERRQ(ierr);

  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = Sigma3.begin_access(); CHKERRQ(ierr);
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = Tnew3.begin_access(); CHKERRQ(ierr);

  ierr = Tb3.begin_access(); CHKERRQ(ierr);

  // counts unreasonably low temperature values; should it be deprecated?
  PetscInt myLowTempCount = 0;  

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
    
      // this is bulge limit constant in K; is max amount by which ice
      //   or bedrock can be lower than surface temperature
      const PetscScalar bulgeMax   = 15.0; 

      // this should *not* be replaced by call to grid.kBelowHeightEQ():
      const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/dzEQ));

      // if isMarginal then only do vertical conduction for ice (i.e. ignore advection
      // and strain heating if isMarginal)
      const bool isMarginal = checkThinNeigh(H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                             H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);
      
      ierr = Tb3.getValColumn(i,j,Mbz,zblevEQ,Tb); CHKERRQ(ierr);

      if (grid.vertical_spacing == EQUAL) {
        ierr = u3.getValColumnPL(i,j,Mz,zlevEQ,u); CHKERRQ(ierr);
        ierr = v3.getValColumnPL(i,j,Mz,zlevEQ,v); CHKERRQ(ierr);
        ierr = w3.getValColumnPL(i,j,Mz,zlevEQ,w); CHKERRQ(ierr);
        ierr = Sigma3.getValColumnPL(i,j,Mz,zlevEQ,Sigma); CHKERRQ(ierr);
        ierr = T3.getValColumnPL(i,j,Mz,zlevEQ,T); CHKERRQ(ierr);
      } else {
        // slower, but right for not-equal spaced
        ierr = u3.getValColumnQUAD(i,j,Mz,zlevEQ,u); CHKERRQ(ierr);
        ierr = v3.getValColumnQUAD(i,j,Mz,zlevEQ,v); CHKERRQ(ierr);
        ierr = w3.getValColumnQUAD(i,j,Mz,zlevEQ,w); CHKERRQ(ierr);
        ierr = Sigma3.getValColumnQUAD(i,j,Mz,zlevEQ,Sigma); CHKERRQ(ierr);
        ierr = T3.getValColumnQUAD(i,j,Mz,zlevEQ,T); CHKERRQ(ierr);
      }

      // go through column and find appropriate lambda for BOMBPROOF
      PetscScalar lambda = 1.0;  // start with centered implicit for more accuracy
      for (PetscInt k = 1; k < ks; k++) {   
        const PetscScalar denom = (PetscAbs(w[k]) + 0.000001/secpera)
                                  * ice->rho * ice->c_p * dzEQ;  
        lambda = PetscMin(lambda, 2.0 * ice->k / denom);
      }
      if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1

      // set up and solve the system for this column; melting not addressed yet
      if (k0+ks>0) {
        ierr = system.tempColumnSetUpAndSolve(
                         i, j, ks, isMarginal, lambda,  T, Tb, u, v, w, Sigma,
                         Ghf[i][j], Ts[i][j], mask[i][j], Tshelfbase[i][j], Rb[i][j],
                         T3, &x);
        if (ierr > 0) {
          SETERRQ3(2,
            "Tridiagonal solve failed at (%d,%d) with zero pivot position %d.\n",
            i, j, ierr);
        } else { CHKERRQ(ierr); }
      }

      // insert bedrock solution; check for too low below
      for (PetscInt k=0; k < k0; k++) {
        Tbnew[k] = x[k];
      }

      // prepare for melting/refreezing
      PetscScalar Hmeltnew = Hmelt[i][j];
      
      // insert solution for generic ice segments
      for (PetscInt k=1; k <= ks; k++) {
        if (allowAboveMelting == PETSC_TRUE) {
          Tnew[k] = x[k0 + k];
        } else {
          const PetscScalar depth = H[i][j] - zlevEQ[k];
          const PetscScalar Tpmp = ice->meltingTemp - ice->beta_CC_grad * depth;
          if (x[k0 + k] > Tpmp) {
            Tnew[k] = Tpmp;
            PetscScalar Texcess = x[k0 + k] - Tpmp; // always positive
            excessToFromBasalMeltLayer(rho_c_I, zlevEQ[k], dzEQ, &Texcess, &Hmeltnew);
            // Texcess  will always come back zero here; ignore it
          } else {
            Tnew[k] = x[k0 + k];
          }
        }
        if (Tnew[k] < globalMinAllowedTemp) {
           ierr = PetscPrintf(PETSC_COMM_SELF,
              "  [[too low (<200) ice segment temp T = %f at %d,%d,%d;"
              " proc %d; mask=%f; w=%f]]\n",
              Tnew[k],i,j,k,grid.rank,mask[i][j],w[k]*secpera); CHKERRQ(ierr);
           myLowTempCount++;
        }
        if (Tnew[k] < Ts[i][j] - bulgeMax) {
          Tnew[k] = Ts[i][j] - bulgeMax;  bulgeCount++;   }
      }
      
      // insert solution for ice/rock interface (or base of ice shelf) segment
      if (ks > 0) {
        if (allowAboveMelting == PETSC_TRUE) {
          Tnew[0] = x[k0];
        } else {  // compute diff between x[k0] and Tpmp; melt or refreeze as appropriate
          const PetscScalar Tpmp = ice->meltingTemp - ice->beta_CC_grad * H[i][j];
          PetscScalar Texcess = x[k0] - Tpmp; // positive or negative
          if (PismModMask(mask[i][j]) == MASK_FLOATING) {
             // when floating, only half a segment has had its temperature raised
             // above Tpmp
             excessToFromBasalMeltLayer(rho_c_I/2, 0.0, dzEQ, &Texcess, &Hmeltnew);
          } else {
             excessToFromBasalMeltLayer(rho_c_av, 0.0, dzEQ, &Texcess, &Hmeltnew);
          }
          Tnew[0] = Tpmp + Texcess;
          if (Tnew[0] > (Tpmp + 0.00001)) {
            SETERRQ(1,"updated temperature came out above Tpmp");
          }
        }
        if (Tnew[0] < globalMinAllowedTemp) {
           ierr = PetscPrintf(PETSC_COMM_SELF,
              "  [[too low (<200) ice/bedrock segment temp T = %f at %d,%d;"
              " proc %d; mask=%f; w=%f]]\n",
              Tnew[0],i,j,grid.rank,mask[i][j],w[0]*secpera); CHKERRQ(ierr);
           myLowTempCount++;
        }
        if (Tnew[0] < Ts[i][j] - bulgeMax) {
          Tnew[0] = Ts[i][j] - bulgeMax;   bulgeCount++;   }
      } else {
        Hmeltnew = 0.0;
      }
      
      // we must agree on redundant values T(z=0) at top of bedrock and at bottom of ice
      if (ks > 0) {
        Tbnew[k0] = Tnew[0];
      } else {
        if (PismModMask(mask[i][j]) == MASK_FLOATING) { // top of bedrock sees ocean
          Tbnew[k0] = Tshelfbase[i][j]; // set by PISMOceanCoupler
        } else { // top of bedrock sees atmosphere
          Tbnew[k0] = Ts[i][j];
        }
      }
      // check bedrock solution        
      for (PetscInt k=0; k <= k0; k++) {
        if (Tbnew[k] < globalMinAllowedTemp) {
           ierr = PetscPrintf(PETSC_COMM_SELF,
              "  [[too low (<200) bedrock segment temp T = %f at %d,%d,%d;"
              " proc %d; mask=%f]]\n",
              Tbnew[k],i,j,k,grid.rank,mask[i][j]); CHKERRQ(ierr);
           myLowTempCount++;
        }
        if (Tbnew[k] < Ts[i][j] - bulgeMax) {
          Tbnew[k] = Ts[i][j] - bulgeMax;   bulgeCount++;   }
      }

      // transfer column into Tb3; neighboring columns will not reference!
      ierr = Tb3.setValColumn(i,j,Mbz,zblevEQ,Tbnew); CHKERRQ(ierr);

      // set to air temp above ice
      for (PetscInt k=ks; k<Mz; k++) {
        Tnew[k] = Ts[i][j];
      }

      // transfer column into Tnew3; communication later
      ierr = Tnew3.setValColumnPL(i,j,Mz,zlevEQ,Tnew); CHKERRQ(ierr);

      // basalMeltRate[][] is rate of mass loss at bottom of ice everywhere;
      //   note massContExplicitStep() calls PISMOceanCoupler separately
      if (PismModMask(mask[i][j]) == MASK_FLOATING) {
        // rate of mass loss at bottom of ice shelf;  can be negative (marine freeze-on)
        basalMeltRate[i][j] = bmr_float[i][j]; // set by PISMOceanCoupler
      } else {
        // rate of change of Hmelt[][];  can be negative (till water freeze-on)
        basalMeltRate[i][j] = (Hmeltnew - Hmelt[i][j]) / dtTempAge;
      }

      if (PismModMask(mask[i][j]) == MASK_FLOATING) {
        // eliminate basal lubrication water if floating; 
        Hmelt[i][j] = 0.0;
      } else {
        // limit Hmelt by default max and store
        Hmelt[i][j] = PetscMin(Hmelt_max, Hmeltnew);
      }

    } 
  }
  
  if (myLowTempCount > maxLowTempCount) { SETERRQ(1,"too many low temps"); }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vHmelt.end_access(); CHKERRQ(ierr);
  ierr = vRb.end_access(); CHKERRQ(ierr);
  ierr = vGhf.end_access(); CHKERRQ(ierr);
  ierr = vbasalMeltRate.end_access(); CHKERRQ(ierr);

  ierr = pccTs->end_access(); CHKERRQ(ierr);
  ierr = pccsbt->end_access();  CHKERRQ(ierr);
  ierr = pccsbmf->end_access();  CHKERRQ(ierr);

  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = Sigma3.end_access(); CHKERRQ(ierr);
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = Tnew3.end_access(); CHKERRQ(ierr);
  
  delete [] x;

  delete [] u;  delete [] v;  delete [] w;  delete [] Sigma;
  delete [] T;  delete [] Tb;  delete [] Tbnew;  delete [] Tnew;

  delete [] zlevEQ;   delete [] zblevEQ;

  return 0;
}


//! Compute the melt water which should go to the base if \f$T\f$ is above pressure-melting.
PetscErrorCode IceModel::excessToFromBasalMeltLayer(
                const PetscScalar rho_c, const PetscScalar z, const PetscScalar dz,
                PetscScalar *Texcess, PetscScalar *Hmelt) {

  const PetscScalar darea = grid.dx * grid.dy,
                    dvol = darea * dz,
                    dE = rho_c * (*Texcess) * dvol,
                    massmelted = dE / ice->latentHeat;

  if (allowAboveMelting == PETSC_TRUE) {
    SETERRQ(1,"excessToBasalMeltLayer() called but allowAboveMelting==TRUE");
  }
  if (*Texcess >= 0.0) {
    if (updateHmelt == PETSC_TRUE) {
      // T is at or above pressure-melting temp, so temp needs to be set to 
      // pressure-melting; a fraction of excess energy is turned into melt water at base
      // note massmelted is POSITIVE!
      const PetscScalar FRACTION_TO_BASE
                           = (z < 100.0) ? 0.2 * (100.0 - z) / 100.0 : 0.0;
      // note: ice-equiv thickness:
      *Hmelt += (FRACTION_TO_BASE * massmelted) / (ice->rho * darea);  
    }
    *Texcess = 0.0;
  } else if (updateHmelt == PETSC_TRUE) {  // neither Texcess nor Hmelt need to change 
                                           // if Texcess < 0.0
    // Texcess negative; only refreeze (i.e. reduce Hmelt) if at base and Hmelt > 0.0
    // note ONLY CALLED IF AT BASE!   note massmelted is NEGATIVE!
    if (z > 0.00001) {
      SETERRQ(1, "excessToBasalMeltLayer() called with z not at base and negative Texcess");
    }
    if (*Hmelt > 0.0) {
      const PetscScalar thicknessToFreezeOn = - massmelted / (ice->rho * darea);
      if (thicknessToFreezeOn <= *Hmelt) { // the water *is* available to freeze on
        *Hmelt -= thicknessToFreezeOn;
        *Texcess = 0.0;
      } else { // only refreeze Hmelt thickness of water; update Texcess
        *Hmelt = 0.0;
        const PetscScalar dTemp = ice->latentHeat * ice->rho * (*Hmelt) / (rho_c * dz);
        *Texcess += dTemp;
      }
    } 
    // note: if *Hmelt == 0 and Texcess < 0.0 then Texcess unmolested; temp will go down
  }
  return 0;
}                           


//! Take a semi-implicit time-step for the age equation.  Also check the horizontal CFL for advection.
/*!
The age equation is\f$d\tau/dt = 1\f$, that is,
    \f[ \frac{\partial \tau}{\partial t} + u \frac{\partial \tau}{\partial x}
        + v \frac{\partial \tau}{\partial y} + w \frac{\partial \tau}{\partial z} = 1\f]
where \f$\tau(t,x,y,z)\f$ is the age of the ice and \f$(u,v,w)\f$  is the three dimensional
velocity field.  This equation is hyperbolic (purely advective).  

The boundary condition is that when the ice fell as snow it had age zero.  
That is, \f$\tau(t,x,y,h(t,x,y)) = 0\f$ in accumulation areas, while there is no 
boundary condition elsewhere, as the characteristics go outward in the ablation zone.

If the velocity in the bottom cell of ice is upward (w[i][j][0]>0) then we apply
an age=0 boundary condition.  This is the case where ice freezes on at the base,
either grounded basal ice or marine basal ice.

\latexonly\index{BOMBPROOF!implementation for age equation}\endlatexonly
The numerical method is first-order upwind but the vertical advection term is computed
implicitly.  Thus there is no CFL-type stability condition for that part.

We use equally-spaced vertical grid in the calculation.  Note that the IceModelVec3 
methods getValColumn() and setValColumn() interpolate back and forth between the grid 
on which calculation is done and the storage grid.  Thus the storage grid can be either 
equally spaced or not.
 */
PetscErrorCode IceModel::ageStep(PetscScalar* CFLviol) {
  PetscErrorCode  ierr;

  PetscInt    Mz, dummyMbz;
  PetscScalar dzEQ, dummydz, *zlevEQ, *dummylev;

  ierr = getMzMbzForTempAge(Mz, dummyMbz); CHKERRQ(ierr);

  zlevEQ = new PetscScalar[Mz];
  dummylev = new PetscScalar[dummyMbz];

  ierr = getVertLevsForTempAge(Mz, dummyMbz, dzEQ, dummydz, zlevEQ, dummylev);
     CHKERRQ(ierr);

  delete [] dummylev;

  // constants associated to CFL checking
  const PetscScalar cflx = grid.dx / dtTempAge,
                    cfly = grid.dy / dtTempAge;

  PetscScalar **H, *tau, *u, *v, *w;
  tau = new PetscScalar[Mz];
  u = new PetscScalar[Mz];
  v = new PetscScalar[Mz];
  w = new PetscScalar[Mz];

  PetscScalar *x;  
  x = new PetscScalar[Mz]; // space for solution

  ageSystemCtx system(Mz); // linear system to solve in each column

  // set constants within system of equations which are independent of column
  ierr = system.ageSetConstants(grid.dx,grid.dy,dtTempAge,dzEQ); CHKERRQ(ierr);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = tau3.begin_access(); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = taunew3.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // this should *not* be replaced by a call to grid.kBelowHeightEQ()
      const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/dzEQ));
      if (ks > Mz-1) {
        SETERRQ3(1,
           "ageStep() ERROR: ks = %d too high in ice column;\n"
           "  H[i][j] = %5.4f exceeds Lz = %5.4f\n",
           ks, H[i][j], grid.Lz);
      }

      if (ks == 0) { // if no ice, set the entire column to zero age
                     // and ignore the velocities in that column
        ierr = taunew3.setColumn(i,j,0.0); CHKERRQ(ierr);
      } else { // general case
        ierr = tau3.getValColumnQUAD(i,j,Mz,zlevEQ,tau); CHKERRQ(ierr);
        ierr = u3.getValColumnQUAD(i,j,Mz,zlevEQ,u); CHKERRQ(ierr);
        ierr = v3.getValColumnQUAD(i,j,Mz,zlevEQ,v); CHKERRQ(ierr);
        ierr = w3.getValColumnQUAD(i,j,Mz,zlevEQ,w); CHKERRQ(ierr);

        // age evolution is pure advection (so provides check on temp calculation):
        //   check horizontal CFL conditions at each point
        for (PetscInt k=0; k<ks; k++) {
          if (PetscAbs(u[k]) > cflx)  *CFLviol += 1.0;
          if (PetscAbs(v[k]) > cfly)  *CFLviol += 1.0;
        }

        // set up and solve the system for this column
        ierr = system.ageColumnSetUpAndSolve(i,j,ks,u,v,w,tau3,&x);
        if (ierr > 0) {
          SETERRQ3(2,
            "Tridiagonal solve failed at (%d,%d) with zero pivot position %d.\n",
            i, j, ierr);
        } else { CHKERRQ(ierr); }

        // x[k] contains age for k=0,...,ks
        for (PetscInt k=ks+1; k<Mz; k++) {
          x[k] = 0.0;  // age of ice above (and at) surface is zero years
        }
        
        // put solution in IceModelVec3
        ierr = taunew3.setValColumnPL(i,j,Mz,zlevEQ,x); CHKERRQ(ierr);
      }
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = tau3.end_access();  CHKERRQ(ierr);
  ierr = u3.end_access();  CHKERRQ(ierr);
  ierr = v3.end_access();  CHKERRQ(ierr);
  ierr = w3.end_access();  CHKERRQ(ierr);
  ierr = taunew3.end_access();  CHKERRQ(ierr);

  delete [] x;  delete [] tau;  delete [] u;  delete [] v;  delete [] w;
  delete [] zlevEQ;

  return 0;
}


bool IceModel::checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                              PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE) {
  const PetscScalar THIN = 100.0;  // thin = (at most 100m thick)
  return (   (E < THIN) || (NE < THIN) || (N < THIN) || (NW < THIN)
          || (W < THIN) || (SW < THIN) || (S < THIN) || (SE < THIN) );
}


/*!
If the storage grid (defined by IceGrid) has equally-spaced vertical, then
the computation in temperatureStep() and ageStep() is done on that grid.  

If IceGrid defines a not equally spaced grid, however, then, internally in temperatureStep()
and ageStep(), we do computation on a fine and equally-spaced grid.  

This method determines the number of levels in the equally-spaced grid used within 
temperatureStep() and ageStep() in either case.  The method getVertLevsForTempAge() sets 
the spacing and the actual levels.

The storage grid may have quite different levels.  The mapping to the storage grid occurs in 
getValColumn(), setValColumn() for the IceModelVec3 or IceModelVec3Bedrock.
 */
PetscErrorCode IceModel::getMzMbzForTempAge(PetscInt &ta_Mz, PetscInt &ta_Mbz) {

  if (grid.vertical_spacing == EQUAL) {
    ta_Mbz = grid.Mbz;
    ta_Mz = grid.Mz;
  } else {
    ta_Mz = 1 + static_cast<PetscInt>(ceil(grid.Lz / grid.dzMIN));
    ta_Mbz = 1 + static_cast<PetscInt>(ceil(grid.Lbz / grid.dzMIN));
  }
  return 0;
}


/*!
See comments for getMzMbzForTempAge().  The arrays ta_zlevEQ and ta_zblevEQ must 
already be allocated arrays of length ta_Mz, ta_Mbz, respectively.
 */
PetscErrorCode IceModel::getVertLevsForTempAge(const PetscInt ta_Mz, const PetscInt ta_Mbz,
                            PetscScalar &ta_dzEQ, PetscScalar &ta_dzbEQ, 
                            PetscScalar *ta_zlevEQ, PetscScalar *ta_zblevEQ) {

  if (grid.vertical_spacing == EQUAL) {
    ta_dzEQ = grid.dzMIN;
    ta_dzbEQ = grid.dzMIN;
    for (PetscInt k = 0; k < ta_Mz; k++) {
      ta_zlevEQ[k] = grid.zlevels[k];
    }
    for (PetscInt k = 0; k < ta_Mbz; k++) {
      ta_zblevEQ[k] = grid.zblevels[k];
    }
  } else {
    // exactly Mz-1 steps for [0,Lz]:
    ta_dzEQ = grid.Lz / ((PetscScalar) (ta_Mz - 1));  
    for (PetscInt k = 0; k < ta_Mz-1; k++) {
      ta_zlevEQ[k] = ((PetscScalar) k) * ta_dzEQ;
    }
    ta_zlevEQ[ta_Mz-1] = grid.Lz;  // make sure it is right on
    if (ta_Mbz > 1) {
      // exactly Mbz-1 steps for [-Lbz,0]:
      ta_dzbEQ = grid.Lbz / ((PetscScalar) (ta_Mbz - 1));  
      for (PetscInt kb = 0; kb < ta_Mbz-1; kb++) {
        ta_zblevEQ[kb] = - grid.Lbz + ta_dzbEQ * ((PetscScalar) kb);
      }
    } else {
      ta_dzbEQ = ta_dzEQ;
    }
    ta_zblevEQ[ta_Mbz-1] = 0.0;  // make sure it is right on
  }  
  return 0;
}

