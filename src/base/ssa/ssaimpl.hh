// Copyright (C) 2009 Jed Brown and Ed Bueler
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

#ifndef _PISM_SSAIMPL_HH
#define _PISM_SSAIMPL_HH

#include "pismssa.hh"
#include "petscda.h"
#include "petscsnes.h"

struct SSANode {
  PetscScalar x,y;
};

struct _SSAOps {
  PetscErrorCode (*Solve)(SSA);
  PetscErrorCode (*View)(SSA,PetscViewer);
  PetscErrorCode (*Destroy)(SSA);
  PetscErrorCode (*SetFromOptions)(SSA);
  PetscErrorCode (*SetUp)(SSA);
};

enum PismSetupState { SETUP_GREEN, SETUP_STALE, SETUP_CURRENT };

#define PISM_SSA_USE_DIMENSIONAL 0

// Manage nondimensionalization
class PismRef {
public:
  PismRef() { SetUp(); }
  void SetUp() {
#if PISM_SSA_USE_DIMENSIONAL
    length = 1;
    height = 1;
    time = 1;
    pressure = 1;
#else
    // The values here have no effect on results, they are chosen so that initial residuals are O(1)
    length = 1e4;
    height = 1e2;
    time = 10 * secpera;        // This cancels out and has no effect on residual norms
    //stress = 910 * 9.81 * height * (height/length); // You can think of this as choosing units of force
    pressure = 910 * 9.81 * height;
#endif
  }
  PetscReal Length() const { return length; }
  PetscReal Area() const { return length*length; }
  PetscReal Height() const { return height; }
  PetscReal Time() const { return time; }
  PetscReal Velocity() const { return length/time; }
  PetscReal VerticalVelocity() const { return height/time; }
  PetscReal StrainRate() const { return 1/time; }
#if PISM_SSA_USE_DIMENSIONAL
  PetscReal Velocity2() const { PetscReal v = Velocity(); return 1.0; }
  PetscReal StrainRate2() const { PetscReal s = StrainRate(); return 1.0; }
#else
  PetscReal Velocity2() const { PetscReal v = Velocity(); return v*v; }
  PetscReal StrainRate2() const { PetscReal s = StrainRate(); return s*s; }
#endif
  PetscReal Slope() const { return height / length; }
  PetscReal Pressure() const { return pressure; }
  PetscReal DrivingStress() const { return Pressure() * Slope(); }
  PetscReal IntegratedViscosity() const { return DrivingStress() * Length() / StrainRate(); }
  PetscReal Drag() const { return DrivingStress() / Velocity(); }
  PetscErrorCode View(PetscViewer viewer) const; // implemented in ssa/ssa.cc
private:
  PetscReal length,height,time,pressure;
};

struct _p_SSA {
  PETSCHEADER(struct _SSAOps);
  IceGrid          *grid;
  PismRef           ref;
  DA                da;
  IceType          *ice;
  PlasticBasalType *basal;
  SeaWaterType     *ocean;
  IceModelVec2     *siaVel;     // Points at start of array of two vectors with SIA values on staggered grid (uvbar)
  IceModelVec2     *mask,*H,*h,*tauc;
  IceModelVec3     *T;
  PismSetupState    setupcalled; // 0 the first time around, 1 when field are updated, 2 when everything is current
  PetscReal         fictitious_nuH,cutoff_thickness;
  PetscReal         regularizingVelocitySchoof,regularizingLengthSchoof,regSchoof;
  PetscTruth        initialGuessNonzero;
  DAPeriodicType    wrap;
  Vec               x,r;
  Mat               J;
  void             *data;
};

// The following are macros (instead of inline functions) so that error handling is less cluttered.  They should be
// replaced with empty macros when in optimized mode.

#define PismValidVelocity(U) do {                               \
    if (!(-   1e+5 < (U).x && (U).x < 1e+5                      \
          && -1e+5 < (U).y && (U).y < 1e+5))                    \
      SETERRQ3(1,"Invalid velocity (%g,%g) not within %g",      \
               (U).x,(U).y,1e+5);                               \
  } while (0)

#define PismValidStrainRate(Du) do {                                    \
    if (!(-   1e+5 < (Du)[0] && (Du)[0] < 1e+5                          \
          && -1e+5 < (Du)[1] && (Du)[1] < 1e+5                          \
          && -1e+5 < (Du)[2] && (Du)[2] < 1e+5))                        \
      SETERRQ4(1,"Invalid Strain Rate (%g,%g,%g) not within %g",        \
               (Du)[0],(Du)[1],(Du)[2],1e+5);                           \
  } while (0)

#define PismValidStress2(f) do {                                        \
    if (!(-   1e4 < (f).x && (f).x < 1e4                                \
          && -1e4 < (f).y && (f).y < 1e4))                              \
      SETERRQ3(1,"Invalid Stress residual (%g,%g) not within %g",       \
               (f).x,(f).y,1e4);                                        \
  } while (0)

#define PismValidFriction(b) do {                                       \
    if (!(0 <= (b) && (b) < 1e25))                                      \
      SETERRQ2(1,"Invalid friction %g not within [0,%g]",(b),1e25);     \
  } while (0)

#endif /* _PISM_SSAIMPL_HH */
