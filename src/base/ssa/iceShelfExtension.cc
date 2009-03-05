// Copyright (C) 2009 Jed Brown
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

#include "iceShelfExtension.hh"
#include "materials.hh"
#include "pism_const.hh"

IceShelfExtension::IceShelfExtension(MPI_Comm c,const char *pre)
  : comm(c) {
  PetscMemzero(prefix,sizeof(prefix));
  if (pre) {
    PetscStrncpy(prefix,pre,sizeof(prefix));
  }
  // defaults for reference temperature, thickness, strain rates
  T = 263.15;
  H = 50;
  Du = (100/secpera) / 1e5;     // 100 m/a per 100 km
  force_nuH = -1;               // Do not force nuH by default
  cached_viscosity = -1;        // Mark the cached viscosity as stale
  // We set a default ice type here so that we always have ice on which to base our extension.  The user will usually
  // want to set our ice type to whatever they are using in the rest of the model, but this is not critical.
  char pre2[256];
  PetscMemcpy(pre2,prefix,sizeof(pre2));
  PetscStrncat(pre2,"shelfext_",sizeof(pre2));
  private_ice = new CustomGlenIce(c,pre2);
  ice = private_ice;
}

IceShelfExtension::~IceShelfExtension() {
  delete private_ice;
}

PetscErrorCode IceShelfExtension::setIce(const IceType *i) {
  ice = i;
  cached_viscosity = -1;
  return 0;
}

PetscErrorCode IceShelfExtension::setTemperature(PetscScalar temp) {
  T = temp;
  cached_viscosity = -1;
  return 0;
}

PetscErrorCode IceShelfExtension::setThickness(PetscScalar thk) {
  H = thk;
  cached_viscosity = -1;
  return 0;
}

PetscErrorCode IceShelfExtension::setStrainRate(PetscScalar du) {
  Du = du;
  cached_viscosity = -1;
  return 0;
}

PetscErrorCode IceShelfExtension::forceNuH(PetscScalar nuH) {
  force_nuH = nuH;
  cached_viscosity = -1;
  return 0;
}

PetscErrorCode IceShelfExtension::setFromOptions() {
  PetscErrorCode ierr;
  PetscTruth pice;

  pice = (PetscTruth)(ice == private_ice);
  ierr = PetscOptionsBegin(comm,prefix,"IceShelfExtension options",NULL);CHKERRQ(ierr);
  {
    ierr = PetscOptionsReal("-shelfext_T","extension temperature (K)","setTemperature",T,&T,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-shelfext_H","extension thickness (m)","setThickness",H,&H,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-shelfext_Du","extension strain rate (s^{-1})","setStrainRate",Du,&Du,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-shelfext_force_nuH","Just force a particular strength (Pa s, not used if negative)","forceNuH",force_nuH,&force_nuH,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsTruth("-shelfext_use_private_ice","Don't use any ice we are given, instead use the private one","",pice,&pice,NULL);CHKERRQ(ierr);
  }
  if (pice) ice = private_ice;
  if (ice == private_ice) {
    // If we are not using private_ice then someone else is respensible for setting options
    ierr = private_ice->setFromOptions();CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  cached_viscosity = -1;        // Mark the cached viscosity as stale
  return 0;
}

PetscErrorCode IceShelfExtension::printInfo(PetscInt verb) {
  PetscErrorCode ierr;

  if (verb <= getVerbosityLevel()) {
    ierr = view(NULL);CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode IceShelfExtension::view(PetscViewer viewer) {
  PetscErrorCode ierr;
  PetscTruth iascii;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(comm,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"IceShelfExtension object (%s)\n",prefix);CHKERRQ(ierr);
    if (force_nuH > 0) {
      ierr = PetscViewerASCIIPrintf(viewer,"  Forcing integrated viscosity nuH = %g, everything below is being ignored!\n",force_nuH);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,
                                  "  temperature (K) %g    thickness (m) %g   strain rate (s^{-1}) %g\n"
                                  "  We are extending a %s IceType, details on which follow\n",
                                  T,H,Du,ice == private_ice?"private":"non-private (it's likely compatible)");CHKERRQ(ierr);
    {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = ice->view(viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,
                                  "  IceShelfExtension is producing an integrated viscosity of %g (Pa s m)\n"
                                  "  Compare this to the Ritz et al (2001) value of 30 MPa yr for \\bar\\nu = 9.45e14 (Pa s)\n",
                                  viscosity());CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this object");
  }
  return 0;
}

PetscScalar IceShelfExtension::viscosity() {
  if (cached_viscosity < 0) {
    if (force_nuH > 0) {
      cached_viscosity = force_nuH;
    } else {
      // Don't drive the flow law with a degenerate case, I don't know if it's a concern
      const PetscScalar Tvec[3] = {T,T,T};
      const PetscScalar zlevels[3] = {0, 0.5*H, H};
      cached_viscosity = ice->effectiveViscosityColumn(H,3,zlevels,Du,0,0,0,Tvec,Tvec);
    }
  }
  return cached_viscosity;
}

PetscScalar IceShelfExtension::thickness() const { return H; }

