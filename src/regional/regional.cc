// Copyright (C) 2010, 2011, 2012, 2013, 2014 Ed Bueler, Daniella DellaGiustina, Constantine Khroulev, and Andy Aschwanden
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

#include "regional.hh"
#include "enthalpyConverter.hh"

namespace pism {

PetscErrorCode SIAFD_Regional::init(Vars &vars) {
  PetscErrorCode ierr;

  ierr = SIAFD::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"  using the regional version of the SIA solver...\n"); CHKERRQ(ierr);

  no_model_mask = vars.get_2d_mask("no_model_mask");
  usurfstore    = vars.get_2d_scalar("usurfstore");

  return 0;
}

PetscErrorCode SIAFD_Regional::compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;

  ierr = SIAFD::compute_surface_gradient(h_x, h_y); CHKERRQ(ierr);

  IceModelVec2Int &nmm = *no_model_mask;
  IceModelVec2S &hst = *usurfstore; // convenience

  const int Mx = grid.Mx, My = grid.My;
  const double dx = grid.dx, dy = grid.dy;  // convenience

  IceModelVec::AccessList list;
  list.add(h_x);
  list.add(h_y);
  list.add(nmm);
  list.add(hst);

  for (PointsWithGhosts p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // x-component, i-offset
    if (nmm(i, j) > 0.5 || nmm(i + 1, j) > 0.5) {

      if (i < 0 || i + 1 > Mx - 1) {
        h_x(i, j, 0) = 0.0;
      } else {
        h_x(i, j, 0) = (hst(i + 1, j) - hst(i, j)) / dx;
      }
    }

    // x-component, j-offset
    if (nmm(i - 1, j + 1) > 0.5 || nmm(i + 1, j + 1) > 0.5 ||
        nmm(i - 1, j)     > 0.5 || nmm(i + 1, j)     > 0.5) {

      if (i - 1 < 0 || j + 1 > My - 1 || i + 1 > Mx - 1) {
        h_x(i, j, 1) = 0.0;
      } else {
        h_x(i, j, 1) = ( + hst(i + 1, j + 1) + hst(i + 1, j)
                         - hst(i - 1, j + 1) - hst(i - 1, j)) / (4.0 * dx);
      }

    }

    // y-component, i-offset
    if (nmm(i, j + 1) > 0.5 || nmm(i + 1, j + 1) > 0.5 ||
        nmm(i, j - 1) > 0.5 || nmm(i + 1, j - 1) > 0.5) {
      if (i < 0 || j + 1 > My - 1 || i + 1 > Mx - 1 || j - 1 < 0) {
        h_y(i, j, 0) = 0.0;
      } else {
        h_y(i, j, 0) = ( + hst(i + 1, j + 1) + hst(i, j + 1)
                         - hst(i + 1, j - 1) - hst(i, j - 1)) / (4.0 * dy);
      }
    }

    // y-component, j-offset
    if (nmm(i, j) > 0.5 || nmm(i, j + 1) > 0.5) {
        
      if (j < 0 || j + 1 > My - 1) {
        h_y(i, j, 1) = 0.0;
      } else {
        h_y(i, j, 1) = (hst(i, j + 1) - hst(i, j)) / dy;
      }
    }

  }

  return 0;
}

SSAFD_Regional::SSAFD_Regional(IceGrid &g, EnthalpyConverter &e, const Config &c)
  : SSAFD(g, e, c) {
  // empty
}

SSAFD_Regional::~SSAFD_Regional() {
  // empty
}

PetscErrorCode SSAFD_Regional::init(Vars &vars) {
  PetscErrorCode ierr;

  ierr = SSAFD::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"  using the regional version of the SSA solver...\n"); CHKERRQ(ierr);

  if (config.get_flag("ssa_dirichlet_bc")) {
    ierr = verbPrintf(2,grid.com,"  using stored SSA velocities as Dirichlet B.C. in the no_model_strip...\n"); 
    CHKERRQ(ierr);
  }

  no_model_mask = vars.get_2d_mask("no_model_mask");
  usurfstore    = vars.get_2d_scalar("usurfstore");
  thkstore      = vars.get_2d_scalar("thkstore");

  return 0;
}

PetscErrorCode SSAFD_Regional::compute_driving_stress(IceModelVec2V &result) {
  PetscErrorCode ierr;

  ierr = SSAFD::compute_driving_stress(result); CHKERRQ(ierr);

  IceModelVec2Int &nmm = *no_model_mask;

  IceModelVec::AccessList list;
  list.add(result);
  list.add(nmm);
  list.add(*usurfstore);
  list.add(*thkstore);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double pressure = EC.getPressureFromDepth((*thkstore)(i,j));
    if (pressure <= 0) {
      pressure = 0;
    }

    if (nmm(i, j) > 0.5 || nmm(i - 1, j) > 0.5 || nmm(i + 1, j) > 0.5) {
      if (i - 1 < 0 || i + 1 > grid.Mx - 1) {
        result(i, j).u = 0;
      } else {
        result(i, j).u = - pressure * usurfstore->diff_x(i,j);
      }
    }

    if (nmm(i, j) > 0.5 || nmm(i, j - 1) > 0.5 || nmm(i, j + 1) > 0.5) {
      if (j - 1 < 0 || j + 1 > grid.My - 1) {
        result(i, j).v = 0;
      } else {
        result(i, j).v = - pressure * usurfstore->diff_y(i,j);
      }
    }
  }

  return 0;
}

PetscErrorCode RegionalDefaultYieldStress::init(Vars &vars) {
  PetscErrorCode ierr;
  int v = getVerbosityLevel(); // turn off second, redundant init message
  ierr = setVerbosityLevel(1); CHKERRQ(ierr);
  ierr = MohrCoulombYieldStress::init(vars); CHKERRQ(ierr);
  ierr = setVerbosityLevel(v); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
    "  using the regional version with strong till in no_model_mask==1 area ...\n");
    CHKERRQ(ierr);

  no_model_mask = vars.get_2d_mask("no_model_mask");

  return 0;
}


PetscErrorCode RegionalDefaultYieldStress::basal_material_yield_stress(IceModelVec2S &result) {
  PetscErrorCode ierr;
  
  // do whatever you normally do
  ierr = MohrCoulombYieldStress::basal_material_yield_stress(result); CHKERRQ(ierr);

  // now set result=tauc to a big value in no_model_strip
  IceModelVec::AccessList list;
  list.add(*no_model_mask);
  list.add(result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*no_model_mask)(i,j) > 0.5) {
      result(i,j) = 1000.0e3;  // large yield stress of 1000 kPa = 10 bar
    }
  }
  return 0;
}


} // end of namespace pism
