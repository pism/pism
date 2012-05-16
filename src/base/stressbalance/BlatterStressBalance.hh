// Copyright (C) 2010-2012 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _BLATTERSTRESSBALANCE_H_
#define _BLATTERSTRESSBALANCE_H_

#include "PISMStressBalance.hh"
#include "THI.hh"
#include "IceGrid.hh"

//! Blatter-Pattyn stress balance based on Jed Brown's PETSc tutorial ex48.c (Brown et al. 2011).
/*!
Toy hydrostatic ice flow with multigrid in 3D

Solves the hydrostatic (aka Blatter/Pattyn/First Order) equations for ice sheet
flow using multigrid.  The ice uses a power-law rheology with Glen exponent 3
(corresponds to p=4/3 in a p-Laplacian).  The focus is on ISMIP-HOM experiments
which assume periodic boundary conditions in the x- and y-directions.

Equations are rescaled so that the domain size and solution are O(1), details of
this scaling can be controlled by the options -units_meter, -units_second, and
-units_kilogram.

A VTK StructuredGrid output file can be written using the option
\c -ex48_o \c filename.vts

The equations for horizontal velocity \f$(u,v)\f$ are
\f{align*}
  - [\eta (4 u_x + 2 v_y)]_x - [\eta (u_y + v_x)]_y - [\eta u_z]_z + \rho g s_x &= 0 \\
  - [\eta (4 v_y + 2 u_x)]_y - [\eta (u_y + v_x)]_x - [\eta v_z]_z + \rho g s_y &= 0
\f}
where
  \f[\eta = B/2 (\epsilon + \gamma)^{(p-2)/2}.\f]
is the nonlinear effective viscosity with regularization epsilon and hardness
parameter B, written in terms of the second invariant
  \f[\gamma = u_x^2 + v_y^2 + u_x v_y + (1/4) (u_y + v_x)^2 + (1/4) u_z^2 + (1/4) v_z^2.\f]

The surface boundary conditions are the natural conditions.  The basal boundary
conditions are either no-slip, or Navier (linear) slip with spatially variant
friction coefficient \f$\beta^2\f$.

In the code, the equations for \f$(u,v)\f$ are multiplied through by \f$1/(\rho g)\f$
so that residuals are O(1).

The discretization is Q1 finite elements, managed by a DA.  The grid is never
distorted in the map \f$(x,y)\f$ plane, but the bed and surface may be bumpy.
This is handled as usual in FEM, through the Jacobian of the coordinate
transformation from a reference element to the physical element.

Since ice-flow is tightly coupled in the z-direction (within columns), the DA is
managed specially so that columns are never distributed, and are always
contiguous in memory.  This amounts to reversing the meaning of X,Y,Z compared
to the DA's internal interpretation, and then indexing as vec[i][j][k].  The
exotic coarse spaces require 2D DAs which are made to use compatible domain
decomposition relative to the 3D DAs.

See the source code $PETSC_DIR/src/snes/examples/tutorials/ex48.c for
compile-time options.
 */
class BlatterStressBalance : public PISMStressBalance
{
public:
  BlatterStressBalance(IceGrid &g, PISMOceanModel *ocean_model, const NCConfigVariable &conf);

  virtual ~BlatterStressBalance();

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode set_boundary_conditions(IceModelVec2Int &/*locations*/,
                                                 IceModelVec2V &/*velocities*/)
  { SETERRQ(grid.com, 1,"not clear yet how to do this"); return 0; }

  virtual PetscErrorCode update(bool fast);

  virtual PetscErrorCode get_advective_2d_velocity(IceModelVec2V* &result)
  { result = &vertically_averaged_velocity; return 0; }

  virtual PetscErrorCode get_diffusive_flux(IceModelVec2Stag* &result)
  { PetscErrorCode ierr = result->set(0.0); CHKERRQ(ierr); return 0; }

  virtual PetscErrorCode get_max_diffusivity(PetscReal &Dmax)
  { Dmax = 0; return 0; }

  virtual PetscErrorCode get_max_2d_velocity(PetscReal &maxu, PetscReal &maxv);

  virtual PetscErrorCode get_3d_velocity(IceModelVec3* &u_out, IceModelVec3* &v_out,
                                         IceModelVec3* &w_out);

  virtual PetscErrorCode get_max_3d_velocity(PetscReal &maxu, PetscReal &maxv, PetscReal &maxw);

  virtual PetscErrorCode get_basal_frictional_heating(IceModelVec2S* &/*result*/)
  { SETERRQ(grid.com, 2,"not implemented; implement by doing it"); return 0; }

  virtual PetscErrorCode get_volumetric_strain_heating(IceModelVec3* &/*result*/)
  { SETERRQ(grid.com, 3,"not implemented; implement by getting max"); return 0; }

  virtual PetscErrorCode get_principal_strain_rates(IceModelVec2S &/*result_e1*/,
                                                    IceModelVec2S &/*result_e2*/)
  { SETERRQ(grid.com, 4,"not implemented; implement by vertical average"); return 0; }

  virtual PetscErrorCode extend_the_grid(PetscInt old_Mz);

  virtual PetscErrorCode stdout_report(string &/*result*/)
  { return 0; }  // FIXME: implementation needed

  virtual void add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &/*result*/)
  {  }  // FIXME: implementation needed

  virtual PetscErrorCode define_variables(set<string> /*vars*/, const PIO &/*nc*/,
                                          PISM_IO_Type /*nctype*/)
  { return 0; }                 // FIXME: implementation needed

  virtual PetscErrorCode write_variables(set<string> /*vars*/, string /*filename*/)
  {  return 0; }  // FIXME: implementation needed

protected:
  IceModelVec3 u, v, Sigma;
  IceModelVec2V vertically_averaged_velocity;
  IceModelVec2S basal_frictional_heating;

  IceModelVec2S *topg, *usurf, *tauc;

  THI thi;
  DMMG *dmmg;

  PetscErrorCode allocate_blatter();
  PetscErrorCode deallocate_blatter();
  PetscErrorCode mesh_to_regular_grid();
  PetscErrorCode setup();

  PetscErrorCode compute_basal_frictional_heating();
  PetscErrorCode compute_volumetric_strain_heating();
};

#endif /* _BLATTERSTRESSBALANCE_H_ */

