#! /usr/bin/env python
#
# Copyright (C) 2012 David Maxwell
# 
# This file is part of PISM.
# 
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
# 
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
# 
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np
import os, math

import PISM
import PISM.invert_ssa

def adjustTauc(mask,tauc):
  """Where ice is floating or land is ice-free, tauc should be adjusted to have some preset default values."""

  grid = mask.get_grid()
  high_tauc = grid.config.get("high_tauc")

  with PISM.util.Access(comm=tauc,nocomm=mask):
    mq = PISM.MaskQuery(mask)
    for (i,j) in grid.points():
      if mq.ocean(i,j):
        tauc[i,j] = 0;
      elif mq.ice_free(i,j):
        tauc[i,j] = high_tauc


## Main code starts here
if __name__ == "__main__":
  context = PISM.Context()
  config = context.config
  com = context.com
  PISM.set_abort_on_sigint(True)

  append_mode = False
  PISM.setVerbosityLevel(1)
  for o in PISM.OptionsGroup(context.com,"","test_ssaforward"):
    input_filename = PISM.optionsString("-i","input file")
    inv_data_filename = PISM.optionsString("-inv_data","inverse data file",default=input_filename)
    verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)
    use_tauc_prior = PISM.optionsFlag("-inv_use_tauc_prior","Use tauc_prior from inverse data file as initial guess.",default=False)

  ssarun = PISM.invert_ssa.InvSSAFromInputFile(input_filename,inv_data_filename)
  ssarun.setup()
  
  vecs = ssarun.modeldata.vecs
  grid = ssarun.grid

  # Determine the prior guess for tauc. This can be one of 
  # a) tauc from the input file (default)
  # b) tauc_prior from the inv_datafile if -use_tauc_prior is set
  tauc_prior = PISM.util.standardYieldStressVec(grid,'tauc_prior')
  tauc_prior.set_attrs("diagnostic", "initial guess for (pseudo-plastic) basal yield stress in an inversion", "Pa", "");
  tauc = PISM.util.standardYieldStressVec(grid)
  if use_tauc_prior:
    tauc_prior.regrid(inv_data_filename,critical=True)
  else:
    if not PISM.util.fileHasVariable(input_filename,"tauc"):
      PISM.verbPrintf(1,com,"Initial guess for tauc is not available as 'tauc' in %s.\nYou can provide an initial guess as 'tauc_prior' using the command line option -use_tauc_prior." % input_filename)
      exit(1)
    tauc.regrid(input_filename,True)
    tauc_prior.copy_from(tauc)

  adjustTauc(vecs.ice_mask,tauc_prior)

  # Convert tauc_prior -> zeta_prior
  zeta1 = PISM.IceModelVec2S();
  zeta1.create(grid, "", PISM.kHasGhosts, PISM.util.WIDE_STENCIL)
  ssarun.tauc_param.convertFromTauc(tauc_prior,zeta1)

  ssarun.ssa.linearize_at(zeta1)
    
  (seed,seed_set) = PISM.optionsIntWasSet("-inv_seed","")
  if seed_set:
    np.random.seed(seed+PISM.Context().rank)


######################################################################################################################
# Jacobian design 

  d = PISM.util.randVectorS(grid,1)
  d_proj = PISM.IceModelVec2S()
  d_proj.create(grid,"",PISM.kNoGhosts)
  d_proj.copy_from(d)
  if vecs.has('zeta_fixed_mask'):
    zeta_fixed_mask = vecs.zeta_fixed_mask
    with PISM.util.Access(nocomm=[d_proj,zeta_fixed_mask]):
      for (i,j) in grid.points():
        if zeta_fixed_mask[i,j] != 0:
          d_proj[i,j] = 0;

  r = PISM.util.randVectorV(grid, grid.convert(1.0, "m/year", "m/second"))

  u1 = PISM.IceModelVec2V();
  u1.create(grid,"",PISM.kHasGhosts,PISM.util.WIDE_STENCIL)
  u1.copy_from(ssarun.ssa.solution())

  rhs1 = PISM.IceModelVec2V();
  rhs1.create(grid,"",PISM.kNoGhosts)
  ssarun.ssa.assemble_residual(u1,rhs1)
  
  eps = 1e-8
  zeta2 = PISM.IceModelVec2S();
  zeta2.create(grid, "zeta_prior", PISM.kHasGhosts, PISM.util.WIDE_STENCIL)
  zeta2.copy_from(d_proj)
  zeta2.scale(eps)
  zeta2.add(1,zeta1)
  ssarun.ssa.set_zeta(zeta2)

  rhs2 = PISM.IceModelVec2V();
  rhs2.create(grid,"",PISM.kNoGhosts)
  ssarun.ssa.assemble_residual(u1,rhs2)

  drhs_fd = PISM.IceModelVec2V();
  drhs_fd.create(grid,"",PISM.kNoGhosts)
  drhs_fd.copy_from(rhs2)
  drhs_fd.add(-1,rhs1)
  drhs_fd.scale(1./eps)

  drhs = PISM.IceModelVec2V();
  drhs.create(grid,"",PISM.kNoGhosts)
  ssarun.ssa.apply_jacobian_design(u1,d,drhs)

  d_drhs = PISM.IceModelVec2V();
  d_drhs.create(grid,"",PISM.kNoGhosts)
  
  d_drhs.copy_from(drhs)
  d_drhs.add(-1,drhs_fd)

  n_drhs_fd = drhs_fd.norm(PETSc.NormType.NORM_2)
  n_drhs_l2 = drhs.norm(PETSc.NormType.NORM_2)
  n_drhs_l1 = drhs.norm(PETSc.NormType.NORM_1)
  n_drhs_linf = drhs.norm(PETSc.NormType.NORM_INFINITY)

  n_d_drhs_l2 = d_drhs.norm(PETSc.NormType.NORM_2)
  n_d_drhs_l1 = d_drhs.norm(PETSc.NormType.NORM_1)
  n_d_drhs_linf = d_drhs.norm(PETSc.NormType.NORM_INFINITY)

  PISM.verbPrintf(1,grid.com,"\nTest Jacobian Design (Comparison with finite differences):\n")
  PISM.verbPrintf(1,grid.com,"jacobian_design_transpose(d): l2 norm %.10g; finite difference %.10g\n" % (n_drhs_l2,n_drhs_fd) )
  PISM.verbPrintf(1,grid.com,"relative difference: l2 norm %.10g l1 norm %.10g linf norm %.10g\n" % (n_d_drhs_l2/n_drhs_l2,n_d_drhs_l1/n_drhs_l1,n_d_drhs_linf/n_drhs_linf) )

######################################################################################################################
# Jacobian design transpose

  stencil_width = 1
  u = ssarun.ssa.solution();
  d = PISM.util.randVectorS(grid,1,stencil_width)
  r = PISM.util.randVectorV(grid,1,stencil_width)

  Jd = PISM.IceModelVec2V();
  Jd.create(grid,"",PISM.kNoGhosts)

  JStarR = PISM.IceModelVec2S();
  JStarR.create(grid,"",PISM.kNoGhosts)

  ssarun.ssa.apply_jacobian_design(u,d,Jd)
  ssarun.ssa.apply_jacobian_design_transpose(u,r,JStarR)

  r_global = PISM.IceModelVec2V();
  r_global.create(grid,"",PISM.kNoGhosts)
  r_global.copy_from(r)

  d_global = PISM.IceModelVec2S();
  d_global.create(grid,"",PISM.kNoGhosts)
  d_global.copy_from(d)
  
  ip1 = Jd.get_vec().dot(r_global.get_vec())
  ip2 = JStarR.get_vec().dot(d_global.get_vec())

  PISM.verbPrintf(1,grid.com,"\nTest Jacobian Design Transpose (comparison of r^T*(J*d) versus (J^T*r)^T*d):\n")
  PISM.verbPrintf(1,grid.com,"ip1 %.10g ip2 %.10g\n" % (ip1,ip2) )
  PISM.verbPrintf(1,grid.com,"relative error %.10g\n",abs((ip1-ip2)/ip1))

  ######################################################################################################################
  # Linearization transpose
  
  d = PISM.util.randVectorS(grid,1)
  r = PISM.util.randVectorV(grid,1)

  Td = PISM.IceModelVec2V()
  Td.create(grid,'',PISM.kNoGhosts)
  TStarR = PISM.IceModelVec2S()
  TStarR.create(grid,'',PISM.kNoGhosts)
  
  ssarun.ssa.apply_linearization(d,Td)
  ssarun.ssa.apply_linearization_transpose(r,TStarR)

  ip1 = Td.get_vec().dot(r.get_vec())
  ip2 = TStarR.get_vec().dot(d.get_vec())

  PISM.verbPrintf(1,grid.com,"\nTest Linearization Transpose (comparison of r^t*(T*d) versus (T^t*r)^t*d):\n")
  PISM.verbPrintf(1,grid.com,"ip1 %.10g ip2 %.10g\n" % (ip1,ip2) )
  PISM.verbPrintf(1,grid.com,"relative error %.10g\n",abs((ip1-ip2)/ip1))

  ######################################################################################################################
  # Linearization

  d = PISM.util.randVectorS(grid,1)
  d_proj = PISM.IceModelVec2S()
  d_proj.create(grid,"",PISM.kNoGhosts)
  d_proj.copy_from(d)
  if vecs.has('zeta_fixed_mask'):
    zeta_fixed_mask = vecs.zeta_fixed_mask
    with PISM.util.Access(nocomm=[d_proj,zeta_fixed_mask]):
      for (i,j) in grid.points():
        if zeta_fixed_mask[i,j] != 0:
          d_proj[i,j] = 0;

  ssarun.ssa.linearize_at(zeta1)
  u1 = PISM.IceModelVec2V();
  u1.create(grid,"",PISM.kHasGhosts,stencil_width)
  u1.copy_from(ssarun.ssa.solution())

  Td = PISM.IceModelVec2V()
  Td.create(grid,"",PISM.kHasGhosts,stencil_width)
  ssarun.ssa.apply_linearization(d,Td)

  eps = 1e-8
  zeta2 = PISM.IceModelVec2S();
  zeta2.create(grid, "", PISM.kHasGhosts, PISM.util.WIDE_STENCIL)
  zeta2.copy_from(d_proj)
  zeta2.scale(eps)
  zeta2.add(1,zeta1)
  ssarun.ssa.linearize_at(zeta2)
  u2 = PISM.IceModelVec2V();
  u2.create(grid,"",PISM.kHasGhosts,stencil_width)
  u2.copy_from(ssarun.ssa.solution())
  
  Td_fd = PISM.IceModelVec2V();
  Td_fd.create(grid,"",PISM.kHasGhosts,stencil_width)
  Td_fd.copy_from(u2)
  Td_fd.add(-1,u1)
  Td_fd.scale(1./eps)

  d_Td = PISM.IceModelVec2V()
  d_Td.create(grid,"",PISM.kHasGhosts,stencil_width)
  d_Td.copy_from(Td_fd)
  d_Td.add(-1,Td)

  n_Td_fd = Td_fd.norm(PETSc.NormType.NORM_2)
  n_Td_l2 = Td.norm(PETSc.NormType.NORM_2)
  n_Td_l1 = Td.norm(PETSc.NormType.NORM_1)
  n_Td_linf = Td.norm(PETSc.NormType.NORM_INFINITY)

  n_d_Td_l2 = d_Td.norm(PETSc.NormType.NORM_2)
  n_d_Td_l1 = d_Td.norm(PETSc.NormType.NORM_1)
  n_d_Td_linf = d_Td.norm(PETSc.NormType.NORM_INFINITY)

  PISM.verbPrintf(1,grid.com,"\nTest Linearization (Comparison with finite differences):\n")
  PISM.verbPrintf(1,grid.com,"apply_linearization(d): l2 norm %.10g; finite difference %.10g\n" % (n_Td_l2,n_Td_fd) )
  PISM.verbPrintf(1,grid.com,"relative difference: l2 norm %.10g l1 norm %.10g linf norm %.10g\n" % (n_d_Td_l2/n_Td_l2,n_d_Td_l1/n_Td_l1,n_d_Td_linf/n_Td_linf) )
  
  PISM.verbPrintf(1,grid.com,"\n")

  d_Td_global = PISM.IceModelVec2V()
  d_Td_global.create(grid,"",PISM.kNoGhosts)
  d_Td_global.copy_from(d_Td)
  d_Td_global.scale(1./n_Td_linf)
  d_Td_global.get_vec().view(PETSc.Viewer.DRAW())
