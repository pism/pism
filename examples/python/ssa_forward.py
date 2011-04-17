#! /usr/bin/env python
#
# Copyright (C) 2011 Ed Bueler and Constantine Khroulev and David Maxwell
# 
# This file is part of PISM.
# 
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
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
import PISM, math

class ssa_forward(PISM.ssa.SSATestCase):
  def __init__(self,Mx,My,boot_file):
    self.boot_file = boot_file
    PISM.ssa.SSATestCase.__init__(self,Mx,My)

  def initGrid(self,Mx,My):
    PISM.util.init_grid_from_file(self.grid,self.boot_file,
                                  periodicity=PISM.XY_PERIODIC);

  def initPhysics(self):
    config = self.config
    self.basal = PISM.IceBasalResistancePlasticLaw(
           config.get("plastic_regularization") / PISM.secpera,
           config.get_flag("do_pseudo_plastic_till"),
           config.get("pseudo_plastic_q"),
           config.get("pseudo_plastic_uthreshold") / PISM.secpera);

    if PISM.optionsIsSet("-ssa_glen"):
      self.ice = PISM.CustomGlenIce(self.grid.com,"",config)
      B_schoof = 3.7e8;     # Pa s^{1/3}; hardness 
      self.ice.setHardness(B_schoof)
      self.config.set("epsilon_ssafd", 0.0);  # don't use this lower bound
    else:
      self.ice =  PISM.GPBLDIce(self.grid.com, "", config)
    self.ice.setFromOptions()
    
    self.enthalpyconverter = PISM.EnthalpyConverter(config)
    if PISM.getVerbosityLevel() >3:
      self.enthalpyconverter.viewConstants(PETSc.Viewer.STDOUT())

  def initSSACoefficients(self):
    # Read PISM SSA related state variables
    solver = self.solver


    (thickness,wasSet) =  PISM.optionsRealWasSet("-ssa_thickness","constant thickness")
    if wasSet:
      solver.thickness.set(thickness)
    else:
      solver.thickness.regrid(self.boot_file,True)

    (bed,wasSet) =  PISM.optionsRealWasSet("-ssa_bed", "constant bed elevation")
    if wasSet:
      solver.bed.set(bed)
    else:
      solver.bed.regrid(self.boot_file,True)



    (enth,wasSet) = PISM.optionsRealWasSet("-ssa_enth","constant enthalpy")
    if wasSet:
      solver.enthalpy.set(enth)
    else:
      solver.enthalpy.regrid(self.boot_file,True)

    # Compute mask directly from bed, thickness, and sea level
    seaLevel = 0; # FIXME
    ocean_rho = self.config.get("sea_water_density");
    ice_rho = self.ice.rho

    (mask,wasSet) = PISM.optionsIntWasSet("-ssa_mask","constant mask")
    if wasSet:
      solver.ice_mask.set(mask)
    else:
      update_mask(seaLevel, ocean_rho, ice_rho, self.grid, 
                  solver.thickness, solver.bed, solver.ice_mask)

    
    (surface,wasSet) =  PISM.optionsRealWasSet("-ssa_surface","constant surface elevation")
    if wasSet:
      print "setting surface %g" % surface
      solver.surface.set(surface)
    else:
      # Compute surface elevation directly from bed, thickness, mask
      is_dry = self.config.get_flag("is_dry_simulation");
      update_surface_elevation(seaLevel, ocean_rho, ice_rho, is_dry,
                               self.grid, solver.thickness, solver.bed,
                               solver.ice_mask, solver.surface)

    (tauc,wasSet) =  PISM.optionsRealWasSet("-ssa_tauc","constant tauc")
    if wasSet:
      solver.tauc.set(tauc)
    else:
      # Compute yield stress from PISM state variables
      # (basal melt rate, tillphi, and basal water height)
      grid = self.grid
      bmr   = PISM.util.standardBasalMeltRateVec(grid)
      tillphi = PISM.util.standardTillPhiVec(grid)
      Hmelt = PISM.util.standardBasalWaterVec(grid)
      for v in [bmr,tillphi,Hmelt]:
        v.regrid(self.boot_file,True)

      standard_gravity = self.config.get("standard_gravity")
      updateYieldStress(grid, ice_rho, standard_gravity,
                        solver.ice_mask, solver.tauc, solver.thickness,
                        Hmelt, bmr, tillphi)


def update_mask(currentSeaLevel, ocean_rho, ice_rho, grid, thickness, 
                bed, mask ):

  with PISM.util.Access([thickness,bed,mask]):
    for (i,j) in grid.points_with_ghosts(nGhosts=2):
      H_ij = thickness[i,j]

      hgrounded = bed[i,j] + H_ij # FIXME task #7297      
      hfloating = currentSeaLevel + (1.0 - ice_rho/ocean_rho) * H_ij;

      is_floating = hfloating > hgrounded + 1.0
      # note: the following implies that ice-free cells with bed evelation
      # exactly at sea level are considered grounded
      is_grounded = not is_floating
      ice_free = H_ij < 0.01;

      if mask[i,j] == PISM.MASK_OCEAN_AT_TIME_0:
        continue

      if is_floating:
        if ice_free:
          mask[i,j] = PISM.MASK_ICE_FREE_OCEAN
        else:
          mask[i,j] = PISM.MASK_FLOATING

      if is_grounded:
        if ice_free:
          mask[i,j] = PISM.MASK_ICE_FREE_BEDROCK;
        else:
          mask[i,j] = PISM.MASK_GROUNDED

def update_surface_elevation(seaLevel,ocean_rho,ice_rho,is_dry_simulation,
                             grid,thickness,bed,mask,surface):
  with PISM.util.Access([surface,thickness,bed,mask]):
    for (i,j) in grid.points_with_ghosts(nGhosts=2):
      H_ij = thickness[i,j]
      if H_ij<0:
        raise Exception("Thickness negative at point i=%d j=%d" % (i,j))

      hgrounded = bed[i,j]+H_ij
      hfloating = seaLevel + (1.0-ice_rho/ocean_rho)*H_ij

      if is_dry_simulation:
        surface[i,j] = hgrounded
        continue

      if mask[i,j] == PISM.MASK_OCEAN_AT_TIME_0:
        surface[i,j] = hfloating
        continue

      if mask.is_floating(i,j):
        surface[i,j] = hfloating
      else:
        surface[i,j] = hgrounded


def updateYieldStress(grid,ice_rho,standard_gravity,
                      mask,tauc,thickness,Hmelt,bmr,tillphi):
  config = PISM.global_config()
  till_pw_fraction = config.get("till_pw_fraction")
  till_c_0 = config.get("till_c_0") * 1e3 # convert from kPa to Pa
  hmelt_max = config.get("hmelt_max");

  with PISM.util.Access([mask,tauc,thickness,Hmelt,bmr,tillphi]):
    GHOSTS = grid.max_stencil_width;
    for (i,j) in grid.points_with_ghosts(nGhosts = GHOSTS):
      if mask.is_floating(i,j):
        tauc[i,j] = 0
        continue

      H_ij = thickness[i,j]
      if H_ij == 0:
        tauc[i,j] = 1000.0e3;  #// large yield stress of 1000 kPa = 10 bar if no ice
      else: # grounded and there is some ice
        p_over = ice_rho * standard_gravity * H_ij
        p_w    = getBasalWaterPressure(ice_rho,standard_gravity,H_ij,
                                       Hmelt[i,j],bmr[i,j],till_pw_fraction, 
                                       hmelt_max)
        N = p_over - p_w #  effective pressure on till
        tauc[i,j] = till_c_0 + N * math.tan((math.pi/180.0) * tillphi[i,j])

def getBasalWaterPressure( ice_rho, standard_gravity, thk, bwat, bmr, 
                           frac, hmelt_max ):  
  if (bwat > hmelt_max + 1.0e-6):
    verbPrintf(1,grid.com,
      "PISM ERROR:  bwat = %12.8f exceeds hmelt_max = %12.8f\n" +
      "  in IceModel::getBasalWaterPressure()\n", bwat, hmelt_max );
    PISM.PISMEnd();

  config        = PISM.global_config()
  usebmr        = config.get_flag("bmr_enhance_basal_water_pressure")
  usethkeff     = config.get_flag("thk_eff_basal_water_pressure")

  bmr_scale     = config.get("bmr_enhance_scale"),
  thkeff_reduce = config.get("thk_eff_reduced"),
  thkeff_H_high = config.get("thk_eff_H_high"),
  thkeff_H_low  = config.get("thk_eff_H_low");

  # the model; note  0 <= p_pw <= frac * p_overburden
  #   because  0 <= bwat <= hmelt_max
  p_overburden = ice_rho * standard_gravity * thk; #// FIXME task #7297
  p_pw = frac * (bwat / hmelt_max) * p_overburden;

  if (usebmr):
    # add to pressure from instantaneous basal melt rate;
    #   note  (additional) <= (1.0 - frac) * p_overburden so
    #   0 <= p_pw <= p_overburden
    p_pw += ( 1.0 - math.exp( - max(0.0,bmr) / bmr_scale ) ) \
            * (1.0 - frac) * p_overburden;
  if usethkeff:
    # ice thickness is surrogate for distance to margin; near margin the till
    #   is presumably better drained so we reduce the water pressure
    if (thk < thkeff_H_high):
      if (thk <= thkeff_H_low):
        p_pw *= thkeff_reduce
      else:
        # case Hlow < thk < Hhigh; use linear to connect (Hlow, reduced * p_pw)
        #   to (Hhigh, 1.0 * p_w)
        p_pw *= thkeff_reduce\
                + (1.0 - thkeff_reduce)\
                    * (thk - thkeff_H_low) / (thkeff_H_high - thkeff_H_low);
  return p_pw;


# The main code for a run follows:
if __name__ == '__main__':
  context = PISM.Context()
  com = context.com

  PISM.set_abort_on_sigint(True)

  PISM.verbosityLevelFromOptions()
  PISM.verbPrintf(2,PISM.Context().com,"SSA forward model.\n")
  PISM.stop_on_version_option()
  usage = \
  """  ssa_forward -i IN.nc -Mx number -My number [-o file.nc]
  where:
    -i      IN.nc is input file in NetCDF format: contains PISM-written model state
    -Mx     number of grid points in the x direction
    -My     number of grid points in the y direction
  notes:
    * -i is required
  """

  PISM.show_usage_check_req_opts(com,"ssa_forward",["-i"],usage)

  config = context.config()
  for o in PISM.OptionsGroup(com,"","SSA Forwawrd"):
    Mx = PISM.optionsInt("-Mx","Number of grid points in the X-direction",default=None)
    My = PISM.optionsInt("-My","Number of grid points in the X-direction",default=None)
    boot_file = PISM.optionsString("-i","file to bootstrap from")
    output_file = PISM.optionsString("-o","output file",default="ssa_forward.nc")

  test_case = ssa_forward(Mx,My,boot_file)
  test_case.grid.printInfo(1)
  test_case.grid.periodicity = PISM.XY_PERIODIC
  test_case.solve()
  test_case.write(output_file)
