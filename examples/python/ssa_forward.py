#! /usr/bin/env python
#
# Copyright (C) 2011 David Maxwell and Constantine Khroulev
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
import PISM, math, time

class ssa_from_boot_file(PISM.ssa.SSARun):
  def __init__(self,boot_file,gridsize=None):
    PISM.ssa.SSARun.__init__(self)
    self.grid = PISM.Context().newgrid()
    self.gridsize=gridsize
    self.config = self.grid.config
    self.boot_file = boot_file

  def _setFromOptions(self):
    config = self.config
    
    # FIXME (DAM 4/28/11)
    # These options probably don't belong here.  Seems like IceBasalResistancePlasticLaw
    # should be able to set these for itself.  
    for o in PISM.OptionsGroup(title="Options for pseudo-plastic ice law"):
      # // use pseudo plastic instead of pure plastic; see iMbasal.cc
      config.flag_from_option("pseudo_plastic", "do_pseudo_plastic_till")

      # // power in denominator on pseudo_plastic_uthreshold; typical is q=0.25; q=0 is pure plastic
      config.scalar_from_option("pseudo_plastic_q", "pseudo_plastic_q")
      if PISM.optionsIsSet("-pseudo_plastic_q"):
        config.set_flag("do_pseudo_plastic_till", True)

      # // threshold; at this velocity tau_c is basal shear stress
      config.scalar_from_option("pseudo_plastic_uthreshold", "pseudo_plastic_uthreshold")
      if PISM.optionsIsSet("-pseudo_plastic_uthreshold"):
        config.set_flag("do_pseudo_plastic_till", True);

      # // controls regularization of plastic basal sliding law
      config.scalar_from_option("plastic_reg", "plastic_regularization")

    for o in PISM.OptionsGroup(title="BasalTillStrength"):
      # // plastic_till_c_0 is a parameter in the computation of the till yield stress tau_c
      # // from the thickness of the basal melt water bwat
      # // Note: option is given in kPa.
      config.scalar_from_option("plastic_c0", "till_c_0");

      # // till_pw_fraction is a parameter in the computation of the till yield stress tau_c
      # // from the thickness of the basal melt water bwat
      # // option a pure number (a fraction); no conversion
      config.scalar_from_option("plastic_pwfrac", "till_pw_fraction")


      config.flag_from_option("thk_eff", "thk_eff_basal_water_pressure")

      if PISM.optionsIsSet("-use_ssa_when_grounded"):
        config.scalar_from_option("use_ssa_when_grounded", "use_ssa_when_grounded")
      else:
        # We're using the SSA, and PISM.PISMYieldStress needs to know this
        # to compute yeild stresses.
        config.set_flag("use_ssa_when_grounded",True);

  def _initGrid(self):
    # FIXME: allow specification of Mx and My different from what's
    # in the boot_file.
    PISM.util.init_grid_from_file(self.grid,self.boot_file,
                                  periodicity=PISM.XY_PERIODIC);

  def _initPhysics(self):
    config = self.config
    basal = PISM.IceBasalResistancePlasticLaw(
           config.get("plastic_regularization") / PISM.secpera,
           config.get_flag("do_pseudo_plastic_till"),
           config.get("pseudo_plastic_q"),
           config.get("pseudo_plastic_uthreshold") / PISM.secpera);

    enthalpyconverter = PISM.EnthalpyConverter(config)
    if PISM.getVerbosityLevel() >3:
      enthalpyconverter.viewConstants(PETSc.Viewer.STDOUT())

    if PISM.optionsIsSet("-ssa_glen"):
      ice = PISM.CustomGlenIce(self.grid.com,"",config,enthalpyconverter)
      B_schoof = 3.7e8;     # Pa s^{1/3}; hardness 
      ice.setHardness(B_schoof)
    else:
      ice =  PISM.GPBLDIce(self.grid.com, "", config,enthalpyconverter)
    ice.setFromOptions()

    self.solver.setPhysics(ice,basal,enthalpyconverter)

  def _initSSACoefficients(self):
    # Read PISM SSA related state variables
    solver = self.solver
    solver.allocateCoeffs()

    thickness = solver.thickness; bed = solver.bed; enthalpy = solver.enthalpy
    mask = solver.ice_mask; surface = solver.surface

    # Read in the PISM state variables that are used directly in the SSA solver
    for v in [thickness, bed, enthalpy]:
      v.regrid(self.boot_file,True)
  
    # variables mask and surface are computed from the geometry previously read
    sea_level = 0 # FIXME setFromOption?
    gc = PISM.GeometryCalculator(sea_level,self.solver.ice,self.config)
    gc.compute(bed,thickness,mask,surface)

    # Compute yield stress from PISM state variables
    # (basal melt rate, tillphi, and basal water height)
    grid = self.grid

    bmr   = PISM.util.standardBasalMeltRateVec(grid)
    tillphi = PISM.util.standardTillPhiVec(grid)
    bwat = PISM.util.standardBasalWaterVec(grid)
    for v in [bmr,tillphi,bwat]:
       v.regrid(self.boot_file,True)
    pvars = PISM.PISMVars()
    for v in [thickness,bed,mask,bmr,tillphi,bwat]:
       pvars.add(v)

    yieldstress = PISM.PISMDefaultYieldStress(grid,grid.config)
    yieldstress.init(pvars) 
    yieldstress.basal_material_yield_stress(solver.tauc)



class BasalTillStrength:
  def __init__(self,grid,ice_rho,standard_gravity):
    self.grid = grid

    self.rho_g = ice_rho*standard_gravity

    self.setFromOptions()

    config        = PISM.global_config()
    
    self.till_pw_fraction = config.get("till_pw_fraction")
    self.till_c_0 = config.get("till_c_0") * 1e3 # convert from kPa to Pa
    self.bwat_max = config.get("bwat_max");

    self.usebmr        = config.get_flag("bmr_enhance_basal_water_pressure")
    self.usethkeff     = config.get_flag("thk_eff_basal_water_pressure")

    self.bmr_scale     = config.get("bmr_enhance_scale")
    self.thkeff_reduce = config.get("thk_eff_reduced")
    self.thkeff_H_high = config.get("thk_eff_H_high")
    self.thkeff_H_low  = config.get("thk_eff_H_low")
    

  def setFromOptions(self):
    for o in PISM.OptionsGroup(title="BasalTillStrength"):
      # // plastic_till_c_0 is a parameter in the computation of the till yield stress tau_c
      # // from the thickness of the basal melt water bwat
      # // Note: option is given in kPa.
      config.scalar_from_option("plastic_c0", "till_c_0");

      # // till_pw_fraction is a parameter in the computation of the till yield stress tau_c
      # // from the thickness of the basal melt water bwat
      # // option a pure number (a fraction); no conversion
      config.scalar_from_option("plastic_pwfrac", "till_pw_fraction")


      config.flag_from_option("thk_eff", "thk_eff_basal_water_pressure")

      # # // "friction angle" in degrees
      # config.scalar_from_option("plastic_phi", "default_till_phi")

# The updateYieldStress and getBasalWaterPressure come from iMBasal.

  def updateYieldStress(self,mask,thickness,bwat,bmr,tillphi,tauc):
    config = PISM.global_config()
    till_pw_fraction = self.till_pw_fraction#config.get("till_pw_fraction")
    till_c_0 = self.till_c_0#config.get("till_c_0") * 1e3 # convert from kPa to Pa
    bwat_max = self.bwat_max#config.get("bwat_max");

    rho_g = self.rho_g


    with PISM.util.Access(nocomm=[mask,thickness,bwat,bmr,tillphi],comm=tauc):
      mq = PISM.MaskQuery(mask)
      GHOSTS = self.grid.max_stencil_width;
      for (i,j) in self.grid.points_with_ghosts(nGhosts = GHOSTS):
        if mq.floating_ice(i,j):
          tauc[i,j] = 0
          continue

        H_ij = thickness[i,j]
        if H_ij == 0:
          tauc[i,j] = 1000.0e3;  #// large yield stress of 1000 kPa = 10 bar if no ice
        else: # grounded and there is some ice
          p_over = rho_g * H_ij
          p_w    = self.getBasalWaterPressure( H_ij,
                                         bwat[i,j],bmr[i,j],till_pw_fraction, 
                                         bwat_max)
          N = p_over - p_w #  effective pressure on till
          tauc[i,j] = till_c_0 + N * math.tan((math.pi/180.0) * tillphi[i,j])

  def getBasalWaterPressure( self, thk, bwat, bmr, frac, bwat_max ):  
    if (bwat > bwat_max + 1.0e-6):
      verbPrintf(1,grid.com,
        "PISM ERROR:  bwat = %12.8f exceeds bwat_max = %12.8f\n" +
        "  in IceModel::getBasalWaterPressure()\n", bwat, bwat_max );
      PISM.PISMEnd();

    # the model; note  0 <= p_pw <= frac * p_overburden
    #   because  0 <= bwat <= bwat_max
    p_overburden = self.rho_g * thk; #// FIXME task #7297
    p_pw = frac * (bwat / bwat_max) * p_overburden;

    if (self.usebmr):
      # add to pressure from instantaneous basal melt rate;
      #   note  (additional) <= (1.0 - frac) * p_overburden so
      #   0 <= p_pw <= p_overburden
      p_pw += ( 1.0 - math.exp( - max(0.0,bmr) / self.bmr_scale ) ) \
              * (1.0 - frac) * p_overburden;
    if self.usethkeff:
      # ice thickness is surrogate for distance to margin; near margin the till
      #   is presumably better drained so we reduce the water pressure
      if (thk < self.thkeff_H_high):
        if (thk <= self.thkeff_H_low):
          p_pw *= self.thkeff_reduce
        else:
          # case Hlow < thk < Hhigh; use linear to connect (Hlow, reduced * p_pw)
          #   to (Hhigh, 1.0 * p_w)
          p_pw *= self.thkeff_reduce\
                  + (1.0 - self.thkeff_reduce)\
                      * (thk - self.thkeff_H_low) / (self.thkeff_H_high - self.thkeff_H_low);
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
"""  ssa_forward.py -i IN.nc -Mx number -My number [-o file.nc]
  or (at python prompt)
    run ssa_forward -i IN.nc -Mx number -My number [-o file.nc]
  where:
    -i      IN.nc is input file in NetCDF format: contains PISM-written model state
    -Mx     number of grid points in the x direction
    -My     number of grid points in the y direction
  notes:
    * -i is required
  """

  PISM.show_usage_check_req_opts(com,"ssa_forward",["-i"],usage)

  config = context.config()
  for o in PISM.OptionsGroup(com,"","SSA Forward"):
    boot_file = PISM.optionsString("-i","file to bootstrap from")
    output_file = PISM.optionsString("-o","output file",default="ssa_forward.nc")

  ssa_run = ssa_from_boot_file(boot_file)

  ssa_run.setup()

  solve_t0 = time.clock()
  ssa_run.solve()
  solve_t = time.clock()-solve_t0

  PISM.verbPrintf(2,context.com,"Solve time %g seconds.\n",solve_t)

  ssa_run.write(output_file)
  ssa_run.teardown()
