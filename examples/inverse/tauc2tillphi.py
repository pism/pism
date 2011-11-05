#! /usr/bin/env python
#
# Copyright (C) 2011 David Maxwell
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

import sys, os, math, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import PISM

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


    Nmin = 1e45
    with PISM.util.Access(nocomm=[mask,thickness,bwat,bmr,tillphi],comm=tauc):
      mq = PISM.MaskQuery(mask)
      GHOSTS = self.grid.max_stencil_width;
      for (i,j) in self.grid.points_with_ghosts(nGhosts = GHOSTS):
        if mq.floating_ice(i,j):
          tauc[i,j] = 0

        H_ij = thickness[i,j]
        if H_ij == 0:
          tauc[i,j] = 1000.0e3;  #// large yield stress of 1000 kPa = 10 bar if no ice
        else: # grounded and there is some ice
          p_over = rho_g * H_ij
          p_w    = self.getBasalWaterPressure( H_ij,
                                         bwat[i,j],bmr[i,j],till_pw_fraction, 
                                         bwat_max)
          N = p_over - p_w #  effective pressure on till
          if N<Nmin:
            Nmin = N
          tauc[i,j] = till_c_0 + N * math.tan((math.pi/180.0) * tillphi[i,j])

  def updateTillPhi_algebraic(self,mask,thickness,bwat,bmr,tauc,tillphi,tillphi_prev=None):
    config = PISM.global_config()
    till_pw_fraction = self.till_pw_fraction#config.get("till_pw_fraction")
    till_c_0 = self.till_c_0#config.get("till_c_0") * 1e3 # convert from kPa to Pa
    bwat_max = self.bwat_max#config.get("bwat_max");

    rho_g = self.rho_g

    Nmin = 1. # This is dorky.

    vars = [mask,thickness,bwat,bmr,tauc]
    if not tillphi_prev is None:
      vars.append(tillphi_prev)
    with PISM.util.Access(nocomm=vars,comm=tillphi):
      mq = PISM.MaskQuery(mask)
      GHOSTS = self.grid.max_stencil_width;
      for (i,j) in self.grid.points_with_ghosts(nGhosts = GHOSTS):
        if mq.floating_ice(i,j):
          if not tillphi_prev is None:
            tillphi[i,j] = tillphi_prev[i,j]
          continue
          
        H_ij = thickness[i,j]
        if H_ij == 0:
          if not tillphi_prev is None:
            tillphi[i,j] = tillphi_prev[i,j]
        else: # grounded and there is some ice
          p_over = rho_g * H_ij
          p_w    = self.getBasalWaterPressure( H_ij,
                                         bwat[i,j],bmr[i,j],till_pw_fraction, 
                                         bwat_max)
          N = p_over - p_w #  effective pressure on till
          if abs(N) < Nmin:
            if not tillphi_prev is None:
              tillphi[i,j] = tillphi_prev[i,j]
          else:
            tillphi[i,j] = (180./math.pi)*math.atan((tauc[i,j]-till_c_0)/N)

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


context = PISM.Context()
config = context.config()

PISM.set_abort_on_sigint(True)

usage = \
"""  sia.py -i IN.nc [-o file.nc]
  where:
    -i      IN.nc is input file in NetCDF format: contains PISM-written model state
  notes:
    * -i is required
  """

PISM.show_usage_check_req_opts(context.com,"sia.py",["-i"],usage)

for o in PISM.OptionsGroup(context.com,"","tauc2tillphi"):
  bootfile = PISM.optionsString("-i","input file")
  output_file = PISM.optionsString("-o","output file",default="tauc2tillphi_"+os.path.basename(bootfile))

  verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)
  PISM.set_config_from_options(context.com,config)

grid = PISM.Context().newgrid()
PISM.util.init_grid_from_file(grid,bootfile,
                                periodicity=PISM.XY_PERIODIC);

enthalpyconverter = PISM.EnthalpyConverter(config)
if PISM.getVerbosityLevel() >3:
  enthalpyconverter.viewConstants(PETSc.Viewer.STDOUT())

if PISM.optionsIsSet("-ssa_glen"):
  ice = PISM.CustomGlenIce(com,"",config,enthalpyconverter)
  B_schoof = 3.7e8;     # Pa s^{1/3}; hardness 
  ice.setHardness(B_schoof)
else:
  ice =  PISM.GPBLDIce(grid.com, "", config,enthalpyconverter)
ice.setFromOptions()

surface    = PISM.util.standardIceSurfaceVec( grid )
thickness  = PISM.util.standardIceThicknessVec( grid )
bed        = PISM.util.standardBedrockElevationVec( grid )
enthalpy   = PISM.util.standardEnthalpyVec( grid )
ice_mask   = PISM.util.standardIceMask( grid )
v = [surface,thickness,bed,enthalpy,ice_mask]
for var in v:
  var.regrid(bootfile,True)
tauc       = PISM.util.standardYieldStressVec( grid )

# Compute yield stress from PISM state variables
# (basal melt rate, tillphi, and basal water height)
bmr   = PISM.util.standardBasalMeltRateVec(grid)
tillphi = PISM.util.standardTillPhiVec(grid)
bwat = PISM.util.standardBasalWaterVec(grid)
for v in [bmr,tillphi,bwat]:
  v.regrid(bootfile,True)

standard_gravity = config.get("standard_gravity")
ice_rho = ice.rho
basal_till = BasalTillStrength(grid,ice_rho,standard_gravity)

basal_till.updateYieldStress(ice_mask, thickness, bwat, bmr, tillphi,tauc)
tillphi2 = PISM.util.standardTillPhiVec(grid,name="tillphi_2")
tillphi2.set(0.)
basal_till.updateTillPhi_algebraic(ice_mask, thickness, bwat, bmr, tauc, tillphi2, tillphi_prev=tillphi)


pio = PISM.PISMIO(grid)
pio.open_for_writing(output_file,False,True)
pio.append_time(grid.config.get_string("time_dimension_name"),0.0)
pio.close()

# Save time & command line
PISM.util.writeProvenance(output_file)
tillphi.write(output_file)
tillphi2.write(output_file)


