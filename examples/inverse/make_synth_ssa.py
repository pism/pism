#! /usr/bin/env python
#
# Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2018 David Maxwell and Constantine Khroulev
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


import PISM
from PISM.util import convert
import sys, time, math

design_prior_scale = 0.2
design_prior_const = None


def groundedIceMisfitWeight(modeldata):
    grid = modeldata.grid
    weight = PISM.model.createVelocityMisfitWeightVec(grid)
    mask = modeldata.vecs.mask
    with PISM.vec.Access(comm=weight, nocomm=mask):
        weight.set(0.)
        grounded = PISM.MASK_GROUNDED
        for (i, j) in grid.points():
            if mask[i, j] == grounded:
                weight[i, j] = 1
    return weight


def fastIceMisfitWeight(modeldata, vel_ssa, threshold):
    grid = modeldata.grid
    weight = PISM.model.createVelocityMisfitWeightVec(grid)
    mask = modeldata.vecs.ice_mask
    threshold = threshold * threshold
    with PISM.vec.Access(comm=weight, nocomm=[vel_ssa, mask]):
        weight.set(0.)
        grounded = PISM.MASK_GROUNDED
        for (i, j) in grid.points():
            u = vel_ssa[i, j].u
            v = vel_ssa[i, j].v
            if mask[i, j] == grounded:
                if u * u + v * v > threshold:
                    weight[i, j] = 1
    return weight


# The main code for a run follows:
if __name__ == '__main__':
    context = PISM.Context()
    com = context.com

    PISM.set_abort_on_sigint(True)

    PISM.verbPrintf(2, PISM.Context().com, "SSA forward model.\n")
    if PISM.OptionBool("-version", "stop after printing PISM version"):
        sys.exit(0)

    usage = \
        """  %s -i IN.nc -Mx number -My number [-o file.nc]
  or (at python prompt)
    run %s -i IN.nc -Mx number -My number [-o file.nc]
  where:
    -i      IN.nc is input file in NetCDF format: contains PISM-written model state
    -Mx     number of grid points in the x direction
    -My     number of grid points in the y direction
  notes:
    * -i is required
  """ % (sys.argv[0], sys.argv[0])

    PISM.show_usage_check_req_opts(context.log, sys.argv[0], ["-i"], usage)

    config = context.config
    if not PISM.OptionString("-ssa_method", "").is_set():
        config.set_string("stress_balance.ssa.method", "fem")

    input_file_name = config.get_string("input.file")

    config.set_string("output.file_name", "make_synth_ssa.nc", PISM.CONFIG_DEFAULT)
    output_file_name = config.get_string("output.file_name")

    design_prior_scale = PISM.OptionReal("-design_prior_scale",
                                          "initial guess for design variable to be this factor of the true value",
                                          design_prior_scale)

    design_prior_const = PISM.OptionReal("-design_prior_const",
                                          "initial guess for design variable to be this constant",
                                          0.0)
    design_prior_const = design_prior_const.value() if design_prior_const.is_set() else None

    noise = PISM.OptionReal("-rms_noise", "pointwise rms noise to add (in m/a)", 0.0)
    noise = noise.value() if noise.is_set() else None

    misfit_weight_type = PISM.OptionKeyword("-misfit_type",
                                            "Choice of misfit weight function",
                                            "grounded,fast",
                                            "grounded").value()

    fast_ice_speed = PISM.OptionReal("-fast_ice_speed",
                                      "Threshold in m/a for determining if ice is fast",
                                      500.0)

    generate_ssa_observed = PISM.OptionBool("-generate_ssa_observed",
                                             "generate observed SSA velocities")

    is_regional = PISM.OptionBool("-regional",
                                   "Compute SIA/SSA using regional model semantics")

    design_var = PISM.OptionKeyword("-inv_ssa",
                                    "design variable for inversion",
                                    "tauc,hardav",
                                    "tauc").value()

    ssa_run = PISM.ssa.SSAFromInputFile(input_file_name)

    ssa_run.setup()

    modeldata = ssa_run.modeldata
    grid = modeldata.grid
    vecs = modeldata.vecs

    # add everything we need to "vecs" *before* it is locked in ssa_run.setup()
    if design_var == 'tauc':
        # Generate a prior guess for tauc
        tauc_prior = PISM.model.createYieldStressVec(grid, name='tauc_prior',
                                                     desc="initial guess for (pseudo-plastic)"
                                                     " basal yield stress in an inversion")
        vecs.add(tauc_prior, writing=True)
    elif design_var == 'hardav':
        # Generate a prior guess for hardav
        vecs.add(PISM.model.createAveragedHardnessVec(grid))

        hardav_prior = PISM.model.createAveragedHardnessVec(grid,
                                                            name='hardav_prior',
                                                            desc="initial guess for vertically averaged"
                                                            " ice hardness in an inversion")
        vecs.add(hardav_prior, writing=True)

    solve_t0 = time.clock()
    vel_ssa = ssa_run.solve()
    solve_t = time.clock() - solve_t0

    PISM.verbPrintf(2, context.com, "Solve time %g seconds.\n", solve_t)

    if design_var == 'tauc':
        if design_prior_const is not None:
            vecs.tauc_prior.set(design_prior_const)
        else:
            vecs.tauc_prior.copy_from(modeldata.vecs.tauc)
            vecs.tauc_prior.scale(design_prior_scale)

        tauc_true = modeldata.vecs.tauc
        tauc_true.metadata(0).set_name('tauc_true')
        tauc_true.set_attrs("diagnostic",
                            "value of basal yield stress used to generate synthetic SSA velocities",
                            "Pa", "")
        vecs.markForWriting(tauc_true)
    elif design_var == 'hardav':
        # Generate a prior guess for hardav

        EC = PISM.EnthalpyConverter(config)
        ice_factory = PISM.IceFlowLawFactory(grid.com, "stress_balance.ssa.", config, EC)
        ice_factory.removeType(PISM.ICE_GOLDSBY_KOHLSTEDT)
        ice_factory.setType(config.get_string("stress_balance.ssa.flow_law"))
        ice_factory.setFromOptions()
        flow_law = ice_factory.create()
        averaged_hardness_vec(flow_law, vecs.land_ice_thickness, vecs.enthalpy, vecs.hardav)

        if design_prior_const is not None:
            vecs.hardav_prior.set(design_prior_const)
        else:
            vecs.hardav_prior.copy_from(vecs.hardav)
            vecs.hardav_prior.scale(hardav_prior_scale)

        hardav_true = vecs.hardav
        hardav_true.metadata(0).set_name('hardav_true')
        hardav_true.set_attrs("diagnostic",
                              "vertically averaged ice hardness used to generate synthetic SSA velocities",
                              "Pa s^0.33333", "")
        vecs.markForWriting(hardav_true)

    vel_ssa_observed = vel_ssa    # vel_ssa = ssa_run.solve() earlier

    vel_ssa_observed.metadata(0).set_name("u_ssa_observed")
    vel_ssa_observed.metadata(0).set_string("long_name", "x-component of 'observed' SSA velocities")

    vel_ssa_observed.metadata(1).set_name("v_ssa_observed")
    vel_ssa_observed.metadata(1).set_string("long_name", "y-component of 'observed' SSA velocities")

    if generate_ssa_observed:
        vecs.markForWriting(vel_ssa_observed)
        final_velocity = vel_ssa_observed
    else:
        sia_solver = PISM.SIAFD
        if is_regional:
            sia_solver = PISM.SIAFD_Regional
        vel_sia_observed = PISM.sia.computeSIASurfaceVelocities(modeldata, sia_solver)

        vel_sia_observed.metadata(0).set_name('u_sia_observed')
        vel_sia_observed.metadata(0).set_string('long_name', "x-component of the 'observed' SIA velocities")

        vel_sia_observed.metadata(1).set_name('v_sia_observed')
        vel_sia_observed.metadata(1).set_string('long_name', "y-component of the 'observed' SIA velocities")

        vel_surface_observed = PISM.model.create2dVelocityVec(grid, "_surface_observed",
                                                              "observed surface velocities",
                                                              stencil_width=1)
        vel_surface_observed.copy_from(vel_sia_observed)
        vel_surface_observed.add(1., vel_ssa_observed)
        vecs.markForWriting(vel_surface_observed)
        final_velocity = vel_surface_observed

    # Add the misfit weight.
    if misfit_weight_type == "fast":
        misfit_weight = fastIceMisfitWeight(modeldata, vel_ssa,
                                            convert(fast_ice_speed, "m/year", "m/second"))
    else:
        misfit_weight = groundedIceMisfitWeight(modeldata)
    modeldata.vecs.add(misfit_weight, writing=True)

    if not noise is None:
        u_noise = PISM.vec.randVectorV(grid, noise / math.sqrt(2), final_velocity.stencil_width())
        final_velocity.add(convert(1.0, "m/year", "m/second"),
                           u_noise)

    pio = PISM.util.prepare_output(output_file_name)
    pio.close()

    vecs.write(output_file_name)

    # Save time & command line
    PISM.util.writeProvenance(output_file_name)
