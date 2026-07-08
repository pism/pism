#! /usr/bin/env python3
#
# Copyright (C) 2026 Andy Aschwanden
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

"""Generate synthetic observed surface velocities using the Blatter solver.

Analogous to make_synth_ssa.py but uses the Blatter (higher-order) stress
balance. The output file contains u_observed/v_observed (surface velocities),
tauc_prior, tauc_true, and vel_misfit_weight — everything needed to run
pismi_blatter.py.
"""

import PISM
import PISM.invert.blatter as blinv
from PISM.util import convert
import sys
import time
import math


design_prior_scale = 0.2
design_prior_const = None


def groundedIceMisfitWeight(modeldata):
    """Create a misfit weight that is 1 on grounded ice, 0 elsewhere."""
    grid = modeldata.grid
    weight = PISM.model.createVelocityMisfitWeightVec(grid)
    mask = modeldata.vecs.mask
    with PISM.vec.Access(comm=weight, nocomm=mask):
        weight.set(0.0)
        grounded = PISM.MASK_GROUNDED
        for (i, j) in grid.points():
            if mask[i, j] == grounded:
                weight[i, j] = 1
    return weight


if __name__ == '__main__':
    context = PISM.Context()
    com = context.com

    PISM.set_abort_on_sigint(True)

    PISM.verbPrintf(2, com, "Blatter forward model (synthetic observations).\n")
    if PISM.OptionBool("-version", "stop after printing PISM version"):
        sys.exit(0)

    config = context.config

    input_file_name = config.get_string("input.file")
    if len(input_file_name) == 0:
        PISM.verbPrintf(1, com, "Error: No input file. Use -i FILE.nc\n")
        sys.exit(1)

    config.set_string("output.file", "make_synth_blatter.nc", PISM.CONFIG_DEFAULT)
    output_file_name = config.get_string("output.file")

    sys_units = context.unit_system
    design_prior_scale = PISM.OptionReal(
        sys_units, "-design_prior_scale",
        "initial guess for design variable as this factor of the true value",
        "1", design_prior_scale).value()

    design_prior_const = PISM.OptionReal(
        sys_units, "-design_prior_const",
        "initial guess for design variable to be this constant",
        "m / s", 0.0)
    design_prior_const = design_prior_const.value() if design_prior_const.is_set() else None

    noise = PISM.OptionReal(
        sys_units, "-rms_noise",
        "pointwise rms noise to add (in m/a)",
        "m / year", 0.0)
    noise = noise.value() if noise.is_set() else None

    misfit_weight_type = PISM.OptionKeyword(
        "-misfit_type", "Choice of misfit weight function",
        "grounded", "grounded").value()

    # ── Set up the Blatter forward run ────────────────────────────────

    blatter_run = blinv.BlatterForwardRunFromInputFile(
        input_file_name, input_file_name, 'tauc')
    blatter_run.setup()

    modeldata = blatter_run.modeldata
    grid = modeldata.grid
    vecs = modeldata.vecs

    # Prior guess for tauc
    tauc_prior = PISM.model.createYieldStressVec(
        grid, name='tauc_prior',
        desc="initial guess for basal yield stress in an inversion")
    vecs.add(tauc_prior, writing=True)

    # ── Solve ─────────────────────────────────────────────────────────

    # Convert tauc -> zeta, then do a forward Blatter solve
    design_param = blatter_run.designVariableParameterization()
    zeta = PISM.Scalar2(grid, "zeta")
    design_param.convertFromDesignVariable(vecs.tauc, zeta)

    solve_t0 = time.time()
    reason = blatter_run.solver.linearize_at(zeta)
    solve_t = time.time() - solve_t0

    PISM.verbPrintf(2, com, "Blatter solve time %g seconds. %s\n",
                    solve_t, reason.description())

    if reason.failed():
        PISM.verbPrintf(1, com, "Blatter solve FAILED: %s\n",
                        reason.description())
        sys.exit(1)

    # ── Collect outputs ───────────────────────────────────────────────

    # Prior
    if design_prior_const is not None:
        vecs.tauc_prior.set(design_prior_const)
    else:
        vecs.tauc_prior.copy_from(vecs.tauc)
        vecs.tauc_prior.scale(design_prior_scale)

    # True tauc
    tauc_true = vecs.tauc
    tauc_true.metadata().set_name('tauc_true')
    tauc_true.metadata().long_name(
        "value of basal yield stress used to generate synthetic velocities"
    ).units("Pa")
    vecs.markForWriting(tauc_true)

    # Surface velocity from the Blatter solve
    vel_observed = blatter_run.solver.solution()

    vel_observed.metadata(0).set_name("u_observed")
    vel_observed.metadata(0).long_name(
        "x-component of 'observed' surface velocities from Blatter solver")

    vel_observed.metadata(1).set_name("v_observed")
    vel_observed.metadata(1).long_name(
        "y-component of 'observed' surface velocities from Blatter solver")

    vecs.add(vel_observed, writing=True)

    # Misfit weight
    misfit_weight = groundedIceMisfitWeight(modeldata)
    vecs.add(misfit_weight, writing=True)

    # Add noise if requested
    if noise is not None:
        u_noise = PISM.vec.randVectorV(
            grid, noise / math.sqrt(2), vel_observed.stencil_width())
        vel_observed.add(convert(1.0, "m/year", "m/second"), u_noise)

    # ── Write ─────────────────────────────────────────────────────────

    F = PISM.util.prepare_output(output_file_name)
    F.close()

    vecs.write(output_file_name)
    PISM.util.writeProvenance(output_file_name)
