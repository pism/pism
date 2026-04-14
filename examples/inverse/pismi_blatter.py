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

"""Blatter-based inversion for basal yield stress (tauc).

Usage:
  mpirun -np N python pismi_blatter.py -i STATE.nc -inv_data OBS.nc -o OUT.nc \\
      -inv_method tikhonov_lmvm -inv_target_misfit 50 [-blatter_Mz 11]

Analogous to pismi.py but uses the Blatter (higher-order) stress balance
instead of SSA. The observed velocities are matched against 2D surface
velocities extracted from the 3D Blatter solution.
"""

import PISM
import PISM.invert.blatter as blinv
from PISM.invert.ssa import createTikhonovFunctionals, printIteration, \
    printTikhonovProgress, MisfitLogger, ZetaSaver
from PISM.logging import logMessage
from PISM.util import convert
import sys
import os


def adjustTauc(mask, tauc):
    """Set tauc to zero where floating and to a high value where ice-free."""
    logMessage("  Adjusting initial estimate of 'tauc' to match PISM model"
               " for floating ice and ice-free bedrock.\n")
    grid = mask.grid()
    high_tauc = grid.ctx().config().get_number("basal_yield_stress.ice_free_bedrock")
    with PISM.vec.Access(comm=tauc, nocomm=mask):
        for (i, j) in grid.points():
            if mask.ocean(i, j):
                tauc[i, j] = 0
            elif mask.ice_free(i, j):
                tauc[i, j] = high_tauc


def createDesignVec(grid, design_var, name=None, **kwargs):
    if name is None:
        name = design_var
    if design_var == "tauc":
        return PISM.model.createYieldStressVec(grid, name=name, **kwargs)
    elif design_var == "hardav":
        return PISM.model.createAveragedHardnessVec(grid, name=name, **kwargs)
    else:
        raise ValueError("Unknown design variable %s" % design_var)


def run():
    context = PISM.Context()
    config = context.config
    com = context.com
    PISM.set_abort_on_sigint(True)

    # ── Parse command-line options ─────────────────────────────────────

    input_filename = config.get_string("input.file")
    if len(input_filename) == 0:
        input_filename = None

    append = PISM.OptionString("-a", "append file")
    append_filename = append.value() if append.is_set() else None

    output_filename = config.get_string("output.file")
    if len(output_filename) == 0:
        output_filename = None

    append_mode = False

    if (input_filename is None) and (append_filename is None):
        PISM.verbPrintf(1, com,
                        "\nError: No input file specified."
                        " Use -i [file.nc] or -a [file.nc].\n")
        sys.exit(1)

    if (input_filename is not None) and (append_filename is not None):
        PISM.verbPrintf(1, com, "\nError: Only one of -i/-a is allowed.\n")
        sys.exit(1)

    if (output_filename is not None) and (append_filename is not None):
        PISM.verbPrintf(1, com, "\nError: Only one of -a/-o is allowed.\n")
        sys.exit(1)

    if append_filename is not None:
        input_filename = append_filename
        output_filename = append_filename
        append_mode = True

    inv_data_filename = PISM.OptionString(
        "-inv_data", "inverse data file", input_filename).value()

    design_var = PISM.OptionKeyword(
        "-inv_blatter", "design variable for inversion",
        "tauc,hardav", "tauc").value()

    do_restart = PISM.OptionBool("-inv_restart",
                                  "Restart a stopped computation.")
    use_design_prior = config.get_flag("inverse.use_design_prior")
    using_zeta_fixed_mask = config.get_flag("inverse.use_zeta_fixed_mask")
    inv_method = config.get_string("inverse.stress_balance.method")

    # Blatter-specific options
    blatter_Mz = 11
    Mz_opt = PISM.OptionInteger("-blatter_Mz",
                                "Number of Blatter vertical levels", blatter_Mz)
    if Mz_opt.is_set():
        blatter_Mz = int(Mz_opt.value())

    coarsening = 2
    coarsening_opt = PISM.OptionInteger(
        "-blatter_coarsening",
        "Blatter multigrid coarsening factor", coarsening)
    if coarsening_opt.is_set():
        coarsening = int(coarsening_opt.value())

    if output_filename is None:
        output_filename = "pismi_blatter_" + os.path.basename(input_filename)

    saving_inv_data = (inv_data_filename != output_filename)

    # ── Set up the Blatter forward problem ────────────────────────────

    PISM.verbPrintf(1, com, "Blatter inverse model.\n")

    forward_run = blinv.BlatterForwardRunFromInputFile(
        input_filename, inv_data_filename, design_var,
        Mz=blatter_Mz, coarsening_factor=coarsening)
    forward_run.setup()

    design_param = forward_run.designVariableParameterization()
    solver = blinv.createInvBlatterSolver(forward_run)

    modeldata = forward_run.modeldata
    vecs = modeldata.vecs
    grid = modeldata.grid

    # ── Design variable prior ─────────────────────────────────────────

    design_prior = createDesignVec(grid, design_var, '%s_prior' % design_var)
    long_name = design_prior.metadata().get_string("long_name")
    units = design_prior.metadata().get_string("units")
    design_prior.metadata().long_name(
        "best prior estimate for %s (used for inversion)" % long_name).units(units)

    if (PISM.util.fileHasVariable(inv_data_filename, "%s_prior" % design_var)
            and use_design_prior):
        logMessage("  Reading '%s_prior' from inverse data file %s.\n"
                   % (design_var, inv_data_filename))
        design_prior.regrid(inv_data_filename, critical=True)
        vecs.add(design_prior, writing=saving_inv_data)
    else:
        if not PISM.util.fileHasVariable(input_filename, design_var):
            PISM.verbPrintf(1, com,
                            "Initial guess for design variable is not available"
                            " as '%s' in %s.\n" % (design_var, input_filename))
            exit(1)
        logMessage("Reading '%s_prior' from '%s' in input file.\n"
                   % (design_var, design_var))
        design = createDesignVec(grid, design_var)
        design.regrid(input_filename, True)
        design_prior.copy_from(design)
        vecs.add(design_prior, writing=True)

    # ── Fixed tauc mask ───────────────────────────────────────────────

    if using_zeta_fixed_mask:
        if PISM.util.fileHasVariable(inv_data_filename, "zeta_fixed_mask"):
            zeta_fixed_mask = PISM.model.createZetaFixedMaskVec(grid)
            zeta_fixed_mask.regrid(inv_data_filename)
            vecs.add(zeta_fixed_mask)
        else:
            if design_var == 'tauc':
                logMessage("  Computing 'zeta_fixed_mask'.\n")
                zeta_fixed_mask = PISM.model.createZetaFixedMaskVec(grid)
                zeta_fixed_mask.set(1)
                mask = vecs.mask
                with PISM.vec.Access(comm=zeta_fixed_mask, nocomm=mask):
                    for (i, j) in grid.points():
                        if mask.grounded_ice(i, j):
                            zeta_fixed_mask[i, j] = 0
                vecs.add(zeta_fixed_mask)
                adjustTauc(vecs.mask, design_prior)
            elif design_var == 'hardav':
                pass
            else:
                raise NotImplementedError(
                    "Unable to build 'zeta_fixed_mask' for %s." % design_var)

    # ── Convert prior -> zeta ─────────────────────────────────────────

    zeta_prior = PISM.Scalar2(grid, "zeta_prior")
    design_param.convertFromDesignVariable(design_prior, zeta_prior)
    vecs.add(zeta_prior, writing=True)

    # ── Initial guess for zeta ────────────────────────────────────────

    zeta = PISM.Scalar2(grid, "zeta_inv")
    zeta.metadata(0).long_name("zeta_inv").units("1").output_units("1")

    if do_restart:
        if not PISM.util.fileHasVariable(output_filename, 'zeta_inv'):
            PISM.verbPrintf(1, com,
                            "Unable to restart: file %s missing 'zeta_inv'\n"
                            % output_filename)
            exit(1)
        logMessage("  Restarting from 'zeta_inv' in %s\n" % output_filename)
        zeta.regrid(output_filename, True)
    elif PISM.util.fileHasVariable(inv_data_filename, 'zeta_inv'):
        logMessage("  Starting from 'zeta_inv' in %s\n" % inv_data_filename)
        zeta.regrid(inv_data_filename, True)
    else:
        zeta.copy_from(zeta_prior)

    # ── Load observed velocities ──────────────────────────────────────

    vel_observed = PISM.model.create2dVelocityVec(
        grid, '_observed', stencil_width=2)

    if PISM.util.fileHasVariable(inv_data_filename, "u_observed"):
        vel_observed.regrid(inv_data_filename, True)
        vecs.add(vel_observed, writing=saving_inv_data)
    elif PISM.util.fileHasVariable(inv_data_filename, "u_ssa_observed"):
        vel_observed.metadata(0).set_name("u_ssa_observed")
        vel_observed.metadata(1).set_name("v_ssa_observed")
        vel_observed.regrid(inv_data_filename, True)
        vel_observed.metadata(0).set_name("u_observed")
        vel_observed.metadata(1).set_name("v_observed")
        vecs.add(vel_observed, writing=saving_inv_data)
    else:
        PISM.verbPrintf(1, com,
                        "No observed velocities (u/v_observed or"
                        " u/v_ssa_observed) found in %s.\n"
                        % inv_data_filename)
        exit(1)

    # ── Load true design variable if available (synthetic tests) ──────

    if PISM.util.fileHasVariable(inv_data_filename, "%s_true" % design_var):
        design_true = createDesignVec(grid, design_var,
                                       '%s_true' % design_var)
        design_true.regrid(inv_data_filename, True)
        vecs.add(design_true, writing=saving_inv_data)

    # ── Prepare output file ───────────────────────────────────────────

    message_logger = PISM.logging.CaptureLogger(output_filename,
                                                 'pismi_blatter_log')
    PISM.logging.add_logger(message_logger)
    if append_mode or do_restart:
        message_logger.readOldLog()

    if not append_mode:
        F = PISM.util.prepare_output(output_filename)
        F.close()
    zeta.write(output_filename)
    PISM.util.writeProvenance(output_filename)

    # ── Attach iteration listeners ────────────────────────────────────

    solver.addIterationListener(printIteration)

    misfit_logger = MisfitLogger()
    solver.addIterationListener(misfit_logger)

    if inv_method.startswith('tikhonov'):
        solver.addIterationListener(printTikhonovProgress)

    solver.addDesignUpdateListener(ZetaSaver(output_filename))

    # ── Run the inverse solver ────────────────────────────────────────

    if do_restart:
        logMessage('************** Restarting Blatter inversion. ****************\n')
    else:
        logMessage('============== Starting Blatter inversion. ==================\n')

    reason = solver.solveInverse(zeta_prior, vel_observed, zeta)
    if reason.failed():
        PISM.logging.logError("Inverse solve ended: %s\n"
                              "Writing partial results.\n"
                              % reason.nested_description(1))
    else:
        logMessage("Inverse solve success (%s)!\n" % reason.description())

    (zeta, u) = solver.inverseSolution()

    # If the state solution is empty (e.g. after DIVERGED_MAXITS),
    # do a final forward solve with the current zeta to get velocities.
    if u is not None and u.norm(PISM.PETSc.NormType.NORM_INFINITY)[0] == 0.0:
        logMessage("  State solution is zero; running final forward solve.\n")
        solver.solveForward(zeta, u)

    # ── Convert zeta back to design variable ──────────────────────────

    if vecs.has(design_var):
        design = vecs.get(design_var)
        design_param.convertToDesignVariable(zeta, design)
    else:
        design = createDesignVec(grid, design_var)
        design_param.convertToDesignVariable(zeta, design)
        vecs.add(design, writing=True)

    vecs.add(zeta, writing=True)

    u.metadata(0).set_name("u_inv")
    u.metadata(0).set_string("long_name",
                              "x-component of velocity computed by Blatter inversion")
    u.metadata(1).set_name("v_inv")
    u.metadata(1).set_string("long_name",
                              "y-component of velocity computed by Blatter inversion")
    vecs.add(u, writing=True)

    # ── Compute residual ──────────────────────────────────────────────

    residual = PISM.model.create2dVelocityVec(grid, name='_inv_residual')
    residual.copy_from(u)
    residual.add(-1, vel_observed)

    r_mag = PISM.Scalar(grid, "inv_residual")
    r_mag.metadata(0).long_name(
        "magnitude of mismatch between observed and inverted velocities"
    ).units("m s^-1").output_units("m year^-1")
    r_mag.metadata().set_number("_FillValue", convert(-0.01, 'm/year', 'm/s'))
    r_mag.metadata().set_number("valid_min", 0.0)

    PISM.compute_magnitude(residual, r_mag)
    PISM.apply_mask(vecs.land_ice_thickness, 0.0, r_mag)

    vecs.add(residual, writing=True)
    vecs.add(r_mag, writing=True)

    # ── Write output ──────────────────────────────────────────────────

    forward_run.write(output_filename, append=True)
    if not append_mode:
        message_logger.write(output_filename)
    misfit_logger.write(output_filename)


if __name__ == "__main__":
    run()
