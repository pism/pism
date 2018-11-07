#! /usr/bin/env python
#
# Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 David Maxwell and Constantine Khroulev
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

# try to start coverage
try:                            # pragma: no cover
    import coverage
    cov = coverage.coverage(branch=True)
    try:
        # try to load coverage data and ignore failures
        cov.load()
    except:
        pass
    cov.start()
except ImportError:             # pragma: no cover
    pass

import PISM
import PISM.invert.ssa
from PISM.logging import logMessage
from PISM.util import convert
import numpy as np
import sys, os, math


class SSAForwardRun(PISM.invert.ssa.SSAForwardRunFromInputFile):

    def write(self, filename, append=False):
        if not append:
            PISM.invert.ssa.SSAForwardRunFromInputFile.write(self, filename)
        else:
            grid = self.grid
            vecs = self.modeldata.vecs

            pio = PISM.PIO(grid.com, "netcdf3", filename, PISM.PISM_READWRITE)  # append mode!

            self.modeldata.vecs.write(filename)
            pio.close()


class InvSSAPlotListener(PISM.invert.listener.PlotListener):

    def __init__(self, grid, Vmax):
        PISM.invert.listener.PlotListener.__init__(self, grid)
        self.Vmax = Vmax
        self.l2_weight = None
        self.l2_weight_init = False

    def __call__(self, inverse_solver, count, data):

        if not self.l2_weight_init:
            vecs = inverse_solver.ssarun.modeldata.vecs
            if vecs.has('vel_misfit_weight'):
                self.l2_weight = self.toproczero(vecs.vel_misfit_weight)
            self.l2_weight_init = True

        method = inverse_solver.method

        r = self.toproczero(data.residual)
        Td = None
        if 'T_zeta_step' in data:
            Td = self.toproczero(data.T_zeta_step)
        TStarR = None
        if 'TStar_residual' in data:
            TStarR = self.toproczero(data.TStar_residual)
        d = None
        if 'zeta_step' in data:
            d = self.toproczero(data.zeta_step)
        zeta = self.toproczero(data.zeta)

        secpera = convert(1.0, "year", "second")

        if self.grid.rank() == 0:
            import matplotlib.pyplot as pp

            pp.figure(self.figure())

            l2_weight = self.l2_weight

            pp.clf()

            V = self.Vmax

            pp.subplot(2, 3, 1)
            if l2_weight is not None:
                rx = l2_weight * r[0, :, :] * secpera
            else:
                rx = r[0, :, :] * secpera
            rx = np.maximum(rx, -V)
            rx = np.minimum(rx, V)
            pp.imshow(rx, origin='lower', interpolation='nearest')
            pp.colorbar()
            pp.title('r_x')
            pp.jet()

            pp.subplot(2, 3, 4)
            if l2_weight is not None:
                ry = l2_weight * r[1, :, :] * secpera
            else:
                ry = r[1, :, :] * secpera
            ry = np.maximum(ry, -V)
            ry = np.minimum(ry, V)
            pp.imshow(ry, origin='lower', interpolation='nearest')
            pp.colorbar()
            pp.title('r_y')
            pp.jet()

            if method == 'ign':
                pp.subplot(2, 3, 2)
                Tdx = Td[0, :, :] * secpera
                pp.imshow(Tdx, origin='lower', interpolation='nearest')
                pp.colorbar()
                pp.title('Td_x')
                pp.jet()

                pp.subplot(2, 3, 5)
                Tdy = Td[1, :, :] * secpera
                pp.imshow(Tdy, origin='lower', interpolation='nearest')
                pp.colorbar()
                pp.title('Td_y')
                pp.jet()
            elif method == 'sd' or method == 'nlcg':
                pp.subplot(2, 3, 2)
                pp.imshow(TStarR, origin='lower', interpolation='nearest')
                pp.colorbar()
                pp.title('TStarR')
                pp.jet()

            if d is not None:
                d *= -1
                pp.subplot(2, 3, 3)
                pp.imshow(d, origin='lower', interpolation='nearest')

                # colorbar does a divide by zero if 'd' is all zero,
                # as it will be at the start of iteration zero.
                # The warning message is a distraction, so we suppress it.
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    pp.colorbar()
                    pp.jet()

                pp.title('-zeta_step')

            pp.subplot(2, 3, 6)
            pp.imshow(zeta, origin='lower', interpolation='nearest')
            pp.colorbar()
            pp.jet()
            pp.title('zeta')

            pp.ion()
            pp.draw()
            pp.show()


class InvSSALinPlotListener(PISM.invert.listener.PlotListener):

    def __init__(self, grid, Vmax):
        PISM.invert.listener.PlotListener.__init__(self, grid)
        self.Vmax = Vmax
        self.l2_weight = None
        self.l2_weight_init = False

    def __call__(self, inverse_solver, count, data):
        # On the first go-around, extract the l2_weight vector onto
        # processor zero.
        if self.l2_weight_init == False:
            vecs = inverse_solver.ssarun.modeldata.vecs
            self.l2_weight = self.toproczero(vecs.vel_misfit_weight)
            self.l2_init = True

        l2_weight = self.l2_weight
        r = self.toproczero(data.r)
        d = self.toproczero(data.d)

        if self.grid.rank() == 0:
            import matplotlib.pyplot as pp
            pp.figure(self.figure())
            pp.clf()

            V = self.Vmax
            pp.subplot(1, 3, 1)
            rx = l2_weight * r[0, :, :]
            rx = np.maximum(rx, -V)
            rx = np.minimum(rx, V)
            pp.imshow(rx, origin='lower', interpolation='nearest')
            pp.colorbar()
            pp.title('ru')
            pp.jet()

            pp.subplot(1, 3, 2)
            ry = l2_weight * r[1, :, :]
            ry = np.maximum(ry, -V)
            ry = np.minimum(ry, V)
            pp.imshow(ry, origin='lower', interpolation='nearest')
            pp.colorbar()
            pp.title('rv')
            pp.jet()

            d *= -1
            pp.subplot(1, 3, 3)
            pp.imshow(d, origin='lower', interpolation='nearest')
            pp.colorbar()
            pp.jet()
            pp.title('-d')

            pp.ion()
            pp.show()


def adjustTauc(mask, tauc):
    """Where ice is floating or land is ice-free, tauc should be adjusted to have some preset default values."""

    logMessage("  Adjusting initial estimate of 'tauc' to match PISM model for floating ice and ice-free bedrock.\n")

    grid = mask.grid()
    high_tauc = grid.ctx().config().get_double("basal_yield_stress.ice_free_bedrock")

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
        design_vec = PISM.model.createYieldStressVec(grid, name=name, **kwargs)
    elif design_var == "hardav":
        design_vec = PISM.model.createAveragedHardnessVec(grid, name=name, **kwargs)
    else:
        raise ValueError("Unknown design variable %s" % design_var)
    return design_vec


# Main code starts here
def run():
    context = PISM.Context()
    config = context.config
    com = context.com
    PISM.set_abort_on_sigint(True)

    WIDE_STENCIL = int(config.get_double("grid.max_stencil_width"))

    usage = \
        """  pismi.py [-i IN.nc [-o OUT.nc]]/[-a INOUT.nc] [-inv_data inv_data.nc] [-inv_forward model]
                [-inv_design design_var] [-inv_method meth]
    where:
    -i            IN.nc       is input file in NetCDF format: contains PISM-written model state
    -o            OUT.nc      is output file in NetCDF format to be overwritten
    -a            INOUT.nc    is input/output file in NetCDF format to be appended to
    -inv_data     inv_data.nc is data file containing extra inversion data (e.g. observed surface velocities)
    -inv_forward  model       forward model: only 'ssa' supported
    -inv_design   design_var  design variable name; one of 'tauc'/'hardav' for SSA inversions
    -inv_method   meth        algorithm for inversion [sd,nlcg,ign,tikhonov_lmvm]

    notes:
      * only one of -i/-a is allowed; both specify the input file
      * only one of -o/-a is allowed; both specify the output file
      * if -o is used, only the variables involved in inversion are written to the output file.
      * if -a is used, the varaibles involved in inversion are appended to the given file. No
        original variables in the file are changed.
   """

    append_mode = False

    input_filename = config.get_string("input.file")
    if len(input_filename) == 0:
        input_filename = None

    append = PISM.OptionString("-a", "append file")
    append_filename = append.value() if append.is_set() else None

    output_filename = config.get_string("output.file_name")
    if len(output_filename) == 0:
        output_filename = None

    if (input_filename is None) and (append_filename is None):
        PISM.verbPrintf(1, com, "\nError: No input file specified. Use one of -i [file.nc] or -a [file.nc].\n")
        sys.exit(0)

    if (input_filename is not None) and (append_filename is not None):
        PISM.verbPrintf(1, com, "\nError: Only one of -i/-a is allowed.\n")
        sys.exit(0)

    if (output_filename is not None) and (append_filename is not None):
        PISM.verbPrintf(1, com, "\nError: Only one of -a/-o is allowed.\n")
        sys.edit(0)

    if append_filename is not None:
        input_filename = append_filename
        output_filename = append_filename
        append_mode = True

    inv_data_filename = PISM.OptionString("-inv_data", "inverse data file", input_filename).value()

    do_plotting = PISM.OptionBool("-inv_plot", "perform visualization during the computation")
    do_final_plot = PISM.OptionBool("-inv_final_plot",
                                     "perform visualization at the end of the computation")
    Vmax = PISM.OptionReal("-inv_plot_vmax", "maximum velocity for plotting residuals", 30)

    design_var = PISM.OptionKeyword("-inv_ssa",
                                    "design variable for inversion",
                                    "tauc,hardav", "tauc").value()
    do_pause = PISM.OptionBool("-inv_pause", "pause each iteration")

    do_restart = PISM.OptionBool("-inv_restart", "Restart a stopped computation.")
    use_design_prior = config.get_boolean("inverse.use_design_prior")

    prep_module = PISM.OptionString("-inv_prep_module",
                                    "Python module used to do final setup of inverse solver")
    prep_module = prep_module.value() if prep_module.is_set() else None

    is_regional = PISM.OptionBool("-regional", "Compute SIA/SSA using regional model semantics")

    using_zeta_fixed_mask = config.get_boolean("inverse.use_zeta_fixed_mask")

    inv_method = config.get_string("inverse.ssa.method")

    if output_filename is None:
        output_filename = "pismi_" + os.path.basename(input_filename)

    saving_inv_data = (inv_data_filename != output_filename)

    forward_run = SSAForwardRun(input_filename, inv_data_filename, design_var)
    forward_run.setup()
    design_param = forward_run.designVariableParameterization()
    solver = PISM.invert.ssa.createInvSSASolver(forward_run)

    modeldata = forward_run.modeldata
    vecs = modeldata.vecs
    grid = modeldata.grid

    # Determine the prior guess for tauc/hardav. This can be one of
    # a) tauc/hardav from the input file (default)
    # b) tauc/hardav_prior from the inv_datafile if -inv_use_design_prior is set
    design_prior = createDesignVec(grid, design_var, '%s_prior' % design_var)
    long_name = design_prior.metadata().get_string("long_name")
    units = design_prior.metadata().get_string("units")
    design_prior.set_attrs("", "best prior estimate for %s (used for inversion)" % long_name, units, "")
    if PISM.util.fileHasVariable(inv_data_filename, "%s_prior" % design_var) and use_design_prior:
        PISM.logging.logMessage("  Reading '%s_prior' from inverse data file %s.\n" % (design_var, inv_data_filename))
        design_prior.regrid(inv_data_filename, critical=True)
        vecs.add(design_prior, writing=saving_inv_data)
    else:
        if not PISM.util.fileHasVariable(input_filename, design_var):
            PISM.verbPrintf(1, com, "Initial guess for design variable is not available as '%s' in %s.\nYou can provide an initial guess in the inverse data file.\n" % (
                design_var, input_filename))
            exit(1)
        PISM.logging.logMessage("Reading '%s_prior' from '%s' in input file.\n" % (design_var, design_var))
        design = createDesignVec(grid, design_var)
        design.regrid(input_filename, True)
        design_prior.copy_from(design)
        vecs.add(design_prior, writing=True)

    if using_zeta_fixed_mask:
        if PISM.util.fileHasVariable(inv_data_filename, "zeta_fixed_mask"):
            zeta_fixed_mask = PISM.model.createZetaFixedMaskVec(grid)
            zeta_fixed_mask.regrid(inv_data_filename)
            vecs.add(zeta_fixed_mask)
        else:
            if design_var == 'tauc':
                logMessage(
                    "  Computing 'zeta_fixed_mask' (i.e. locations where design variable '%s' has a fixed value).\n" % design_var)
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
                PISM.logging.logPrattle(
                    "Skipping 'zeta_fixed_mask' for design variable 'hardav'; no natural locations to fix its value.")
                pass
            else:
                raise NotImplementedError("Unable to build 'zeta_fixed_mask' for design variable %s.", design_var)

    # Convert design_prior -> zeta_prior
    zeta_prior = PISM.IceModelVec2S()
    zeta_prior.create(grid, "zeta_prior", PISM.WITH_GHOSTS, WIDE_STENCIL)
    design_param.convertFromDesignVariable(design_prior, zeta_prior)
    vecs.add(zeta_prior, writing=True)

    # Determine the initial guess for zeta.  If we are restarting, load it from
    # the output file.  Otherwise, if 'zeta_inv' is in the inverse data file, use it.
    # If none of the above, copy from 'zeta_prior'.
    zeta = PISM.IceModelVec2S()
    zeta.create(grid, "zeta_inv", PISM.WITH_GHOSTS, WIDE_STENCIL)
    zeta.set_attrs("diagnostic", "zeta_inv", "1", "zeta_inv")
    if do_restart:
        # Just to be sure, verify that we have a 'zeta_inv' in the output file.
        if not PISM.util.fileHasVariable(output_filename, 'zeta_inv'):
            PISM.verbPrintf(
                1, com, "Unable to restart computation: file %s is missing variable 'zeta_inv'", output_filename)
            exit(1)
        PISM.logging.logMessage("  Inversion starting from 'zeta_inv' found in %s\n" % output_filename)
        zeta.regrid(output_filename, True)

    elif PISM.util.fileHasVariable(inv_data_filename, 'zeta_inv'):
        PISM.logging.logMessage("  Inversion starting from 'zeta_inv' found in %s\n" % inv_data_filename)
        zeta.regrid(inv_data_filename, True)
    else:
        zeta.copy_from(zeta_prior)

    vel_ssa_observed = None
    vel_ssa_observed = PISM.model.create2dVelocityVec(grid, '_ssa_observed', stencil_width=2)
    if PISM.util.fileHasVariable(inv_data_filename, "u_ssa_observed"):
        vel_ssa_observed.regrid(inv_data_filename, True)
        vecs.add(vel_ssa_observed, writing=saving_inv_data)
    else:
        if not PISM.util.fileHasVariable(inv_data_filename, "u_surface_observed"):
            PISM.verbPrintf(
                1, context.com, "Neither u/v_ssa_observed nor u/v_surface_observed is available in %s.\nAt least one must be specified.\n" % inv_data_filename)
            exit(1)
        vel_surface_observed = PISM.model.create2dVelocityVec(grid, '_surface_observed', stencil_width=2)
        vel_surface_observed.regrid(inv_data_filename, True)
        vecs.add(vel_surface_observed, writing=saving_inv_data)

        sia_solver = PISM.SIAFD
        if is_regional:
            sia_solver = PISM.SIAFD_Regional
        vel_sia_observed = PISM.sia.computeSIASurfaceVelocities(modeldata, sia_solver)

        vel_sia_observed.metadata(0).set_name('u_sia_observed')
        vel_sia_observed.metadata(0).set_string('long_name', "x-component of the 'observed' SIA velocities")

        vel_sia_observed.metadata(1).set_name('v_sia_observed')
        vel_sia_observed.metadata(1).set_string('long_name', "y-component of the 'observed' SIA velocities")

        vel_ssa_observed.copy_from(vel_surface_observed)
        vel_ssa_observed.add(-1, vel_sia_observed)
        vecs.add(vel_ssa_observed, writing=True)

    # If the inverse data file has a variable tauc/hardav_true, this is probably
    # a synthetic inversion.  We'll load it now so that it will get written
    # out, if needed, at the end of the computation in the output file.
    if PISM.util.fileHasVariable(inv_data_filename, "%s_true" % design_var):
        design_true = createDesignVec(grid, design_var, '%s_true' % design_var)
        design_true.regrid(inv_data_filename, True)
        design_true.read_attributes(inv_data_filename)
        vecs.add(design_true, writing=saving_inv_data)

    # Establish a logger which will save logging messages to the output file.
    message_logger = PISM.logging.CaptureLogger(output_filename, 'pismi_log')
    PISM.logging.add_logger(message_logger)
    if append_mode or do_restart:
        message_logger.readOldLog()

    # Prep the output file from the grid so that we can save zeta to it during the runs.
    if not append_mode:
        pio = PISM.util.prepare_output(output_filename)
        pio.close()
    zeta.write(output_filename)

    # Log the command line to the output file now so that we have a record of
    # what was attempted
    PISM.util.writeProvenance(output_filename)

    # Attach various iteration listeners to the solver as needed for:

    # Iteration report.
    solver.addIterationListener(PISM.invert.ssa.printIteration)

    # Misfit reporting/logging.
    misfit_logger = PISM.invert.ssa.MisfitLogger()
    solver.addIterationListener(misfit_logger)

    if inv_method.startswith('tikhonov'):
        solver.addIterationListener(PISM.invert.ssa.printTikhonovProgress)

    # Saving the current iteration
    solver.addDesignUpdateListener(PISM.invert.ssa.ZetaSaver(output_filename))

    # Plotting
    if do_plotting:
        solver.addIterationListener(InvSSAPlotListener(grid, Vmax))
        if solver.method == 'ign':
            solver.addLinearIterationListener(InvSSALinPlotListener(grid, Vmax))

    # Solver is set up.  Give the user's prep module a chance to do any final
    # setup.

    if prep_module is not None:
        if prep_module.endswith(".py"):
            prep_module = prep_module[0:-2]
        exec("import %s as user_prep_module" % prep_module)
        user_prep_module.prep_solver(solver)

    # Pausing (add this after the user's listeners)
    if do_pause:
        solver.addIterationListener(PISM.invert.listener.pauseListener)

    # Run the inverse solver!
    if do_restart:
        PISM.logging.logMessage('************** Restarting inversion. ****************\n')
    else:
        PISM.logging.logMessage('============== Starting inversion. ==================\n')

    # Try solving
    reason = solver.solveInverse(zeta_prior, vel_ssa_observed, zeta)
    if reason.failed():
        PISM.logging.logError("Inverse solve FAILURE:\n%s\n" % reason.nested_description(1))
        quit()

    PISM.logging.logMessage("Inverse solve success (%s)!\n" % reason.description())

    (zeta, u) = solver.inverseSolution()

    # It may be that a 'tauc'/'hardav' was read in earlier.  We replace it with
    # our newly generated one.
    if vecs.has(design_var):
        design = vecs.get(design_var)
        design_param.convertToDesignVariable(zeta, design)
    else:
        # Convert back from zeta to tauc or hardav
        design = createDesignVec(grid, design_var)
        design_param.convertToDesignVariable(zeta, design)
        vecs.add(design, writing=True)

    vecs.add(zeta, writing=True)

    u.metadata(0).set_name("u_ssa_inv")
    u.metadata(0).set_string("long_name", "x-component of SSA velocity computed by inversion")

    u.metadata(1).set_name("v_ssa_inv")
    u.metadata(1).set_string("long_name", "y-component of SSA velocity computed by inversion")

    vecs.add(u, writing=True)

    residual = PISM.model.create2dVelocityVec(grid, name='_inv_ssa_residual')
    residual.copy_from(u)
    residual.add(-1, vel_ssa_observed)

    r_mag = PISM.IceModelVec2S()
    r_mag.create(grid, "inv_ssa_residual", PISM.WITHOUT_GHOSTS, 0)

    r_mag.set_attrs("diagnostic", "magnitude of mismatch between observed surface velocities and their reconstrution by inversion",
                    "m s-1", "inv_ssa_residual", 0)
    r_mag.metadata().set_double("_FillValue", convert(-0.01, 'm/year', 'm/s'))
    r_mag.metadata().set_double("valid_min", 0.0)
    r_mag.metadata().set_string("glaciological_units", "m year-1")

    r_mag.set_to_magnitude(residual)
    r_mag.mask_by(vecs.land_ice_thickness)

    vecs.add(residual, writing=True)
    vecs.add(r_mag, writing=True)

    # Write solution out to netcdf file
    forward_run.write(output_filename, append=append_mode)
    # If we're not in append mode, the previous command just nuked
    # the output file.  So we rewrite the siple log.
    if not append_mode:
        message_logger.write(output_filename)

    # Save the misfit history
    misfit_logger.write(output_filename)


if __name__ == "__main__":
    run()

# try to stop coverage and save a report:
try:                            # pragma: no cover
    cov.stop()
    report = PISM.OptionBool("-report_coverage", "save coverage information and a report")
    if report:
        cov.save()
        cov.html_report(include=["pismi.py"], directory="pismi_coverage")
except:                         # pragma: no cover
    pass
