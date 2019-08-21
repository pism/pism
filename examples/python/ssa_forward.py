#! /usr/bin/env python
#
# Copyright (C) 2011, 2014, 2015, 2016, 2018 David Maxwell and Constantine Khroulev
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
import time

# The main code for a run follows:
if __name__ == '__main__':
    context = PISM.Context()
    com = context.com

    PISM.set_abort_on_sigint(True)

    PISM.verbPrintf(2, PISM.Context().com, "SSA forward model.\n")
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

    PISM.show_usage_check_req_opts(context.log, "ssa_forward", ["-i"], usage)

    input_file = config.get_string("input.file")
    if len(input_file) == 0:
        import sys
        sys.exit(1)

    config.set_string("output.file_name", "ssa_forward.nc")

    ssa_run = PISM.ssa.SSAFromInputFile(input_file)

    ssa_run.setup()

    solve_t0 = time.clock()
    vel_ssa = ssa_run.solve()
    solve_t = time.clock() - solve_t0

    PISM.verbPrintf(2, context.com, "Solve time %g seconds.\n", solve_t)

    ssa_run.write(config.get_string("output.file_name"))
