#! /usr/bin/env python
#
# Copyright (C) 2012, 2014 David Maxwell
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

help_description = """Verifies that the end result of an SSA inversion using 'pismi.py'
terminates after an expected number of iterations and with an expected final misfit."""
help_epilog = """E.g.

verify_ssa_inv.py inversion_result.nc --desired_misfit 100 --iter_max 50 --morozov

will succeed if the inversion terminated with a misfit less than 100 and terminated using no more
than 50 iterations.  The --morozov flag indicates that the next-to-last iteration should
checked to have a misfit of at least 100.  On the other hand

verify_ssa_inv.py inversion_result.nc --desired_misfit 100 --iter_max 50 --misfit_tolerance 2.5

will succeed if the inversion terminated with a misfit within 2.5 of 100 and terminated using no more
than 50 iterations.

"""


def fail(msg):
    print("Test Failed: %s" % msg)
    exit(1)


def success():
    print("Test Succeeded.")
    exit(0)


netCDF = None
try:
    import netCDF4 as netCDF
except:
    fail('Unable to import netCDF4.')
    exit(1)

import argparse

parser = argparse.ArgumentParser(description=help_description,
                                 epilog=help_epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('inv_filename', help='NC file with inversion results', type=str)
parser.add_argument('-m', '--desired_misfit', required=True, type=float, help='desired final misfit')
parser.add_argument('-i', '--iter_max', required=True, type=int, help='maximum number of iterations')
parser.add_argument('-z', '--morozov', required=False, action='store_true',
                    help='verify that next-to-last misfit is larger than desired misfit')
parser.add_argument('-t', '--misfit_tolerance', required=False, type=float,
                    help='maximum difference between desired and computed misfit')

args = parser.parse_args()

inv_filename = args.inv_filename
try:
    data = netCDF.Dataset(inv_filename)
except:
    fail('Unable to open inversion file "%s"\nPerhaps the inversion run failed to converge, or terminated unexpectedly.' % inv_filename)

# Grab the misfit history from the file.
misfit = data.variables.get('inv_ssa_misfit')
if misfit is None:
    fail('Inversion file %s missing misfit history variable "inv_ssa_misfit".\nPerhaps the inversion run failed to converge, or terminated unexpectedly.' % inv_filename)

desired_misfit = args.desired_misfit

iter_count = len(misfit)
last_misfit = misfit[-1]

if iter_count > args.iter_max:
    fail("Inversion took an excessive number of iterations: %d > %d" % (iter_count, args.iter_max))

if args.misfit_tolerance:
    d_misfit = abs(last_misfit - desired_misfit)
    if d_misfit > args.misfit_tolerance:
        fail("Final misfit too far from desired misfit: |final-desired| = |%0.4g-%0.4g| = %0.4g > %0.4g = desired" %
             (last_misfit, desired_misfit, d_misfit, args.misfit_tolerance))
else:
    if last_misfit >= desired_misfit:
        fail("Desired final misfit not met: computed = %0.4g >= %0.4g = desired" % (last_misfit, desired_misfit))

if args.morozov:
    penultimate_misfit = misfit[-2]
    if penultimate_misfit < desired_misfit:
        fail("Penultimate misfit too small: computed = %0.4g < %0.4g = desired" % (penultimate_misfit, desired_misfit))


success()
