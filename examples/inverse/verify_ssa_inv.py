#! /usr/bin/env python
#
# Copyright (C) 2012 David Maxwell
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

help_description = """Verifies that the end result of an SSA inversion using 'vel2tauc.py'
terminates after an expected number of iterations and with an expected final misfit."""
help_epilog = """E.g.

verify_ssa_inv.py inversion_result.nc --desired_misfit 100 --desired_misfit_tolerance 5 --iter_max 50

will succeed if the inversion terminated with a misfit of 100 +/- 5 m/a and terminated using no more
than 50 iterations."""

def fail(msg):
  print "Test Failed: %s" % msg
  exit(1)

def success():
  print "Test Succeeded."
  exit(0)

netCDF = None
try:
  import netCDF3 as netCDF
except:
  try:
    import netCDF4 as netCDF
  except:
    fail('Unable to import netCDF3/netCDF4.')
    exit(1)

import argparse

parser = argparse.ArgumentParser(description=help_description,
epilog=help_epilog,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-m','--desired_misfit', required=True,type=float,help='desired final misfit')
parser.add_argument('-e','--desired_misfit_tolerance',required=True, type=float,help='acceptable margin of error for final misfit')
parser.add_argument('-i','--iter_max', required=True, type=int,help='maximum number of iterations')
args = parser.parse_args()

inv_filename = 'tiny_inv.nc'
try:
  data = netCDF.Dataset(inv_filename)
except:
  fail('Unable to open inversion file "%s"\nPerhaps the inversion run failed to converge, or terminated unexpectedly.' % inv_filename)

# Grab the misfit history from the file.
misfit = data.variables.get('inv_ssa_misfit')
if misfit is None:
  fail('Inversion file %s missing misfit history variable "inv_ssa_misfit".\nPerhaps the inversion run failed to converge, or terminated unexpectedly.' % inv_filename)

desired_misfit = args.desired_misfit
eps_misfit = args.desired_misfit_tolerance

iter_count = len(misfit)
last_misfit = misfit[-1]

if iter_count > args.iter_max :
  fail("Inversion took an excessive number of iterations: %d > %d" % (iter_count,args.iter_max))

if abs(last_misfit-desired_misfit) > eps_misfit:
  fail("Final misfit is unexpected: computed %0.4g and desired %0.4g +/- %0.4g" % (last_misfit,desired_misfit,eps_misfit))

success()
