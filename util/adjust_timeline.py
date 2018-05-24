#!/usr/bin/env python

# @package adjust_timeline
# \author Andy Aschwanden, University of Alaska Fairbanks, USA
# \brief Script adjusts a time axis of a file.
# \details Script adjusts the time axis of a file.

# Say you have monthly climate forcing from 1980-1-1 through 2001-1-1 in
# the forcing file foo_1980-1999.nc to be used with, e.g. -surface_given_file,
# but you want the model to run from 1991-1-1 through 2001-1-1.
#
# Usage:
#
# \verbatim $ adjust_timeline.py --start_date '1991-1-1'
# time_1991-2000.nc \endverbatim

import os
from argparse import ArgumentParser
from dateutil import rrule
from dateutil.parser import parse
import time
import numpy as np

try:
    import netCDF4 as netCDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)
NC = netCDF.Dataset
from netcdftime import utime, datetime

# Set up the option parser
parser = ArgumentParser()
parser.description = '''Script adjusts the time file with time and time
bounds that can be used to determine to force PISM via command line
option -time_file or adjust the time axis for postprocessing.'''
parser.add_argument("FILE", nargs='*')
parser.add_argument("-p", "--periodicity", dest="periodicity",
                    help='''periodicity, e.g. monthly, daily, etc. Default=monthly''',
                    default="monthly")
parser.add_argument("-a", "--start_date", dest="start_date",
                    help='''Start date in ISO format. Default=1989-1-1''',
                    default='1989-1-1')
parser.add_argument("-c", "--calendar", dest="calendar",
                    choices=['standard', 'gregorian', 'no_leap', '365_day', '360_day', 'julian'],
                    help='''Sets the calendar. Default="standard".''',
                    default='standard')
parser.add_argument("-i", "--interval_type", dest="interval_type",
                    choices=['start', 'mid', 'end'],
                    help='''Defines whether the time values t_k are the end points of the time bounds tb_k or the mid points 1/2*(tb_k -tb_(k-1)). Default="mid".''',
                    default='mid')
parser.add_argument("-u", "--ref_unit", dest="ref_unit",
                    help='''Reference unit. Default=days. Use of months or
                    years is NOT recommended.''', default='days')
parser.add_argument("-d", "--ref_date", dest="ref_date",
                    help='''Reference date. Default=1960-1-1''',
                    default='1960-1-1')

options = parser.parse_args()
interval_type = options.interval_type
periodicity = options.periodicity.upper()
start_date = parse(options.start_date)
ref_unit = options.ref_unit
ref_date = options.ref_date
args = options.FILE
infile = args[0]

time1 = time.time()


nc = NC(infile, 'a')
nt = len(nc.variables['time'])

time_units = ("%s since %s" % (ref_unit, ref_date))
time_calendar = options.calendar

cdftime = utime(time_units, time_calendar)

# create a dictionary so that we can supply the periodicity as a
# command-line argument.
pdict = {}
pdict['SECONDLY'] = rrule.SECONDLY
pdict['MINUTELY'] = rrule.MINUTELY
pdict['HOURLY'] = rrule.HOURLY
pdict['DAILY'] = rrule.DAILY
pdict['WEEKLY'] = rrule.WEEKLY
pdict['MONTHLY'] = rrule.MONTHLY
pdict['YEARLY'] = rrule.YEARLY
prule = pdict[periodicity]

# reference date from command-line argument
r = time_units.split(' ')[2].split('-')
refdate = datetime(int(r[0]), int(r[1]), int(r[2]))

# create list with dates from start_date for nt counts
# periodicity prule.
bnds_datelist = list(rrule.rrule(prule, dtstart=start_date, count=nt+1))

# calculate the days since refdate, including refdate, with time being the
bnds_interval_since_refdate = cdftime.date2num(bnds_datelist)
if interval_type == 'mid':
    # mid-point value:
    # time[n] = (bnds[n] + bnds[n+1]) / 2
    time_interval_since_refdate = (bnds_interval_since_refdate[0:-1] +
                                   np.diff(bnds_interval_since_refdate) / 2)
elif interval_type == 'start':
    time_interval_since_refdate = bnds_interval_since_refdate[:-1]
else:
    time_interval_since_refdate = bnds_interval_since_refdate[1:]

# create a new dimension for bounds only if it does not yet exist
time_dim = "time"
if time_dim not in list(nc.dimensions.keys()):
    nc.createDimension(time_dim)

# create a new dimension for bounds only if it does not yet exist
bnds_dim = "nb2"
if bnds_dim not in list(nc.dimensions.keys()):
    nc.createDimension(bnds_dim, 2)

# variable names consistent with PISM
time_var_name = "time"
bnds_var_name = "time_bnds"

# create time variable
if time_var_name not in nc.variables:
    time_var = nc.createVariable(time_var_name, 'd', dimensions=(time_dim))
else:
    time_var = nc.variables[time_var_name]
time_var[:] = time_interval_since_refdate
time_var.bounds = bnds_var_name
time_var.units = time_units
time_var.calendar = time_calendar
time_var.standard_name = time_var_name
time_var.axis = "T"

# create time bounds variable
if bnds_var_name not in nc.variables:
    time_bnds_var = nc.createVariable(bnds_var_name, 'd', dimensions=(time_dim, bnds_dim))
else:
    time_bnds_var = nc.variables[bnds_var_name]
time_bnds_var[:, 0] = bnds_interval_since_refdate[0:-1]
time_bnds_var[:, 1] = bnds_interval_since_refdate[1::]

# writing global attributes
script_command = ' '.join([time.ctime(), ':', __file__.split('/')[-1],
                           ' '.join([str(x) for x in args])])
nc.history = script_command
nc.Conventions = "CF 1.5"
nc.close()

time2 = time.time()
print('adjust_timeline.py took {}s'.format(time2 - time1))
