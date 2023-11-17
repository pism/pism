#!/usr/bin/env python3

# @package create_timeline
# \author Andy Aschwanden, University of Alaska Fairbanks, USA
# \brief Script creates a timeline for enviornmental forcing.
# \details Script creates a timeline file that can be used with
# the -time_file option and -*_given_file to control the applicaton
# of forcing data, and to control start and end date of a PISM simulation.
#
# Say you have monthly climate forcing from 1980-1-1 through 2001-1-1 in
# the forcing file foo_1980-1999.nc to be used with, e.g. -surface_given_file,
# but you want the model to run from 1991-1-1 through 2001-1-1.
#
# Usage:
#
# \verbatim $ create_timeline.py --start_date '1991-1-1' -end_date '2001-1-1'
# time_1991-2000.nc \endverbatim

from argparse import ArgumentParser
from datetime import datetime

import cftime
import numpy as np
import xarray as xr
from dateutil import rrule
from dateutil.parser import parse

# Set up the option parser
parser = ArgumentParser()
parser.description = """Script creates a time file with time and time
bounds that can be used to determine to force PISM via command line
option -time_file"""
parser.add_argument("FILE", nargs="*")
parser.add_argument(
    "-c",
    "--calendar",
    dest="calendar",
    choices=["standard", "gregorian", "no_leap", "365_day", "360_day", "julian"],
    help="""Sets the calendar. Default="standard".""",
    default="standard",
)
parser.add_argument(
    "-p",
    "--periodicity",
    dest="periodicity",
    help="""periodicity, e.g. monthly, daily, etc. Default=monthly""",
    default="monthly",
)
parser.add_argument(
    "-a",
    "--start_date",
    dest="start_date",
    help="""Start date in ISO format. Default=1980-1-1""",
    default="1980-1-1",
)
parser.add_argument(
    "-e",
    "--end_date",
    dest="end_date",
    help="""End date in ISO format. Default=2020-1-1""",
    default="2020-1-1",
)
parser.add_argument(
    "-i",
    "--interval_type",
    dest="interval_type",
    choices=["start", "mid", "end"],
    help="""Defines whether the time values t_k are the end points of the time bounds tb_k or the mid points 1/2*(tb_k -tb_(k-1)). Default="mid".""",
    default="mid",
)
parser.add_argument(
    "-u",
    "--ref_unit",
    dest="ref_unit",
    help="""Reference unit. Default=seconds. Use of months or
                    years is NOT recommended.""",
    default="seconds",
)
parser.add_argument(
    "-d",
    "--ref_date",
    dest="ref_date",
    help="""Reference date. Default=1980-1-1""",
    default="1980-1-1",
)

options = parser.parse_args()
interval_type = options.interval_type
periodicity = options.periodicity.upper()
start_date = parse(options.start_date)
end_date = parse(options.end_date)
ref_unit = options.ref_unit
ref_date = options.ref_date
args = options.FILE
infile = args[0]

"""
Generate offset file using xarray
"""
time_units = f"{ref_unit} since {ref_date}"
calendar = options.calendar

# create a dictionary so that we can supply the periodicity as a
# command-line argument.
pdict = {}
pdict["SECONDLY"] = rrule.SECONDLY
pdict["MINUTELY"] = rrule.MINUTELY
pdict["HOURLY"] = rrule.HOURLY
pdict["DAILY"] = rrule.DAILY
pdict["WEEKLY"] = rrule.WEEKLY
pdict["MONTHLY"] = rrule.MONTHLY
pdict["YEARLY"] = rrule.YEARLY
prule = pdict[periodicity]

# reference date from command-line argument
r = time_units.split(" ")[2].split("-")
refdate = datetime(int(r[0]), int(r[1]), int(r[2]))

# create list with dates from start_date until end_date with
# periodicity prule.
bnds_datelist = list(rrule.rrule(prule, dtstart=start_date, until=end_date))

# calculate the days since refdate, including refdate, with time being the
bnds_interval_since_refdate = cftime.date2num(
    bnds_datelist, time_units, calendar=calendar
)
if interval_type == "mid":
    # mid-point value:
    # time[n] = (bnds[n] + bnds[n+1]) / 2
    time_interval_since_refdate = (
        bnds_interval_since_refdate[0:-1] + np.diff(bnds_interval_since_refdate) / 2
    )
elif interval_type == "start":
    time_interval_since_refdate = bnds_interval_since_refdate[:-1]
else:
    time_interval_since_refdate = bnds_interval_since_refdate[1:]

time_array = np.array(time_interval_since_refdate)
time_bounds = np.zeros((len(time_array), 2))
time_bounds[:, 0] = bnds_interval_since_refdate[0:-1]
time_bounds[:, 1] = bnds_interval_since_refdate[1::]

ds = xr.Dataset(
    data_vars=dict(  # pylint: disable=R1735
        time_bounds=(["time", "bnds"], time_bounds, {"_FillValue": False}),
    ),
    coords=dict(  # pylint: disable=R1735
        time=(
            "time",
            time_array,
            {
                "units": time_units,
                "axis": "T",
                "calendar": calendar,
                "bounds": "time_bounds",
                "_FillValue": False,
            },
        )
    ),
)
ds.to_netcdf(infile)
