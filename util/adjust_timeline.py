#!/usr/bin/env python3

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

from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

import cftime
import numpy as np
import xarray as xr
from dateutil.parser import parse

# Set up the option parser
parser = ArgumentParser()
parser.description = """Script adjusts the time file with time and time
bounds that can be used to determine to force PISM via command line
option -time_file or adjust the time axis for postprocessing."""
parser.add_argument("FILE", nargs="*")
parser.add_argument(
    "-p",
    "--periodicity",
    dest="periodicity",
    help="""periodicity, e.g. monthly, daily, etc. Default=MS""",
    default="MS",
)
parser.add_argument(
    "-a",
    "--start_date",
    dest="start_date",
    help="""Start date in ISO format. Default=1980-1-1""",
    default="1980-1-1",
)
parser.add_argument(
    "-c",
    "--calendar",
    dest="calendar",
    choices=["standard", "gregorian", "no_leap", "365_day", "360_day", "julian"],
    help="""Sets the calendar. Default="standard".""",
    default="standard",
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
    help="""Reference unit. Default=days. Use of months or
                    years is NOT recommended.""",
    default="seconds",
)
parser.add_argument(
    "-d",
    "--ref_date",
    dest="ref_date",
    help="""Reference date. Default=1960-1-1""",
    default="1960-1-1",
)

options = parser.parse_args()
interval_type = options.interval_type
periodicity = options.periodicity.upper()
start_date = parse(options.start_date)
ref_unit = options.ref_unit
ref_date = options.ref_date
args = options.FILE
infile = Path(args[0])

time_var_name = "time"
time_bnds_var_name = "time_bounds"
bnds_var_name = "bnds"

time_units = "%s since %s" % (ref_unit, ref_date)
calendar = options.calendar

with xr.open_dataset(infile) as ds:
    nt = len(ds[time_var_name])

date_range = xr.cftime_range(
    start_date, periods=nt + 1, freq=periodicity, calendar=calendar
)
dates = np.array(cftime.date2num(date_range, time_units, calendar=calendar))
# reference date from command-line argument
r = time_units.split(" ")[2].split("-")
refdate = datetime(int(r[0]), int(r[1]), int(r[2]))

if interval_type == "mid":
    # mid-point value:
    # time[n] = (bnds[n] + bnds[n+1]) / 2
    time = dates[0:-1] + (dates[1::] - dates[0:-1]) / 2
elif interval_type == "start":
    time = dates[0:-1]
else:
    time = dates[1::]

time_bounds = np.zeros((nt, 2))
time_bounds[:, 0] = dates[0:-1]
time_bounds[:, 1] = dates[1::]

time = np.floor(time)
time_bounds = np.floor(time_bounds)

ds[time_var_name] = (
    [time_var_name],
    time,
    {
        "bounds": time_bnds_var_name,
        "_FillValue": False,
        "units": time_units,
        "calendar": calendar,
        "axis": "T",
        "long_name": "time",
        "standard_name": "time",
    },
)
if hasattr(ds[time_var_name], "bounds"):
    bounds_var = ds[time_var_name].bounds
    ds[bounds_var] = (
        [time_var_name, bnds_var_name],
        time_bounds,
        {"_FillValue": False},
    )
else:
    ds[time_bnds_var_name] = (
        [time_var_name, bnds_var_name],
        time_bounds,
        {"_FillValue": False},
    )
ds.to_netcdf(infile, mode="a")
