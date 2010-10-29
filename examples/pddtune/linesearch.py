#!/usr/bin/env python

## @package linesearch
## \author Ed Bueler, University of Alaska Fairbanks, USA
## \brief A script for doing a line search over stddev parameter of PDD.
## Copyright (C) 2010 Ed Bueler
##
## see README for role of this script
## This script uses NCO (http://nco.sourceforge.net/).

import commands
from numpy import array, double, int
from getopt import getopt, GetoptError
from sys import argv, exit

# FIXME : need to actually implement a line search!  for now this is just
#   "runcase.sh" re-written into python, with weights.py incorporated

usage="""
LINESEARCH.PY  Do bisection search on stddev parameter of PDD.  Abort search
if stddev=0 and stddev=10 do not cause mean smb error to bracket zero.

Imports computeobjective() method from stand-alone objective.py.

Examples:
  ./linesearch.py --thresh=268 --snow=0.001 --refreeze=0.4 --stddev=1.0 \
      --diffsfile=foo.txt --startfile=start.nc
  
  ./linesearch.py --help   # print this message
"""

def usagefailure(message):
    print message
    print
    print usage
    exit(2)

# FIXME:  make this controllable at command line
# name of PISM file with Greenland geometry and precip,smb from Ettema et al.
#   and other needed info to run pclimate:
DATANAME = "Greenland_5km_v1.1.nc"
PISMDATA = "pism_" + DATANAME


if __name__ == "__main__":
    try:
      opts, args = getopt(argv[1:], "", 
                          ["thresh=", "snow=", "refreeze=","stddev=",
                           "diffsfile=", "startfile=", "deletenc", "help","usage"])
    except GetoptError:
      usagefailure('ERROR: INCORRECT COMMAND LINE ARGUMENTS FOR linesearch.py')
    threshold = 273.15
    ddfsnow = 0.003
    refreeze = 0.6
    stddev = 2.53
    diffsfile = ""
    startfile = "start.nc"
    deletencfiles = False
    for (opt, optarg) in opts:
        if opt in ("--thresh"):
            threshhold = float(optarg)
        if opt in ("--snow"):
            ddfsnow = float(optarg)
        if opt in ("--refreeze"):
            refreeze = float(optarg)
        if opt in ("--stddev"):
            stddev = float(optarg)
        if opt in ("--diffsfile"):
            diffsfile = optarg
        if opt in ("--startfile"):
            startfile = optarg
        if opt in ("--deletenc"):
            deletencfiles = True
        if opt in ("--help", "--usage"):
            print usage
            exit(0)

    # lock together ddfs for snow and ice; FIXME: allow ratio choice
    ddfice = 2.0 * ddfsnow
    
    nameroot = "%.2f_%.4f_%.2f_%.2f" % (threshhold, ddfsnow, refreeze, stddev)
    print "LINESEARCH.PY CASE  (threshhold=%.2f, ddfsnow=%.4f, refreeze=%.2f)" \
              % (threshhold, ddfsnow, refreeze)
    print ""

    configfile = "config_" + nameroot + ".nc"
    configopt = " -config_override  " + configfile
    print "  creating -config_override file %s ..." % configfile

    try:
      from netCDF4 import Dataset as NC
    except:
      from netCDF3 import Dataset as NC

    #FIXME: need this before creation?:  rm -rf configfile

    try:
      nc_config = NC(configfile, 'w')
    except:
      usagefailure("ERROR: NETCDF FILE '%s' CANNOT BE OPENED FOR WRITING" % configfile)

    overs = nc_config.createVariable("pism_overrides", 'b')  # variable type is NC_BYTE

    overs.pdd_positive_threshold_temp = threshhold
    overs.pdd_factor_snow = ddfsnow
    overs.pdd_factor_ice = ddfice
    overs.pdd_refreeze = refreeze
    overs.pdd_std_dev = stddev
    
    overs[:] = 0
    
    nc_config.close()
    
    # this is output of pclimate, which will be evaluated against PISMDATA
    climatefile = "clim_" + nameroot + ".nc"  

    # change this to "mpiexec -n 8" or similar to run on multiple processes
    mpido=""

    # coupler settings: Fausto 2m air temp parameterization, but default PDD
    #   (w/o Greve/Fausto settings of PDD parameters)
    coupleropt = " -atmosphere searise_greenland -surface pdd"

    climstartopt = " -ys 1990"
    climendopt = " -ye 1991"

    #dt = 0.0833333333 # monthly = (1/12) of year
    dtopt = " -dt 1.0"

    if len(diffsfile) > 0:
      print "  will add lines of text to " + diffsfile + " ..."
    if deletencfiles:
      print "  will delete NetCDF files %s and %s when no longer needed ..." % \
              (configfile, climatefile)

    # run PISM:
    command = mpido + " pclimate -i " + startfile + coupleropt + configopt \
              + climstartopt + climendopt + dtopt + " -o " + climatefile
    print "  doing:"
    print "    " + command
    try:
        (status,output) = commands.getstatusoutput(command)
    except KeyboardInterrupt:
        exit(2)
    if status:
        #exit(status)
        print status
        print output

    # run objective.py
    from objective import computeobjective

    countvalid, avdiff, avL2, av2kL2 = computeobjective(
                                         startfile, climatefile, PISMDATA,
                                         "acab", "smb")
    print "  result: "    
    print "    %d locations for valid (thk>0) comparison:" % countvalid
    print "    average of signed differences (whole sheet) is  %12.7f" % avdiff
    print "    average of squared differences (whole sheet) is %12.7f" % avL2
    print "    average of squared differences (H < 2000) is    %12.7f" % av2kL2

    # FIXME:  allow choice of weights
    weighted = 10.0 * avdiff**2.0 + 1.0 * avL2 + 3.0 * av2kL2
    print "    weighted average of above quantities is         %12.7f" % weighted
    
    # write results to file if a file name was given, otherwise stdout
    lineoftxt = "%s %12.7f %12.7f %12.7f %12.7f\n" % \
        (climatefile, avdiff, avL2, av2kL2, weighted)
    if len(diffsfile) > 0:
      diffs = file(diffsfile,'a')
      diffs.write(lineoftxt)
      diffs.close()
    else:
      print "  result in one line:"
      print lineoftxt

    # finalize: if option was given, clean up generated NetCDF files
    if deletencfiles:
      command = "rm -rf " + configfile + " " + climatefile
      print "  doing:"
      print "    " + command
      try:
          (status,output) = commands.getstatusoutput(command)
      except KeyboardInterrupt:
          exit(2)
      if status:
          #exit(status)
          print status
          print output

    print ""

