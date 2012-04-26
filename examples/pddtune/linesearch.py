#!/usr/bin/env python

## @package linesearch
## \author Ed Bueler, University of Alaska Fairbanks, USA
## \brief A script for doing a line search over stddev parameter of PDD.
## Copyright (C) 2010, 2012 Ed Bueler
##
## see README for role of this script
## This script uses NCO (http://nco.sourceforge.net/).

import commands
from numpy import array, double, int
from getopt import getopt, GetoptError
from sys import argv, exit

usage="""
LINESEARCH.PY  Do line search (bisection and false position) on stddev parameter
of the PDD model for surface mass balance.  Abort search if stddev=0.5 and
stddev=10.0 do not make pclimate output total smb bracket the target smb.

Imports computeobjective() method from stand-alone objective.py.

Examples:

1) do a linesearch:
  ./linesearch.py --thresh=268 --snow=0.001 --refreeze=0.4 --lapse=0.3 \
      --diffsfile=diffs.txt --startfile=start.nc

2) reproduce a file from its name:
  ./linesearch.py --reproduce=clim_268_0.001_0.4_0.300_1.0.nc --startfile=start.nc
  
3) get help; print this message:
  ./linesearch.py --help

See:  dotune.sh, boot.sh, objective.py, README.
"""

def usagefailure(message):
    print message
    print
    print usage
    exit(2)

class Case:
    """Holds parameters for one pclimate case.  Has functionality related to
naming scheme for pclimate configuration and output NetCDF files, e.g.
'clim_268.0_0.002_0.80_0.001_1.688.nc' and similar.  For this name, the
'nameroot' is '268.0_0.002_0.80_0.001_1.688'."""
    def __init__(self):
      self.threshhold = 273.15
      self.ddfsnow = 0.003
      self.ddfice = 2.0 * self.ddfsnow
      self.refreeze = 0.6
      self.stddev_lapse = 0.0
      self.stddev = 2.53
      self.annualizepdd = False
    def update_ddfice(self):
      self.ddfice = 2.0 * self.ddfsnow
    def get_nameroot(self):
      self.update_ddfice()
      return "%.1f_%.3f_%.2f_%.3f_%.4f" \
         % (self.threshhold, self.ddfsnow, self.refreeze, self.stddev_lapse, self.stddev)
    def put_nameroot(self,namerootIN):
      p = namerootIN.split('_')
      if len(p) < 5:
        print "ERROR: nameroot provided to Case.put_nameroot has too few params"
        exit(1)
      self.threshhold = float(p[0])
      self.ddfsnow = float(p[1])
      self.refreeze = float(p[2])
      self.stddev_lapse = float(p[3])
      self.stddev = float(p[4])
      self.update_ddfice()      
    def put_pism_overrides(self,nc):
      """input nc is an open NetCDF file"""
      self.update_ddfice()
      # variable type is NC_BYTE
      overs = nc.createVariable("pism_overrides", 'b')  
      overs.pdd_positive_threshold_temp = self.threshhold
      overs.pdd_factor_snow = self.ddfsnow
      overs.pdd_factor_ice = self.ddfice
      overs.pdd_refreeze = self.refreeze
      overs.pdd_std_dev_lapse_lat_rate = self.stddev_lapse
      overs.pdd_std_dev = self.stddev
      overs[:] = 0

class Files:
    def __init__(self):
      # it might be good to have this controllable at command line:
      # name of PISM file with Greenland geometry and precip,smb from Ettema et al.
      #   and other needed info to run pclimate:
      DATANAME = "Greenland_5km_v1.1.nc"
      self.targetfile = "pism_" + DATANAME
      self.startfile = "start.nc"  # effectively allow user to forget --startfile=
      self.diffsfile = ""  # if user forgets --diffsfile=foo.txt then writes to stdout
    def copy(self):
      cpfn = Files()
      cpfn.targetfile = self.targetfile
      cpfn.startfile = self.startfile
      cpfn.diffsfile = self.diffsfile
      return cpfn
    def configfile(self,cp):
      "this is output of pclimate, which will be evaluated against PISMDATA"
      return "config_" + cp.get_nameroot() + ".nc"
    def climatefile(self,cp):
      "this is the configuration file for the -config_override mechanism"
      return "clim_" + cp.get_nameroot() + ".nc"
    def printme(self, cp):
      print "**PRINTING Files OBJECT:**"
      print "  targetfile      = " + self.targetfile
      print "  startfile       = " + self.startfile
      print "  diffsfile       = " + self.diffsfile
      print "  configfile(cp)  = " + self.configfile(cp)
      print "  climatefile(cp) = " + self.climatefile(cp)
      print "**END**"

def evalcase(stddev, cp, fn, deletencfiles):
    """evaluates one pclimate run against smb target in a file"""
    # input cp is of type Case()
    # input fn is of type Filenames()

    cp.stddev = stddev

    try:
      from netCDF4 import Dataset as NC
    except:
      from netCDF3 import Dataset as NC

    from objective import computeobjective

    configopt = " -config_override " + fn.configfile(cp)
    print "  creating -config_override file %s ..." % fn.configfile(cp)
    try:
      nc_config = NC(fn.configfile(cp), 'w')
    except:
      usagefailure("ERROR: NETCDF FILE '%s' CANNOT BE OPENED FOR WRITING" \
                   % fn.configfile(cp) )
    cp.put_pism_overrides(nc_config)
    nc_config.close()

    # change this to "mpiexec -n 8" or similar to run on multiple processes
    mpido=""
    # coupler settings: Fausto 2m air temp parameterization, but default PDD
    #   (w/o Greve/Fausto settings of PDD parameters)
    coupleropts = " -atmosphere searise_greenland -surface pdd"
    if cp.annualizepdd:
      coupleropts += " -pdd_annualize"
    timeopts = " -times 1990:1.0:1991"
    #dt = 0.0833333333 # monthly = (1/12) of year

    if len(fn.diffsfile) > 0:
      print "  will add lines of text to " + fn.diffsfile + " ..."
    if deletencfiles:
      print "  will delete NetCDF files %s and %s when no longer needed ..." % \
              (fn.configfile(cp), fn.climatefile(cp))

    # run PISM:
    command = mpido + " pclimate -i " + fn.startfile + coupleropts + configopt \
                + timeopts + " -o " + fn.climatefile(cp)
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


    countvalid, avdiff, avL2, av2kL2 = computeobjective(
                                           fn.startfile, fn.climatefile(cp), 
                                           fn.targetfile,
                                           "acab", "smb")
    print "  result: "    
    print "    %d locations for valid (thk>0) comparison:" % countvalid
    print "    average of signed differences (whole sheet) is  %12.7f" % avdiff
    print "    average of squared differences (whole sheet) is %12.7f" % avL2
    print "    average of squared differences (H < 2000) is    %12.7f" % av2kL2

    # FIXME:  allow choice of weights
    weighted = 1.0 * avL2 + 3.0 * av2kL2
    print "    weighted average of above quantities is         %12.7f" % weighted

    # write results to file if a file name was given, otherwise stdout
    lineoftxt = "%s %12.7f %12.7f %12.7f %12.7f\n" % \
                  (fn.climatefile(cp), avdiff, avL2, av2kL2, weighted)
    if len(fn.diffsfile) > 0:
      diffs = file(fn.diffsfile,'a')
      diffs.write(lineoftxt)
      diffs.close()
    else:
      print "  result in one line:"
      print lineoftxt

    # finalize: if option was given, clean up generated NetCDF file
    if deletencfiles:
      command = "rm -rf " + fn.configfile(cp) + " " + fn.climatefile(cp)
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
    return avdiff

def reproduce(climfilename, fn):
    """parses a pclimate output file name, and runs evalcase() to reproduce that run"""
    # input fn is of type Filenames()
    
    if not(climfilename.endswith('.nc')):
      print "ERROR: filename needs to be a NetCDF file"
      exit(1)
    if not(climfilename.startswith('clim_')):
      print "WARNING: filename was expected to start with 'clim_'"
      print "         (in any case, filename will be stripped through first '_' to"
      print "         generate the nameroot)"
    nr = climfilename[climfilename.find('_')+1:-3]
    cp = Case()
    cp.put_nameroot(nr)
    # fn.printme(cp)
    result = evalcase(cp.stddev, cp, fn, False)
    print "************ result at stddev=%.5f is %.5f ***********" \
          % (cp.stddev, result)
    

if __name__ == "__main__":

    cp = Case()
    fn = Files()

    try:
      opts, args = getopt(argv[1:], "", 
                          ["thresh=", "snow=", "refreeze=", "lapse=",
                           "sd=",  # <-- if this is an option then don't do linesearch
                           "tol=",
                           "reproduce=",
                           "diffsfile=", "startfile=",
                           "annualize", "deletenc", "help","usage"])
    except GetoptError:
      usagefailure('ERROR: INCORRECT COMMAND LINE ARGUMENTS FOR linesearch.py')
    dolinesearch = True
    cp.annualizepdd = False
    deletencfiles = False
    avdifftol = 0.001
    merely_reproduce = False
    reproname = ""
    for (opt, optarg) in opts:
        if opt in ("--thresh"):
            cp.threshhold = float(optarg)
        if opt in ("--snow"):
            cp.ddfsnow = float(optarg)
        if opt in ("--refreeze"):
            cp.refreeze = float(optarg)
        if opt in ("--lapse"):
            cp.stddev_lapse = float(optarg)
        if opt in ("--sd"):
            dolinesearch = False
            cp.stddev = float(optarg)
        if opt in ("--annualize"):
            cp.annualizepdd = True
        if opt in ("--tol"):
            avdifftol = float(optarg)
        if opt in ("--diffsfile"):
            fn.diffsfile = optarg
        if opt in ("--startfile"):
            fn.startfile = optarg
        if opt in ("--deletenc"):
            deletencfiles = True
        if opt in ("--reproduce"):
            merely_reproduce = True
            reproname = optarg
        if opt in ("--help", "--usage"):
            print usage
            exit(0)

    if merely_reproduce:
      if deletencfiles:
        print "WARNING: option --deletenc is ignored; output files will remain"
      reproduce(reproname, fn.copy())
      exit(0)
    
    print "LINESEARCH.PY CASE  (threshhold=%.2f, ddfsnow=%.4f, refreeze=%.2f, sdlapse=%.3f)" \
              % (cp.threshhold, cp.ddfsnow, cp.refreeze, cp.stddev_lapse)
    print ""

    if not(dolinesearch):
      result = evalcase(cp.stddev, cp, fn, False)
      print "************ result at stddev=%.5f is %.5f ***********" \
                % (cp.stddev, result)
      exit(0)
      
    # line search:  combines bisection with false position at end
    Asd = 0.5
    Bsd = 10.0
    stddevtol = 1.0e-5 * (Bsd - Asd);
    
    F_Asd = evalcase(Asd, cp, fn, deletencfiles)
    F_Bsd = evalcase(Bsd, cp, fn, deletencfiles)

    if F_Asd * F_Bsd > 0:
      print "************ at stddev in [%.5f,%.5f], NO BRACKET ***********" \
               % (Asd,Bsd)
    else:
      count = 2
      maxcount = 20  # max number of function evals per case
      while True:
        if count < 5:
          # bisection for a few steps
          Csd = 0.5 * (Asd + Bsd)
        else:
          # now false position
          if abs(Bsd - Asd) <= stddevtol:
            break
          Csd = Asd - (F_Asd) * (Bsd - Asd) / (F_Bsd - F_Asd);
        result = evalcase(Csd, cp, fn, deletencfiles)
        print "************ result at stddev=%.5f is %.5f ***********" \
                % (Csd, result)
        if abs(result) <= avdifftol:
          break
        count = count + 1
        if count >= maxcount:
          print "************ max number of function evals reached in bisection **********"
          break
        if F_Asd * result > 0.0:
          Asd = Csd
          F_Asd = result
        else:
          Bsd = Csd
          F_Bsd = result

