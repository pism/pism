#!/bin/bash
 
# Script for MISMIP = Marine Ice Sheet Intercomparison Project.  See
#    http://homepages.ulb.ac.be/~fpattyn/mismip/
# and especially
#    http://homepages.ulb.ac.be/~fpattyn/mismip/mismip_4.pdf
# (or latest version of intercomparison description).

# This intercomparison describes a flow line ice stream and ice shelf
# system.  PISM is a 3D model so there is a notable lack of efficiency
# relative to flow line models.  Furthermore there is (for now) no 
# boundary layer or other advanced treatment of dynamics in the vicinity 
# of the grounding line; compare C. Schoof, (2007). "Marine ice-sheet
# dynamics. Part 1. The case of rapid sliding," J. Fluid Mech., vol. 573,
# pp. 27–55.

# See comments at start of src/ismip/iceMISMIPModel.cc regarding grounding
# line flux issues.

# A PISM MISMIP run always saves an ASCII file in the MISMIP specified format:
#   ABC1_1b_M1_A1_t
# and a standard PISM NetCDF file
#   ABC1_1b_M1_A1.nc
# If the run achieves steady state by the MISMIP criterion, then the run also
# saves two more ASCII files in the MISMIP spec format:
#   ABC1_1b_M1_A1_ss
#   ABC1_1b_M1_A1_f
# See intercomparison description (URL above) for format of these files.


NN=2           # set default number of processors here

MPIDO=mpiexec  # change this if "mpirun", etc.

MYINITIALS=EBU # MISMIP says "first character of the first name followed
               #    by the first two characters of the last name"


if [ $# -gt 0 ] ; then  # if user says "mismip.sh 8" then NN = 8
  NN="$1"
fi

SHOWONLY=0
if [ $# -gt 1 ] ; then  # if user says "mismip.sh 8 D" then NN = 8 and only shows; no run 
  SHOWONLY=1
fi

set -e  # exit on error

# function to run "pisms -mismip EXPER -initials MYINITIALS OTHEROPTIONS"
#   on NN processors
mpimismip()
{
    # change this if "bin/pisms", etc.:
    cmd="$MPIDO -np $1 pisms -mismip $2 -initials $MYINITIALS $3"
    
    echo 
    echo "date = '`date`' on host '`uname -n`':"
    echo "trying '$cmd'"
    echo
    if [ $SHOWONLY = 0 ] ; then
      $cmd
    fi
}


#uncomment these (and "#fi" at end) and move around to bypass completed stuff:
#if [ ]; then   # always goes to "else"
#else           # put this before restart location


# EBU1_1a_M1_A1_*
mpimismip $NN 1a "-run 1 -Mx 3 -My 151 -Mz 11 -ksp_rtol 1e-11"

# EBU2_1a_M1_A1_*
mpimismip $NN 1a "-super -run 1 -Mx 3 -My 151 -Mz 11 -ksp_rtol 1e-11"

# EBU1_1a_M2_A1_*
mpimismip $NN 1a "-run 1 -Mx 3 -My 1501 -Mz 11 -ksp_rtol 1e-11"

exit

# EBU1_3a_M1_A1_*
mpimismip $NN 3a "-run 1 -Mx 3 -My 151 -Mz 11 -ksp_rtol 1e-11"

# EBU2_3a_M1_A1_*
mpimismip $NN 3a "-super -run 1 -Mx 3 -My 151 -Mz 11 -ksp_rtol 1e-11"

#fi

