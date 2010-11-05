#!/bin/bash

SUBDIR=pism
SILENT=1
# Write the summary to this file:
SUMMARY=$PWD/"$SUBDIR.txt"
# Remove the old summary file, it is exists.
rm -f $SUMMARY

# Load some functions
source functions.sh

EXEC=`which pismr`
PREFIX=`dirname $EXEC`

# Modify the PATH variable so that we always use the right version of PISM:
PATH=$PREFIX/bin:$PREFIX/util:$PATH

# Get PISM version:
# Note that 'cat' is here so that the whole output of pismr is read and no
# "pipe is broken" error is generated:
VERSION=`$MPIDO -n 1 pismr -pismversion | cat | head -1 | sed -e "s/PISMR //; s/ (.\{1,\}//"`

# Print the summary title.
print "Testing $SUBDIR $VERSION in $PREFIX."

PASSCOUNT=0
FAILCOUNT=0
cd $SUBDIR
for i in "$@";
do
    source test_${i}.sh
done

print "Summary:"
print "PASSED: $PASSCOUNT out of $(( $PASSCOUNT + $FAILCOUNT ))."
print "FAILED: $FAILCOUNT out of $(( $PASSCOUNT + $FAILCOUNT ))."
