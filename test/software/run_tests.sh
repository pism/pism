#!/bin/bash

TOOL=$1
SILENT=1
# Write the summary to this file:
SUMMARY=`pwd`/"$TOOL.txt"
# Remove the old summary file, it is exists.
rm -f $SUMMARY

# Load some functions
source functions.sh

# Print the summary title.
print "Testing $TOOL."

PASSCOUNT=0
FAILCOUNT=0
for each in *;
do
    if [ -d $each ] && [ $each == $TOOL ];
    then
	cd $each
	for i in *.sh;
	do
	    source $i
	done
    fi
done

print "Summary:"
print "PASSED: $PASSCOUNT out of $(( $PASSCOUNT + $FAILCOUNT ))."
print "FAILED: $FAILCOUNT out of $(( $PASSCOUNT + $FAILCOUNT ))."