#!/bin/bash
 
# Simply creates useful .png figures from MISMIP steady state output files,
# i.e. ABC1_1a_M1_A1_ss and ABC1_1a_M1_A1_f.  Applies python script 
# figsMISMIP.py to all such files.

for file in *_ss
do
  prefix=`echo $file | sed 's/_ss//g'`  # strip "_ss" from file name
  ./figsMISMIP.py -p $prefix -e
done

