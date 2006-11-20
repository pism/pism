#!/bin/sh
# VERIFYNOW is a script to do fairly quick verification of the isothermal and 
# thermocoupled shallow ice components of ice model COMMVNISM.  Uses only tests
# C and G.  Intended to do minimal amount of computation to show convergence 
# to continuum results.  Takes about [ABOUT 6 HOURS??! on my laptop]
# ELB 7/3/06; 11/19/06

echo "++++++++ starting with isothermal test C and Mx=31,61,121:"

for myMx in 31 61 121
do
   obj/pismv -test C -Mx $myMx -My $myMx -Mz 31 -y 15208.0 -verbose > _temp_result.txt 
   sed '/  history =/,+1!d' _temp_result.txt | sed 1d
   sed '/Actual ERRORS/,+2!d' _temp_result.txt
   date
done

echo "+++++++++ continuing with thermocoupled test G and Mx=61,91,121:"

for myMx in 61 91 121
do
   obj/pismv -test G -Mx $myMx -My $myMx -Mz $myMx -y 25000.0 -verbose > _temp_result.txt 
   sed '/  history =/,+1!d' _temp_result.txt | sed 1d
   sed '/Actual ERRORS/,+4!d' _temp_result.txt
   date
done

# comment next line to keep temporary file
rm -f _temp_result.txt

# illustrates redirect to both standard out and (append) to a file
#./verify -test F -Mx 61 -My 61 -Mz 41 -y 100.0 |tee result.txt
