#!/bin/sh

# HEINOCATDAT.SH is a shell script to concatenate .dat files generated from
# ISMIP-HEINO runs on PISM.  (I.e. "pisms -ismip H -datprefix FOO ...").
# Execute as 
#     test/catdat ST FOO BAR JOE
# if files are FOO_ST_ts*.dat and BAR_ST_ts*.dat.  Generates files named
# JOE_ST_ts*.dat.

#ELB 11/1/06

echo "attempting to create files $4_$1_ts*.dat:"
echo -n "   writing: "

count=0

for name in ts_iv ts_tba tss_ait tss_ahbt tss_tba
do
  if [ -e $2_$1_$name.dat ]
  then
    if [ -e $3_$1_$name.dat ]
    then
      cat $2_$1_$name.dat $3_$1_$name.dat >> $4_$1_$name.dat
      echo -n $name
      echo -n " "
      count=$(($count + 1))
    else
      echo "   [expected file $3_$1_$name.dat not found]"
    fi
  else
    echo "   [expected file $2_$1_$name.dat not found]"
  fi
done

for nameend in it hbt bfh
do
  for ((myn=1 ; myn <= 7 ; myn++))
  do
    if [ -e $2_$1_tsp$myn\_$nameend.dat ]
    then
      if [ -e $3_$1_tsp$myn\_$nameend.dat ]
      then
        cat $2_$1_tsp$myn\_$nameend.dat $3_$1_tsp$myn\_$nameend.dat\
              >> $4_$1_tsp$myn\_$nameend.dat
        echo -n tsp$myn\_$nameend
        echo -n " "
        count=$(($count + 1))
      else
        echo "   [expected file $3_$1_tsp$myn\_$nameend.dat not found]"
      fi
    else
      echo "   [expected file $2_$1_tsp$myn\_$nameend.dat not found]"
    fi
  done
done

echo
echo "created $count files"
