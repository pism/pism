#!/bin/bash

# This test script ensures that parameters in pism_config.cdl remain sorted.

pism_config_cdl=${1:?Usage: $0 /path/to/pism_config.cdl}

temp_dir=`mktemp -d -t pism_test.XXXX` || exit 1

cat ${pism_config_cdl} | \
  grep "pism_config:" | \
  grep -v "long_name" | \
  grep -v "_doc =" | \
  grep -v "_units =" | \
  grep -v "_type =" | \
  grep -v "_option =" | \
  grep -v "_choices =" | \
  sed "s/^ *pism_config://" | \
  sed "s/ =.*//" > ${temp_dir}/raw_config.txt

# Avoid locale-dependent "sort" results:
export LC_COLLATE=C

cat ${temp_dir}/raw_config.txt | sort > ${temp_dir}/sorted_config.txt

diff -u ${temp_dir}/raw_config.txt ${temp_dir}/sorted_config.txt
