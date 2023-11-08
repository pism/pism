#!/bin/bash

# Run *very* short versions of spinup simulations to make sure that command-line options
# are up to date, etc.

set -e

NOMASS_RUN_LENGTH=10 CONSTANT_CLIMATE_RUN_LENGTH=10 ./antspin-coarse.sh 8

RUN_LENGTH=10 ./antspin-regridtofine.sh 8
