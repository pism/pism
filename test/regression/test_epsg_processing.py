#!/usr/bin/env python

from netCDF4 import Dataset
import os
import shutil
import sys

PISM_PATH = sys.argv[1]
PISMR = os.path.join(PISM_PATH, "pismr")
PISMS = os.path.join(PISM_PATH, "pisms")


def run(command):
    print(command)
    return os.system(command)


files = ["pisms-output.nc",
         "pismr-output.nc",
         "both-consistent.nc",
         "both-string-missing.nc",
         "both-string-mismatch.nc",
         "both-double-missing.nc",
         "both-double-mismatch.nc"]

# create an input file
run(PISMS + " -verbose 1 -Mx 3 -My 3 -Mz 5 -y 10 -o pisms-output.nc")

# add the PROJ.4 string
nc = Dataset("pisms-output.nc", "a")
nc.proj4 = "+init=epsg:3413"
nc.close()

print("Test running PISM initialized from a file w/o mapping  but with proj4...")
assert run(PISMR + " -verbose 1 -i pisms-output.nc -y 10 -o both-consistent.nc") == 0

print("Test that the mapping variable was initialized using the proj4 attribute...")
nc = Dataset("both-consistent.nc", "r")
mapping = nc.variables["mapping"]
assert mapping.grid_mapping_name == "polar_stereographic"
nc.close()

print("Test re-starting PISM with consistent proj4 and mapping...")
assert run(PISMR + " -verbose 1 -i both-consistent.nc -o pismr-output.nc") == 0

# remove a required string attribute
shutil.copy("both-consistent.nc", "both-string-missing.nc")
nc = Dataset("both-string-missing.nc", "a")
del nc.variables["mapping"].grid_mapping_name
nc.close()

print("Test that PISM stops if a required string attribute is missing...")
assert run(PISMR + " -verbose 1 -i both-string-missing.nc -o pismr-output.nc") != 0

# alter a required string sttribute
shutil.copy("both-consistent.nc", "both-string-mismatch.nc")
nc = Dataset("both-string-mismatch.nc", "a")
nc.variables["mapping"].grid_mapping_name = "wrong"
nc.close()

print("Test that PISM stops if a required string attribute has a wrong value...")
assert run(PISMR + " -verbose 1 -i both-string-mismatch.nc -o pismr-output.nc") != 0

# remove a required double attribute
shutil.copy("both-consistent.nc", "both-double-missing.nc")
nc = Dataset("both-double-missing.nc", "a")
del nc.variables["mapping"].standard_parallel
nc.close()

print("Test that PISM stops when a required double attribute is missing...")
assert run(PISMR + " -verbose 1 -i both-double-missing.nc -o pismr-output.nc") != 0

# alter a required double attribute
shutil.copy("both-consistent.nc", "both-double-mismatch.nc")
nc = Dataset("both-double-mismatch.nc", "a")
nc.variables["mapping"].standard_parallel = 45.0
nc.close()

print("Test that PISM stops if a required double attribute has a wrong value...")
assert run(PISMR + " -verbose 1 -i both-double-mismatch.nc -o pismr-output.nc") != 0

# cleanup
for f in files:
    print("Removing %s..." % f)
    os.remove(f)
