#!/usr/bin/env python3

import xarray as xr
import os
import shutil
import sys

PISM_PATH = sys.argv[1]
PISMR = os.path.join(PISM_PATH, "pismr")

def run(command):
    print(command)
    return os.system(command)


files = ["eisII-output.nc",
         "pismr-output.nc",
         "both-consistent.nc",
         "both-string-missing.nc",
         "both-string-mismatch.nc",
         "both-double-missing.nc",
         "both-double-mismatch.nc"]

# create an input file
run(PISMR + " -eisII A -verbose 1 -Mx 3 -My 3 -Mz 5 -y 10 -o eisII-output.nc")

# add the PROJ string
nc = xr.open_dataset("eisII-output.nc")
nc.attrs["proj"] = "epsg:3413"
nc.to_netcdf("eisII-output.nc")

print("Test running PISM initialized from a file w/o mapping  but with proj...")
assert run(PISMR + " -verbose 1 -i eisII-output.nc -y 10 -o both-consistent.nc") == 0

print("Test that the mapping variable was initialized using the proj attribute...")
nc = xr.open_dataset("both-consistent.nc")
mapping = nc.variables["mapping"]
assert mapping.attrs["grid_mapping_name"] == "polar_stereographic"

print("Test re-starting PISM with consistent proj and mapping...")
assert run(PISMR + " -verbose 1 -i both-consistent.nc -o pismr-output.nc") == 0

# remove a required string attribute
shutil.copy("both-consistent.nc", "both-string-missing.nc")
nc = xr.open_dataset("both-string-missing.nc")
del nc.variables["mapping"].attrs["grid_mapping_name"]
nc.to_netcdf("both-string-missing.nc")

print("Test that PISM stops if a required string attribute is missing...")
assert run(PISMR + " -verbose 1 -i both-string-missing.nc -o pismr-output.nc") != 0

# alter a required string sttribute
shutil.copy("both-consistent.nc", "both-string-mismatch.nc")
nc = xr.open_dataset("both-string-mismatch.nc")
nc.variables["mapping"].attrs["grid_mapping_name"] = "wrong"
nc.to_netcdf("both-string-mismatch.nc")

print("Test that PISM stops if a required string attribute has a wrong value...")
assert run(PISMR + " -verbose 1 -i both-string-mismatch.nc -o pismr-output.nc") != 0

# remove a required double attribute
shutil.copy("both-consistent.nc", "both-double-missing.nc")
nc = xr.open_dataset("both-double-missing.nc")
del nc.variables["mapping"].attrs["standard_parallel"]
nc.to_netcdf("both-double-missing.nc")

print("Test that PISM stops when a required double attribute is missing...")
assert run(PISMR + " -verbose 1 -i both-double-missing.nc -o pismr-output.nc") != 0

# alter a required double attribute
shutil.copy("both-consistent.nc", "both-double-mismatch.nc")
nc = xr.open_dataset("both-double-mismatch.nc")
nc.variables["mapping"].attrs["standard_parallel"] = 45.0
nc.to_netcdf("both-double-mismatch.nc")

print("Test that PISM stops if a required double attribute has a wrong value...")
assert run(PISMR + " -verbose 1 -i both-double-mismatch.nc -o pismr-output.nc") != 0

# cleanup
for f in files:
    print("Removing %s..." % f)
    os.remove(f)
