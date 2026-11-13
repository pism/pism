#!/bin/bash

set -e
set -u
set -x

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

pism="${MPIEXEC} -n 4 ${PISM_PATH}/pism"

regular=${PISM_SOURCE_DIR}/test/regression/pico_split/bedmap2_schmidtko14_50km.nc
inverted=ant_inverted.nc

# create a temporary directory and set up automatic cleanup
temp_dir=$(mktemp -d --tmpdir pism-regrid-invertlat-XXXX)
echo $temp_dir
trap 'rm -rf "$temp_dir"' EXIT
cd $temp_dir

cat > invert_y.py <<EOF
#!/usr/bin/env python3
from netCDF4 import Dataset
import sys

def copy_attrs(v1, v2):
    for attr in v1.ncattrs():
        v2.setncattr(attr, v1.getncattr(attr))

def copy_variable(ds1, ds2, name, flip_dim):
    v1 = ds1.variables[name]

    for dim in v1.dimensions:
        dv1 = ds1.variables[dim]

        ds2.createDimension(dim, ds1.dimensions[dim].size)
        dv2 = ds2.createVariable(dim, 'f8', (dim,),)

        copy_attrs(dv1, dv2)

        dv2[:] = dv1[::-1] if dim == flip_dim else dv1[:]

    v2 = ds2.createVariable(name, 'f8', v1.dimensions)
    copy_attrs(v1, v2)

    data = v1[:]
    for k, d in enumerate(v1.dimensions):
        if d == flip_dim:
            slab = len(v1.dimensions) * [slice(None)]
            slab[k] = slice(None, None, -1)
            data = data[tuple(slab)]
            break

    v2[:] = data

with Dataset(sys.argv[1]) as ds1:
    with Dataset(sys.argv[2], "w") as ds2:
        copy_variable(ds1, ds2, "topg", "y")
        copy_variable(ds1, ds2, "mapping", None)
EOF

chmod u+x ./invert_y.py

./invert_y.py ${regular} ${inverted}

common_options="
-config ${PISM_PATH}/pism_config.nc
-atmosphere uniform
-bootstrap
-energy none
-i ${regular}
-no_mass
-regrid_vars topg
-stress_balance none
-surface simple
-y 1s
"

${pism} ${common_options} -regrid_file ${regular}  -o regular_topg.nc
${pism} ${common_options} -regrid_file ${inverted} -o inverted_topg.nc

${PISM_PATH}/pism_nccmp -v topg regular_topg.nc inverted_topg.nc
