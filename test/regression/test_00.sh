#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test # 0: presence of tools and Python modules needed by other tests."
# The list of files to delete when done.
files="foo.txt"

rm -f $files

tools="ncpdq ncap2 ncks ncrename diff rm cp"
modules="numpy sys"
success=1

for tool in $tools;
do
    if [ -z `which $tool` ];
    then
        echo "ERROR: $tool is not on the path"
        success=0
    fi
done

for module in $modules;
do
    echo "import $module" | /usr/bin/env python &> foo.txt
    if [ -s foo.txt ];
    then
        echo "ERROR: Python module \"$module\" is not installed"
        success=0
    fi
done

/usr/bin/env python &> foo.txt <<EOF
try:
    import netCDF3
except:
    import netCDF4
EOF

if [ -s foo.txt ];
then
    echo "ERROR: neither netCDF3 nor netCDF4 module is installed"
    success=0
fi

if [ $success == 0 ];
then
    rm -f $files
    exit 1
fi

rm -f $files; exit 0
