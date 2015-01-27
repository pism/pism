#!/bin/bash

# generate all figures

INITROOT=$1

./genGreenfig.sh $INITROOT Greenland_5km_v1.1
mv -v $INITROOT-velbase_mag.png $INITROOT-velbase-mag.png
mv -v $INITROOT-velsurf_mag.png $INITROOT-velsurf-mag.png
mv -v Greenland_5km_v1.1-surfvelmag.png Greenland-surfvelmag.png

./basemapfigs.py routing-decoupled tillwat
./basemapfigs.py routing-decoupled bwat

./basemapfigs.py distributed-decoupled tillwat
./basemapfigs.py distributed-decoupled bwat
./basemapfigs.py distributed-decoupled bwprel

./genscatfig.sh ex_distributed-decoupled.nc $INITROOT.png

rm -rf listpng.txt
ls *.png > listpng.txt

for name in `cat listpng.txt`; do
  echo "autocropping ${name} ..."
  mogrify -trim +repage $name
done
rm -rf listpng.txt

