#!/bin/bash

# these figures show velocity bins, with coloring by ice thickness,
# where P(W) looks quite different in each bin
# each bin is managably small for even 2km Greenland run, I think
# run as:
#   $ ./genscatfig.sh ex_distributed-decoupled.nc g2km.png
# to generate figure files like
#   bin*-g2km.png

FILENAME=$1
OUTROOT=$2

./showPvsW.py -wmin 0.0 -wmax 0.15 -c thk -cmin 0 -cmax 2000 -s hydrovelbase_mag -smin 300  -smax 600 -o bin300-${OUTROOT} --colorbar $FILENAME
./showPvsW.py -wmin 0.0 -wmax 0.15 -c thk -cmin 0 -cmax 2000 -s hydrovelbase_mag -smin 100  -smax 200  -o bin100-${OUTROOT}  $FILENAME
./showPvsW.py -wmin 0.0 -wmax 0.15 -c thk -cmin 0 -cmax 2000 -s hydrovelbase_mag -smin 30   -smax 60  -o bin30-${OUTROOT}   $FILENAME
./showPvsW.py -wmin 0.0 -wmax 0.15 -c thk -cmin 0 -cmax 2000 -s hydrovelbase_mag -smin 10   -smax 20   -o bin10-${OUTROOT}   $FILENAME
./showPvsW.py -wmin 0.0 -wmax 0.15 -c thk -cmin 0 -cmax 2000 -s hydrovelbase_mag -smin 3    -smax 6   -o bin1-${OUTROOT}    $FILENAME

