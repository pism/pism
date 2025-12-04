#!/bin/bash

set -u
set -x
set -e

# adjust grid size and -output.checkpoint.interval so that it saves a checkpoint file 2 or
# 3 times
pism \
  -eisII A \
  -grid.Lz 5000 \
  -grid.Mx 40 \
  -grid.My 40 \
  -grid.Mz 21 \
  -output.checkpoint.interval 5s \
  -output.file eis2-a.nc \
  -time.run_length 1e4 \
  ""
