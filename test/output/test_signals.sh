#!/bin/bash

set -u
set -x
set -e

eis2_grid="-Mx 60 -My 60 -Mz 11 -Lz 5000"

pism -eisII A -y 1e6 ${eis2_grid} -o test_signals.nc &

pid=$!

sleep 1

kill -s SIGUSR1 ${pid}

sleep 1

kill -s SIGUSR2 ${pid}

sleep 1

kill -s SIGTERM ${pid}

wait ${pid}
