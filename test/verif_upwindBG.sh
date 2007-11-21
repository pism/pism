#!/bin/bash
# verification of tests B and G using upwind method for SIA; use in revision 210 or later

NN=8

echo mpiexec -np $NN pismv -test B -Mx 31 -My 31 -Mz 31 -ys 422.45 -y 25000.0 -upwind -max_dt 16.0
mpiexec -np $NN pismv -test B -Mx 31 -My 31 -Mz 31 -ys 422.45 -y 25000.0 -verbose 1 -upwind -max_dt 16.0

echo mpiexec -np $NN pismv -test B -Mx 41 -My 41 -Mz 31 -ys 422.45 -y 25000.0 -upwind -max_dt 8.0
mpiexec -np $NN pismv -test B -Mx 41 -My 41 -Mz 31 -ys 422.45 -y 25000.0 -verbose 1 -upwind -max_dt 8.0

echo mpiexec -np $NN pismv -test B -Mx 61 -My 61 -Mz 31 -ys 422.45 -y 25000.0 -upwind -max_dt 4.0
mpiexec -np $NN pismv -test B -Mx 61 -My 61 -Mz 31 -ys 422.45 -y 25000.0 -verbose 1 -upwind -max_dt 4.0

echo mpiexec -np $NN pismv -test B -Mx 81 -My 81 -Mz 31 -ys 422.45 -y 25000.0 -upwind -max_dt 2.3
mpiexec -np $NN pismv -test B -Mx 81 -My 81 -Mz 31 -ys 422.45 -y 25000.0 -verbose 1 -upwind -max_dt 2.3

echo mpiexec -np $NN pismv -test B -Mx 121 -My 121 -Mz 31 -ys 422.45 -y 25000.0 -upwind -max_dt 1.0
mpiexec -np $NN pismv -test B -Mx 121 -My 121 -Mz 31 -ys 422.45 -y 25000.0 -verbose 1 -upwind -max_dt 1.0


echo mpiexec -np $NN pismv -test G -Mx 61 -My 61 -Mz 61 -y 25000.0 -upwind -max_dt 22.0
mpiexec -np $NN pismv -test G -Mx 61 -My 61 -Mz 61 -y 25000.0 -verbose 1 -upwind -max_dt 22.0

echo mpiexec -np $NN pismv -test G -Mx 91 -My 91 -Mz 91 -y 25000.0 -upwind -max_dt 9.9
mpiexec -np $NN pismv -test G -Mx 91 -My 91 -Mz 91 -y 25000.0 -verbose 1 -upwind -max_dt 9.9

echo mpiexec -np $NN pismv -test G -Mx 121 -My 121 -Mz 121 -y 25000.0 -upwind -max_dt 5.6
mpiexec -np $NN pismv -test G -Mx 121 -My 121 -Mz 121 -y 25000.0 -verbose 1 -upwind -max_dt 5.6

echo mpiexec -np $NN pismv -test G -Mx 181 -My 181 -Mz 181 -y 25000.0 -upwind -max_dt 2.5
mpiexec -np $NN pismv -test G -Mx 181 -My 181 -Mz 181 -y 25000.0 -verbose 1 -upwind -max_dt 2.5

echo mpiexec -np $NN pismv -test G -Mx 241 -My 241 -Mz 241 -y 25000.0 -upwind -max_dt 1.4
mpiexec -np $NN pismv -test G -Mx 241 -My 241 -Mz 241 -y 25000.0 -verbose 1 -upwind -max_dt 1.4


