#!/bin/bash

H0=800
v0=300
ca=1750e3
N0=601

ln -s ../mismip2d/MISMIP.py
ln -s ../../preprocessing/PISMNC.py

for per in "default" "p1" "p2" "p3"
do
  rf=result_1a_A1_${H0}_${v0}_${ca}_${N0}_${per}.nc
  of=ICESHELF_1a_A1_${H0}_${v0}_${ca}_${N0}_${per}.nc
  python3 prepare_iceshelf.py -c $ca -v $v0 -H $H0 -o $of -N $N0 -p $per

  pism \
    -Lz 6000 \
    -Mx $N0 \
    -My 3 \
    -Mz 3 \
    -bootstrap \
    -cfbc \
    -energy none \
    -i $of \
    -no_mass \
    -o $rf -o_order zyx -o_size big \
    -output.sizes.big taud \
    -ssa_dirichlet_bc \
    -ssa_flow_law isothermal_glen \
    -ssa_method fd \
    -stress_balance ssa \
    -stress_balance.ssa.fd.flow_line_mode \
    -y 1s \
    > pism-${per}.log
    # -ssafd_ksp_rtol 1e-7 \

  echo $rf $of $addopt

done

python3 plot_diff.py
