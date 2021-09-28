
H0=800
v0=300
ca=1750e3
N0=601
#N0=1801
#N0=3001

ln -s ../mismip2d/MISMIP.py
ln -s ../../preprocessing/PISMNC.py
#python get_conf.py -e 1a --mode=1

for scheme in "onesided" "centered" 
do
for per in "default" "p1" "p2" "p3"
do
rf=result_1a_A1_${H0}_${v0}_${ca}_${N0}_${scheme}_${per}.nc
of=ICESHELF_1a_A1_${H0}_${v0}_${ca}_${N0}_${per}.nc
python prepare_iceshelf.py -c $ca -v $v0 -H $H0 -o $of -N $N0 -p $per 2>/dev/null

addopt=""
if [[ "$scheme" == "centered" ]]; then
  addopt="-taudfix"
fi

/home/albrecht/pism21/pism-dev/bin/pismr -i $of -bootstrap -Mx $N0 -My 3 -Mz 15 -Lz 6000 \
-energy none -no_mass -shelf_base_melt_rate 0.0 \
-stress_balance ssa -ssa_dirichlet_bc -ssa_flow_law isothermal_glen \
-ssa_method fd -cfbc -part_grid -ssafd_ksp_rtol 1e-7 \
-yield_stress constant -tauc 7.624000e+06 -pseudo_plastic -pseudo_plastic_q 3.333333e-01 -pseudo_plastic_uthreshold 3.155693e+07 \
-front_retreat_file $of \
-config_override MISMIP_conf_1a_A1.nc \
-ys 0 -ye 0.0001 -options_left \
-o $rf -o_order zyx -o_size big \
$addopt >> pism.out

echo $rf $of $addopt
#python plot_iceshelf.py $rf $of
done
done

python plot_diff.py


