TS=-ts_times 1 -ts_file
GRID= -Mx 201 -My 3 -Mz 3 -Lz 1000
OTHER=-stress_balance prescribed_sliding -prescribed_sliding_file input.nc -o_order zyx -energy none -no_sia -shelf_base_melt_rate 0.0 -max_dt 0.123456789 -float_kill
EXTRA=-extra_times 1 -extra_vars thk,mask,usurf,topg -extra_file
RUN_LENGTH=800

OPTIONS=${GRID} -y ${RUN_LENGTH} ${OTHER}
input.nc: generate_input.py
	python ./generate_input.py -M 201 -o $@

default: input.nc
	pismr -boot_file input.nc ${OPTIONS} ${TS} ts_$@.nc ${EXTRA} ex_$@.nc

part_grid: input.nc
	pismr -boot_file input.nc ${OPTIONS} ${TS} ts_$@.nc ${EXTRA} ex_$@.nc -part_grid

part_redist: input.nc
	pismr -boot_file input.nc ${OPTIONS} ${TS} ts_$@.nc ${EXTRA} ex_$@.nc -part_grid -part_redist

clean:
	@rm -f *.nc*
