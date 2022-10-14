TS=-scalar_times 1 -scalar_file
GRID= -Mx 201 -My 3 -Mz 3 -Lz 1000
OTHER=-stress_balance prescribed_sliding -prescribed_sliding_file input.nc -o_order zyx -energy none -shelf_base_melt_rate 0.0 -max_dt 0.123456789 -calving float_kill
EXTRA=-spatial_times 1 -spatial_vars thk,mask,usurf,topg,Href -spatial_file
RUN_LENGTH=800

OPTIONS=${GRID} -y ${RUN_LENGTH} ${OTHER}
input.nc: generate_input.py
	python ./generate_input.py -M 201 -o $@

run: input.nc
	pismr -i input.nc -bootstrap ${OPTIONS} ${TS} ts_$@.nc ${EXTRA} ex_$@.nc

clean:
	@rm -f *.nc*
