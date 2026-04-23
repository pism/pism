all: isochrones.mpg

isochrones.mpg: spatial.nc plot.py
	python3 plot.py -o $@ spatial.nc

isochrones.png: spatial.nc output.nc plot.py
	python3 plot.py -f -p output.nc -o $@ spatial.nc

spatial.nc output.nc: run.sh
	bash -x $^

clean:
	@rm -f *.nc

spatial_part1.nc spatial_part2.nc o_part1.nc o_final.nc: run_with_restart.sh
	bash -x $^

isochrones_with_restart.mpg: spatial_part1.nc spatial_part2.nc
	python3 plot.py -o $@ $^

spatial_regrid_part1.nc spatial_regrid_part2.nc o_regrid_part1.nc o_regrid_part2.nc: run_with_regrid.sh
	bash -x $^

isochrones_with_regrid.mpg: spatial_regrid_part1.nc spatial_regrid_part2.nc
	python3 plot.py -o $@ $^

spatial_bootstrap.nc: run_with_bootstrap.sh
	bash -x $^

isochrones_with_bootstrap.mpg: spatial_bootstrap.nc
	python3 plot.py -o $@ $^
