all: isochrones.mpg

isochrones.mpg: ex.nc plot.py
	python3 plot.py -o $@ ex.nc

isochrones.png: ex.nc output.nc plot.py
	python3 plot.py -f -p output.nc -o $@ ex.nc

ex.nc output.nc: run.sh
	bash -x $^

clean:
	@rm -f *.nc

ex_part1.nc ex_part2.nc o_part1.nc o_final.nc: run_with_restart.sh
	bash -x $^

isochrones_with_restart.mpg: ex_part1.nc ex_part2.nc
	python3 plot.py -o $@ $^

ex_regrid_part1.nc ex_regrid_part2.nc o_regrid_part1.nc o_regrid_part2.nc: run_with_regrid.sh
	bash -x $^

isochrones_with_regrid.mpg: ex_regrid_part1.nc ex_regrid_part2.nc
	python3 plot.py -o $@ $^

ex_bootstrap.nc: run_with_bootstrap.sh
	bash -x $^

isochrones_with_bootstrap.mpg: ex_bootstrap.nc
	python3 plot.py -o $@ $^
