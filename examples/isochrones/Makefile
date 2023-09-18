all: isochrones.mpg

isochrones.mpg: ex.nc
	python3 plot.py $^ $@

ex.nc output.nc: run.sh
	bash -x $^

clean:
	@rm -f *.nc

ex_with_restart.nc intermediate_result.nc final_result.nc: run_with_restart.sh
	bash -x $^

isochrones_with_restart.mpg: ex_with_restart.nc
	python3 plot.py $^ $@
