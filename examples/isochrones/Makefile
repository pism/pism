all: isochrones.mpg

isochrones.mpg: ex.nc
	python3 plot.py $^ $@

ex.nc output.nc: run.sh
	bash -x $^

clean:
	@rm -f *.nc
