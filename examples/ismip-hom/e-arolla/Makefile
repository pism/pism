all: figures

output-sliding.nc: run.sh arolla100-sliding.nc
	bash $^ $@

output-no-slip.nc: run.sh arolla100-no-slip.nc
	bash $^ $@

arolla100.dat:
	wget -nc "http://homepages.ulb.ac.be/~fpattyn/ismip/arolla100.dat"

arolla100-sliding.nc arolla100-no-slip.nc: arolla100.dat convert.py
	python3 convert.py

clean:
	@rm -f *.nc

figures: output-no-slip.nc output-sliding.nc
	python3 plot.py
