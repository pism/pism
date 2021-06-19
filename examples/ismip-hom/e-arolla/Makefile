all: figures

clean:
	@rm -f *.nc

arolla100.dat:
	wget -nc "http://homepages.ulb.ac.be/~fpattyn/ismip/arolla100.dat"

arolla100-%.nc: convert.py arolla100.dat
	python3 $^

figures: output-sliding-401.nc output-no-slip-401.nc
	python3 plot.py $^

output-sliding-%.nc: run.sh arolla100-sliding.nc
	bash $^ $@ $* 4
	$(extract-flowline)

output-no-slip-%.nc: run.sh arolla100-no-slip.nc
	bash $^ $@ $* 4
	$(extract-flowline)

output-z-refinement-%.nc: run.sh arolla100-no-slip.nc
	bash $^ $@ 401 $*
	$(extract-flowline)

.PHONY: x-refinement-no-slip x-refinement-sliding z-refinement
Mxs=51 101 201 401 801
x-refinement-no-slip: $(foreach Mx,${Mxs},output-no-slip-${Mx}.nc)

x-refinement-sliding: $(foreach Mx,${Mxs},output-sliding-${Mx}.nc)

Ms=2 3 4 5
z-refinement: $(foreach M,${Ms},output-z-refinement-${M}.nc)

define extract-flowline =
ncks -O -d y,1 $@ $@
ncwa -O -a time,y $@ $@
ncpdq -O -a z,x $@ $@
ncpdq -O -a z_sigma,x $@ $@
endef
