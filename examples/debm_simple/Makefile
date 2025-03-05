all: present_day.nc paleo.nc with_albedo.nc

present_day.nc: pism_Greenland_5km_v1.1.nc
	present_day=1 ./run.sh $@ $^

with_albedo.nc: pism_Greenland_5km_v1.1.nc present_day.nc
	albedo=1 ./run.sh $@ $^

paleo.nc: pism_Greenland_5km_v1.1.nc orbital_parameters.nc
	paleo=1 ./run.sh $@ $^

pism_Greenland_5km_v1.1.nc:
	cd ../std-greenland/ && ./preprocess.sh
	ln -v -s ../std-greenland/$@ .

orbital_parameters.nc: orbital_parameters.py
	python3 $< -o $@
	ncap2 -O -s 'defdim("nv",2);time_bnds=make_bounds(time,$$nv,"time_bnds");' $@ $@

.PHONY: clean very_clean

clean:
	rm -f present_day.nc with_albedo.nc paleo.nc

very_clean: clean
	rm -f orbital_parameters.nc pism_Greenland_5km_v1.1.nc

.INTERMEDIATE: INSOLN.LA2004.BTL.ASC
INSOLN.LA2004.BTL.ASC:
	wget "https://vo.imcce.fr/insola/earth/online/earth/La2004/INSOLN.LA2004.BTL.ASC"

# Convert from Fortran to C/C++ floating point format:
laskar2004_insoln.txt: INSOLN.LA2004.BTL.ASC
	sed "s/D+/E+/g;s/D-/E-/g" $^ > $@

berger1978_vs_laskar2004.png: laskar2004_insoln.txt
	python3 ./compare.py $^ $@
