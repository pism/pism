all: pism ismiphom plot

pism:
	mpiexec -n 2 python3 run-ismiphom.py -Mx 101 -My 101 -bp_snes_monitor -bp_ksp_monitor -bp_pc_type mg -bp_pc_mg_levels 2 -stress_balance.blatter.coarsening_factor 4 -bp_mg_coarse_pc_type gamg

ismiphom: ismip_all
	python3 convert-ismiphom.py

tc-2-95-2008-supplement.zip:
	wget -nc "https://tc.copernicus.org/articles/2/95/2008/tc-2-95-2008-supplement.zip"

ismip_all: tc-2-95-2008-supplement.zip
	unzip -n $^
	unzip -n tc-2007-0019-sp2.zip

plot:
	python3 plot-ismiphom.py
