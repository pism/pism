all: pism ismiphom plot

pism:
	mpiexec -n 2 python3 run-ismiphom.py -Mx 101 -My 101 -snes_monitor -ksp_monitor -pc_type mg -pc_mg_levels 2 -stress_balance.blatter.coarsening_factor 4 -mg_coarse_pc_type gamg

ismiphom:
	python3 convert-ismiphom.py

plot:
	python3 plot-ismiphom.py
