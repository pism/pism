# Run "make all" to update lists of diagnostics and configuration parameters.

all: update_diagnostics update_parameters

.INTERMEDIATE: diagnostics.json dummy.nc pism_config.nc

dummy.nc:
	# use a thermal bedrock layer model and run for 1 year to complete the
	# bootstrapping
	pisms -y 1 -Mbz 5 -Lbz 1000 -Mz 5 -Mx 5 -My 5 -o dummy.nc -verbose 1

	# rename climatic_mass_balance: we need precipitation for an atmosphere model
	ncrename -v climatic_mass_balance,precipitation -O dummy.nc dummy.nc

	# create dummy firn and snow depth fields
	ncap2 -s "firn_depth=thk" -O dummy.nc dummy.nc
	ncatted -a standard_name,firn_depth,d,, -O dummy.nc dummy.nc
	ncap2 -s "snow_depth=firn_depth" -O dummy.nc dummy.nc
	# create a dummy (zero) basal water thickness field
	ncap2 -s "bwat=firn_depth" -O dummy.nc dummy.nc

diagnostics.json: dummy.nc
	pismr -i $^ \
	      -energy enthalpy \
	      -bed_def iso \
	      -surface pdd \
	      -atmosphere searise_greenland \
	      -stress_balance ssa+sia \
	      -yield_stress mohr_coulomb \
	      -topg_to_phi 10,30,-50,0 \
	      -hydrology routing \
	      -calving eigen_calving,thickness_calving,frontal_melt,ocean_kill,vonmises_calving \
	      -ocean_kill_file dummy.nc \
	      -verbose 1 -o_size none \
	      -list_diagnostics_json > $@

update_diagnostics: diagnostics.json
	python ./list_diagnostics.py $^ > manual/diagnostics/diagnostics-list.txt

pism_config.nc:
	ncgen -o $@ ../../src/pism_config.cdl

update_parameters: pism_config.nc
	python ./list_parameters.py $^ > manual/parameters/parameter-list.txt
