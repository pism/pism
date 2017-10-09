# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = PISM
SOURCEDIR     = .
BUILDDIR      = _build

OUTDIR = ./docs			# for publishing on GitHub pages

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Targets for publishing on GitHub Pages

github: Makefile clean
	mkdir -p ${OUTDIR}	# create the output directory

	ln -s $(abspath ${OUTDIR}) ${BUILDDIR}/html # create a link to make Sphinx use it

	${MAKE} html

github_clean:
	rm -rf ${BUILDDIR}/html	    # remove the link
	rm -rf $(abspath ${OUTDIR}) # clean out the output directory

clean: Makefile github_clean

# Add the "latex" tag when building PDFs

latex%: Makefile
	echo "Building LaTeX stuff"
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) -t latex

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Targets for updating the list of diagnostics:

.INTERMEDIATE: diagnostics.json dummy.nc

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

PISM_CONFIG=~/local/pism/share/pism/pism_config.nc

update_parameters:
	python ./list_parameters.py ${PISM_CONFIG} > manual/parameters/parameter-list.txt
