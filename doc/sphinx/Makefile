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

github: Makefile clean
	mkdir -p ${OUTDIR}	# create the output directory

	ln -s $(abspath ${OUTDIR}) ${BUILDDIR}/html # create a link to make Sphinx use it

	${MAKE} html

github_clean:
	rm -rf ${BUILDDIR}/html	    # remove the link
	rm -rf $(abspath ${OUTDIR}) # clean out the output directory

clean: Makefile github_clean

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
