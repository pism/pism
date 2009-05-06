ALL : installation userman browser refman

update: svn_update
ifeq ($(shell touch revision.tex; svnversion | diff revision.tex -),)
	@echo "doc/ directory is up to date."
else
	@echo "Rebuilding manuals ..."
	@svnversion > revision.tex
	$(MAKE) ALL
endif

# default manuals and browser:

installation : installation.tex
	pdflatex installation.tex
	pdflatex installation.tex

userman : manual.tex ice_bib.bib
	pdflatex manual.tex
	bibtex manual
	makeindex -l manual
	pdflatex manual.tex
	pdflatex manual.tex

refman : doxyfile_refman
	doxygen doxyfile_refman
	(cd doxy/latex/; ${MAKE} all)
	cp doxy/latex/refman.pdf .

browser : doxyfile_browser doxybib.bbl config_doc.txt
	doxygen doxyfile_browser

doxybib.bbl: ice_bib.bib
	./doxybib.py

config_doc.txt: ../lib/pism_config.nc
	./config_doc.py > config_doc.txt

# non-default targets; sometimes useful:
#   browsenow      rebuild browser and then open it with firefox
#   summary        build a PDF summarizing PISM features; not always up to date
#   fullbib        build a PDF which merely lists every reference in ice_bib.bib
#   vfigs          create PDF figures from a test/verifynow.py output

browsenow : browser
	firefox doxy/html/index.html &

summary : pism_summary.tex
	pdflatex pism_summary
	bibtex pism_summary
	pdflatex pism_summary
	pdflatex pism_summary


PISMREV_VFIGS ?= 589 # revision stamp for verification results in User's Manual figures;
                     # generally not the current trunk revision
VFNOW_OUT=vfnow_r589_L5_5mar09.txt
# list old ones here; keep them, for comparison:
#VFNOW_OUT=vfnow_r327_L5_23jun08.txt
#VFNOW_OUT=vfnow_r193_L5_10oct07.txt
vfigs : verif/$(VFNOW_OUT)
	# producing figures thickerrs_G.pdf, temperrs_G.df, etc. 
	# in subdirectory figs/ for use in manual
	../test/vnreport.py -f verif/$(VFNOW_OUT) -t G -e 'geometry' -x 'prcntVOLrelmaxETAcenterH' \
	   -r ${PISMREV_VFIGS} -g pdf -o figs/thickerrs
	../test/vnreport.py -f verif/$(VFNOW_OUT) -t G -e 'temp' -x 'basecenterT' \
	   -r ${PISMREV_VFIGS} -g pdf -o figs/temperrs
	../test/vnreport.py -f verif/$(VFNOW_OUT) -t G -e 'surf vels' -x 'maxUvecmaxW' \
	   -r ${PISMREV_VFIGS} -g pdf -o figs/surfvelerrs
	../test/vnreport.py -f verif/$(VFNOW_OUT) -t I -e 'velocity' -x 'maxvectorprcntavvecavvector' \
	   -r ${PISMREV_VFIGS} -g pdf -o figs/velerrs
	../test/vnreport.py -f verif/$(VFNOW_OUT) -t K -e 'temp' \
	   -r ${PISMREV_VFIGS} -g pdf -o figs/temperrs

clean :
	@rm -f *.aux *.blg *.idx *.ind *.ilg *.log *.out *.toc texput.* *.bbl
	@rm -f config_doc.txt

svn_update:
	@echo "Running 'svn update' ($(shell svn info |grep 'Repository Root'))..."
	@svn update

