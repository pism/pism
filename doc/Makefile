ALL : installation userman refman browser

update: svn_update
ifeq ($(shell touch revision.tex; svnversion | diff revision.tex -),)
	@echo "doc/ directory is up to date."
else
	@echo "Rebuilding manuals ..."
	@svnversion > revision.tex
	$(MAKE) ALL
endif

# manuals:

installation : installation.tex
	pdflatex installation.tex
	pdflatex installation.tex

userman : manual.tex ice_bib.bib
	pdflatex manual.tex
	bibtex manual
	makeindex -l manual
	pdflatex manual.tex
	pdflatex manual.tex

refman : Doxyfile_tex
	doxygen Doxyfile_tex
	(cd doxy/latex/; ${MAKE} all)
	cp doxy/latex/refman.pdf .

browser : Doxyfile_html
	doxygen Doxyfile_html


# not default targets; generally not needed by users:

summary : pism_summary.tex
	pdflatex pism_summary
	bibtex pism_summary
	pdflatex pism_summary
	pdflatex pism_summary

fullbib : fullbib.tex
	pdflatex fullbib
	bibtex fullbib
	pdflatex fullbib
	pdflatex fullbib


PISMREV ?= 327 # revision stamp for verification results in User's Manual figures;
               # generally not the current trunk revision
VFNOW_OUT=vfnow_r327_L5_23june08.txt
vfigs : $(VFNOW_OUT)
	# producing figures thickerrs_B.pdf, thickerrs_G.df, etc. 
	# in subdirectory figs/ for use in manual
	../test/vnreport.py -f $(VFNOW_OUT) -t B -e 'geometry' -x 'prcntVOLrelmaxETAcenterH' \
	   -r ${PISMREV} -g pdf -o figs/thickerrs
	../test/vnreport.py -f $(VFNOW_OUT) -t G -e 'geometry' -x 'prcntVOLrelmaxETAcenterH' \
	   -r ${PISMREV} -g pdf -o figs/thickerrs
	../test/vnreport.py -f $(VFNOW_OUT) -t G -e 'temp' -x 'basecenterT' \
	   -r ${PISMREV} -g pdf -o figs/temperrs
	../test/vnreport.py -f $(VFNOW_OUT) -t G -e 'surf vels' -x 'maxUvecmaxW' \
	   -r ${PISMREV} -g pdf -o figs/surfvelerrs
	../test/vnreport.py -f $(VFNOW_OUT) -t I -e 'velocity' -x 'maxvectorprcntavvecavvector' \
	   -r ${PISMREV} -g pdf -o figs/velerrs
	../test/vnreport.py -f $(VFNOW_OUT) -t K -e 'temp' -r ${PISMREV} -g pdf -o figs/temperrs

clean :
	@rm -f *.aux *.blg *.idx *.ind *.ilg *.log *.out *.toc
	@rm -f manual.bbl refman.bbl pism_summary.bbl  # I like to keep fullbib.bbl

svn_update:
	@echo "Running 'svn update' ($(shell svn info |grep 'Repository Root'))..."
	@svn update

