ALL : userman summary refman browser

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


PISMREV ?= 205 # revision stamp for verification results in User Manual figures;
               # generally not the current trunk revision
VFNOW_OUT=vfnow.txt
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
	rm -f manual.aux manual.bbl manual.blg manual.dvi manual.idx manual.ind \
	      manual.ilg  manual.log manual.out manual.ps manual.toc
	rm -f fullbib.aux fullbib.log fullbib.blg fullbib.bbl fullbib.out
	rm -f pism_summary.aux pism_summary.log pism_summary.blg pism_summary.bbl \
	      pism_summary.out
	

