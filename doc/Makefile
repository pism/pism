ALL : userman summary refman browser

userman : manual.tex ice_bib.bib figs/*.png
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


# this variable is the revision stamp for verification results figures which appear
#   in the users manual; it is generally not the current trunk revision
PISMREV ?= 205

vfigs : vfnow.txt
	# producing figures thickerrs_B.png, thickerrs_G.png, temperrs_G.png, surfvelerrs_G.png, 
	# velerrs_I.png, and temperrs_K.png in subdirectory figs/ for use in manual
	vnreport.py -f vfnow.txt -t B -e 'geometry' -x 'prcntVOLrelmaxETAcenterH' -r ${PISMREV} \
	   -o figs/thickerrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/thickerrs_B.png figs/thickerrs_B.eps
	rm figs/thickerrs_B.eps
	vnreport.py -f vfnow.txt -t G -e 'geometry' -x 'prcntVOLrelmaxETAcenterH' -r ${PISMREV} \
	   -o figs/thickerrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/thickerrs_G.png figs/thickerrs_G.eps
	rm figs/thickerrs_G.eps
	vnreport.py -f vfnow.txt -t G -e 'temp' -x 'basecenterT' -r ${PISMREV} -o figs/temperrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/temperrs_G.png figs/temperrs_G.eps
	rm figs/temperrs_G.eps
	vnreport.py -f vfnow.txt -t G -e 'surf vels' -x 'maxUvecmaxW' -r ${PISMREV} \
	   -o figs/surfvelerrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/surfvelerrs_G.png figs/surfvelerrs_G.eps
	rm figs/surfvelerrs_G.eps
	vnreport.py -f vfnow.txt -t I -e 'velocity' -x 'maxvectorprcntavvecavvector' -r ${PISMREV} \
	   -o figs/velerrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/velerrs_I.png figs/velerrs_I.eps
	rm figs/velerrs_I.eps
	vnreport.py -f vfnow.txt -t K -e 'temp' -r ${PISMREV} -o figs/temperrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/temperrs_K.png figs/temperrs_K.eps
	rm figs/temperrs_K.eps


clean :
	rm -f manual.aux manual.bbl manual.blg manual.dvi manual.idx manual.ind manual.ilg \
	      manual.log manual.out manual.ps manual.toc
	rm -f fullbib.aux fullbib.log fullbib.blg fullbib.bbl
	rm -f pism_summary.aux pism_summary.log pism_summary.blg pism_summary.bbl pism_summary.out
	

