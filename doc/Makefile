ALL : manual.pdf docs

docs : refman browser

refman : Doxyfile_tex
	doxygen Doxyfile_tex
	(cd doxy/latex/; ${MAKE} all)
	cp doxy/latex/refman.pdf .

browser : Doxyfile_html
	doxygen Doxyfile_html

vfigs : vfnow.txt
	# producing figures thickerrs_C.png, thickerrs_G.png, temperrs_G.png, surfvelerrs_G.png, 
	# velerrs_I.png, and temperrs_K.png in subdirectory figs/ for use in manual
	vnreport.py -f vfnow.txt -t C -e 'geometry' -x 'prcntVOLrelmaxETA' -o figs/thickerrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/thickerrs_C.png figs/thickerrs_C.eps
	rm figs/thickerrs_C.eps
	vnreport.py -f vfnow.txt -t G -e 'geometry' -x 'prcntVOLrelmaxETA' -o figs/thickerrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/thickerrs_G.png figs/thickerrs_G.eps
	rm figs/thickerrs_G.eps
	vnreport.py -f vfnow.txt -t G -e 'temp' -x 'basecenterT' -o figs/temperrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/temperrs_G.png figs/temperrs_G.eps
	rm figs/temperrs_G.eps
	vnreport.py -f vfnow.txt -t G -e 'surf vels' -x 'maxUvecmaxW' -o figs/surfvelerrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/surfvelerrs_G.png figs/surfvelerrs_G.eps
	rm figs/surfvelerrs_G.eps
	vnreport.py -f vfnow.txt -t I -e 'velocity' -x 'maxvectorprcntavvec' -o figs/velerrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/velerrs_I.png figs/velerrs_I.eps
	rm figs/velerrs_I.eps
	vnreport.py -f vfnow.txt -t K -e 'temp' -o figs/temperrs
	gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 \
	   -sOutputFile=figs/temperrs_K.png figs/temperrs_K.eps
	rm figs/temperrs_K.eps

%.pdf : %.tex
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

clean :
	rm -f manual.aux manual.bbl manual.blg manual.dvi manual.idx manual.log manual.out manual.ps manual.toc

