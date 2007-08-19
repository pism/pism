ALL : manual.pdf docs

docs : refman browser

refman : Doxyfile_tex
	doxygen Doxyfile_tex
	(cd doxy/latex/; ${MAKE} all)
	cp doxy/latex/refman.pdf .

browser : Doxyfile_html
	doxygen Doxyfile_html

%.pdf : %.tex
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

#%.dvi : %.tex
#	latex $*
#	bibtex $*
#	latex $*
#	latex $*

#%.ps : %.dvi
#	dvips $*

clean :
	rm -f manual.aux manual.bbl manual.blg manual.dvi manual.idx manual.log manual.out manual.ps manual.toc
#	rm -f eqns3D.aux eqns3D.bbl eqns3D.blg eqns3D.dvi eqns3D.log eqns3D.out eqns3D.ps eqns3D.toc

