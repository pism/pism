ALL : manual.pdf

%.pdf : %.tex
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

%.dvi : %.tex
	latex $*
	bibtex $*
	latex $*
	latex $*

%.ps : %.dvi
	dvips $*

clean :
	rm -f manual.aux manual.bbl manual.blg manual.dvi manual.log manual.out manual.ps manual.toc
	rm -f eqns3D.aux eqns3D.bbl eqns3D.blg eqns3D.dvi eqns3D.log eqns3D.out eqns3D.ps eqns3D.toc

