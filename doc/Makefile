ALL : manual.pdf eqns3D.dvi eqns3D.ps

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
	rm -f {manual,eqns3D}.{aux,bbl,blg,log,toc,dvi,ps,pdf}
