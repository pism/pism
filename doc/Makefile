ALL : manual.dvi manual.ps eqns3D.dvi eqns3D.ps

%.dvi : %.tex
	latex $*
	bibtex $*
	latex $*
	latex $*

%.ps : %.dvi
	dvips $*

clean :
	rm -f {manual,eqns3D}.{aux,bbl,blg,log,toc,dvi,ps}
