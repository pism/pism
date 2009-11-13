ALL: userman

revision:
	@svnversion > revision.tex

revclean:
	@rm -f revision.tex

userman : manual.tex ../ice_bib.bib revision
	pdflatex manual.tex
	BIBINPUTS=".:..:" bibtex manual
	makeindex -l manual
	pdflatex manual.tex
	pdflatex manual.tex

clean : revclean
	@rm -f *.aux *.blg *.idx *.ind *.ilg *.log *.out *.toc texput.* *.bbl
