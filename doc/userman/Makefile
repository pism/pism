ALL: userman

parts := eismint-greenland.tex eismint-ross.tex getting-started.tex \
	 modeling-choices.tex practical-usage.tex simplified-geometry.tex \
	 verification.tex

revision:
	@svnversion > revision.tex

revclean:
	@rm -f revision.tex

userman : manual.tex ../ice_bib.bib revision $(parts)
	pdflatex manual.tex
	BIBINPUTS=".:..:" bibtex manual
	makeindex -l manual
	pdflatex manual.tex
	pdflatex manual.tex

clean : revclean
	@rm -f *.aux *.blg *.idx *.ind *.ilg *.log *.out *.toc texput.* *.bbl
