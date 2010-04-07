ALL: userman

parts := eismint-greenland.tex eismint-ross.tex getting-started.tex \
	 modeling-choices.tex practical-usage.tex simplified-geometry.tex \
	 verification.tex

revision:
	@svnversion > revision.tex

revclean:
	@rm -f revision.tex

# creates a symlink to ice_bib.bib with a different name (so that we can avoid
# underscores in filenames included in the manual)
ice-bib.bib:
	ln -s ../ice_bib.bib ice-bib.bib

userman : manual.tex ice-bib.bib revision $(parts)
	pdflatex manual.tex
	BIBINPUTS=".:..:" bibtex manual
	makeindex -l general
	printf ",s/\\_ /\\_/g\n,wq\n" | ed general.ind 
	makeindex -l options
	printf ",s/\\_ /\\_/g\n,wq\n" | ed options.ind
	pdflatex manual.tex
	pdflatex manual.tex

clean : revclean
	@rm -f *.aux *.blg *.idx *.ind *.ilg *.log *.out *.toc texput.* *.bbl
