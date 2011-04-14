all:
	maxima -b ssa_prettyprint.mac
	maxima -b ssa_code.mac
	pdflatex ssa-discretization.tex
	pdflatex ssa-discretization.tex

clean:
	rm -f *.pdf *.aux *.log *.rel
