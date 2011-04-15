all: pdf

pdf: formulas.tex
	pdflatex ssa-discretization.tex
	pdflatex ssa-discretization.tex

formulas.tex:
	maxima -b ssa_prettyprint.mac

code:
	maxima -b ssa_code.mac

clean:
	rm -f *.pdf *.aux *.log *.rel
