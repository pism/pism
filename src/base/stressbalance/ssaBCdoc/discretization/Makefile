all: pdf

pdf: formulas.tex
	pdflatex ssa-discretization.tex
	pdflatex ssa-discretization.tex

formulas.tex: ssa_prettyprint.mac ssa_coeffs.mac ssa.mac
	maxima -b ssa_prettyprint.mac

code: ssa_code.mac ssa_coeffs.mac ssa.mac
	maxima -b ssa_code.mac

clean:
	rm -f *.pdf *.aux *.log *.rel
