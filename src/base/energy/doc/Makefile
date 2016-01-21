all: discretization.pdf

discretization.pdf: discretization.tex formulas.tex
	pdflatex discretization.tex
	pdflatex discretization.tex
	pdflatex discretization.tex

formulas.tex: discretization.mac generic_equation.mac neumann_bc.mac typesetting.mac
	maxima -q -b typesetting.mac

clean:
	@rm -f formulas.tex *~ *.aux *.toc *.log *.out
