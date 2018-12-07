all: ssafem_notes.pdf

%.pdf: %.tex
	pdflatex $^
	bibtex $(^:.tex=)
	pdflatex $^
	pdflatex $^

clean:
	rm -f *~ *.aux *.log *.rel *.out *.blg *.bbl *.synctex.gz *.toc
