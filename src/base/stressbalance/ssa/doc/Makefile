all: new_notes.pdf

new_notes.pdf: new_notes.tex
	pdflatex new_notes.tex
	bibtex new_notes
	pdflatex new_notes.tex

clean:
	rm -f *~ *.aux *.log *.rel *.out *.blg *.bbl *.synctex.gz
