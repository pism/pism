all: stressbalance-out.png

stressbalance-out.pdf: stressbalance-out.dot
	dot -o stressbalance-out.pdf -Tpdf stressbalance-out.dot 
	pdfcrop stressbalance-out.pdf stressbalance-out.pdf

stressbalance-out.png: stressbalance-out.dot
	dot -o stressbalance-out.png -Tpng stressbalance-out.dot 

