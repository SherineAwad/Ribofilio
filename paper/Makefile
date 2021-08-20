all: ribofilio.pdf

ribofilio.pdf: ribofilio.tex references.bib
	pdflatex ribofilio
	bibtex ribofilio
	pdflatex ribofilio

clean:
	-rm -f ribofilio.aux ribofilio.bbl ribofilio.blg \
		ribofilio.log ribofilio.pdf
