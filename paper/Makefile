.PHONY: all clean

all: main.pdf

main.pdf: main.tex references.bib sn-jnl.cls sn-vancouver.bst
	pdflatex -interaction=nonstopmode main.tex
	bibtex main
	pdflatex -interaction=nonstopmode main.tex
	pdflatex -interaction=nonstopmode main.tex

clean:
	rm -f main.aux main.bbl main.blg main.log main.out main.spl main.toc main.pdf
