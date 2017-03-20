PAPER_NAME = divpoly_pit

all:
	pdflatex $(PAPER_NAME).tex
	bibtex   $(PAPER_NAME)
	pdflatex $(PAPER_NAME).tex
	pdflatex $(PAPER_NAME).tex


