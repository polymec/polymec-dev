all: poisson.pdf

poisson.pdf: *.tex
	pdflatex poisson
	bibtex poisson
	pdflatex poisson

clean:
	@rm *.idx *.aux *.log *.blg *.bbl *.pdf
