all: poisson.pdf solid.pdf

poisson.pdf: *.tex
	pdflatex poisson
	bibtex poisson
	pdflatex poisson
	bibtex poisson
	pdflatex poisson

solid.pdf: *.tex
	pdflatex solid
	bibtex solid
	pdflatex solid
	bibtex solid
	pdflatex solid

clean:
	@rm *.idx *.aux *.log *.blg *.bbl *.pdf
