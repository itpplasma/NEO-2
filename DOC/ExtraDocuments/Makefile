all: 2018_Documentation.pdf neo-2-ql-demo_runs.pdf

%.pdf: %.ps
	ps2pdf $<

%.ps: %.dvi
	dvips $<

%.dvi: %.tex
	latex $*
	latex $*

clean:
	rm -f *.dvi *.ps *.aux *.bbl *.bcf *.blg *.idx *.ilg *.ind *.lof *.log *.nav *.out *.snm *.toc *.run.xml

help:
	@echo "Default makefile for creating a pdf document via latex."
	@echo "Without arguments all pdf's are created."
