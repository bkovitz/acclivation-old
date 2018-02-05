LAT = latex -shell-escape
#TEX = pdflatex --shell-escape
TEXINPUTS=.:./sty:
TEX = TEXINPUTS=.:./sty: latexmk -pdf -xelatex
BIB = bibtex8
TEXFILES = $(wildcard *.tex)
PDFFILES = $(TEXFILES:.tex=.pdf)
BIBFILES = $(wildcard *.bib)
DOT = dot -Tpdf

all: acclivation.pdf data

data: seed1-vfit seed1-pop

seed1-vfit seed1-pop: run.clj acclivation/core.clj
	./run.clj :seed 1 :step 0.005

$(PDFFILES): $(BIBFILES)
%.pdf: %.tex
	$(TEX) $<

%.pdf: %.dot
	$(DOT) < $< >$@

acclivation.pdf: papers.bib

%.pdf: %.data mkplot.gpi mkplots
	./mkplots

run1:
	lein run -m farg.acclivation/run-and-save :generations 20 :n-epochs 10 :seed 1

tags:
	ctags -R src/ checkouts/

clean:
	rm -f *.aux *.bcf *.log *.blg *.bbl *.toc *.dvi *.fls *.fdb_latexmk \
	   *.run.xml acclivation.pdf

.PHONY: tags clean data
