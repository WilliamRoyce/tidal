TEXINPUTS := .:../figures:$(TEXINPUTS)
export TEXINPUTS

MAIN := main
SECTIONS := $(filter-out $(MAIN).tex preamble.tex, $(wildcard *.tex))

.PHONY: all clean watch

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex preamble.tex $(SECTIONS) references.bib
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex

watch:
	latexmk -pdf -pvc -interaction=nonstopmode $(MAIN).tex

clean:
	latexmk -C
	rm -f *.aux *.log *.out *.toc *.bbl *.blg *.fls *.fdb_latexmk
