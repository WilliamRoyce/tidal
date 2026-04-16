TEXINPUTS := ./theme:$(TEXINPUTS)
export TEXINPUTS

MAIN := talk

.PHONY: all clean watch

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex sections/*.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex

watch:
	latexmk -pdf -pvc -interaction=nonstopmode $(MAIN).tex

clean:
	latexmk -C
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb *.bbl *.blg *.fls *.fdb_latexmk
