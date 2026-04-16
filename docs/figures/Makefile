FIGURES := $(patsubst %.tex,%.pdf,$(wildcard *.tex))
# tidal-tikz-styles.sty is a package file, not a figure
FIGURES := $(filter-out tidal-tikz-styles.pdf,$(FIGURES))

.PHONY: all clean

all: $(FIGURES)

%.pdf: %.tex tidal-tikz-styles.sty
	latexmk -pdf -interaction=nonstopmode -halt-on-error $<

clean:
	latexmk -C $(patsubst %.pdf,%.tex,$(FIGURES))
	rm -f *.aux *.log *.out *.fls *.fdb_latexmk *.synctex.gz
