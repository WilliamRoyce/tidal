TEXINPUTS := ./theme:$(TEXINPUTS)
export TEXINPUTS

MAIN := poster

.PHONY: all clean watch figures

all: $(MAIN).pdf

# Gemini theme requires LuaLaTeX (fontspec).
$(MAIN).pdf: $(MAIN).tex sections/*.tex poster-macros.tex theme/*.sty figures/*.pdf
	latexmk -lualatex -interaction=nonstopmode -halt-on-error $(MAIN).tex

# Regenerate the bespoke poster figures (Cherry/Crest palette).
figures:
	cd scripts && python conversion_schematic.py && python pipeline_schematic.py && python result_landscape.py

watch:
	latexmk -lualatex -pvc -interaction=nonstopmode $(MAIN).tex

clean:
	latexmk -C
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb *.bbl *.blg *.fls *.fdb_latexmk
