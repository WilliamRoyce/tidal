TEXINPUTS := ./theme:$(TEXINPUTS)
export TEXINPUTS

MAIN := poster

.PHONY: all clean watch figures

all: $(MAIN).pdf

# Gemini theme requires LuaLaTeX (fontspec).
$(MAIN).pdf: $(MAIN).tex sections/*.tex poster-macros.tex theme/*.sty figures/*.pdf
	latexmk -lualatex -interaction=nonstopmode -halt-on-error $(MAIN).tex

# Regenerate the bespoke poster figures (Cherry/Crest palette).
# Schematics + result corner (Python), then the magnetar graphical abstract (TikZ).
figures:
	cd scripts && python pipeline_schematic.py && python result_landscape.py
	cd figures && latexmk -pdf -interaction=nonstopmode magnetar_conversion.tex && latexmk -c magnetar_conversion.tex

watch:
	latexmk -lualatex -pvc -interaction=nonstopmode $(MAIN).tex

clean:
	latexmk -C
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb *.bbl *.blg *.fls *.fdb_latexmk
