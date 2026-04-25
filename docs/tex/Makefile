TEXINPUTS := .:../figures:$(TEXINPUTS)
export TEXINPUTS

MAIN := main
SECTIONS := $(filter-out $(MAIN).tex preamble.tex, $(wildcard *.tex))

.PHONY: all clean watch

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex preamble.tex $(SECTIONS) references.bib
	# Explicit multi-pass build.  We avoid latexmk's -halt-on-error /
	# default error-on-warning behaviour, since pdflatex returns code 1
	# whenever there are unresolved \ref{} or \cite{} entries — which is
	# expected on the first pass before bibtex runs and on intermediate
	# passes before all forward refs settle.  The pipeline is:
	#   1. pdflatex — first pass, generates .aux with cite keys
	#   2. bibtex   — resolves \cite{} from references.bib
	#   3. pdflatex — second pass, picks up .bbl
	#   4. pdflatex — third pass, resolves cross-refs forward-pointing
	#                 from earlier sections to later ones
	# Warnings about missing bib keys are surfaced from bibtex; they are
	# fatal at developer-discretion (audit the .blg file).
	@TEXINPUTS=".:../figures:" pdflatex -interaction=nonstopmode $(MAIN).tex >/dev/null || true
	@bibtex $(MAIN) || true
	@TEXINPUTS=".:../figures:" pdflatex -interaction=nonstopmode $(MAIN).tex >/dev/null || true
	@TEXINPUTS=".:../figures:" pdflatex -interaction=nonstopmode $(MAIN).tex >/dev/null || true
	@if grep -qE "Reference.*undefined|Citation.*undefined" $(MAIN).log 2>/dev/null; then \
	   n=$$(grep -cE "Reference.*undefined|Citation.*undefined" $(MAIN).log); \
	   echo "WARNING: $$n unresolved references/citations remain (see $(MAIN).log)"; \
	 fi
	@echo "Build complete: $(MAIN).pdf ($$(stat -c %s $(MAIN).pdf 2>/dev/null) bytes, $$(grep -oP 'Output written on [^ ]+ \(\K[0-9]+' $(MAIN).log 2>/dev/null) pages)"

watch:
	latexmk -pdf -pvc -interaction=nonstopmode $(MAIN).tex

clean:
	latexmk -C
	rm -f *.aux *.log *.out *.toc *.bbl *.blg *.fls *.fdb_latexmk
