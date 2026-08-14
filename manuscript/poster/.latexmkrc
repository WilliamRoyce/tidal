# Ensure the Gemini + Cambridge theme is on the search path for any latexmk
# invocation (make, LaTeX Workshop, CLI). Gemini needs LuaLaTeX.
$ENV{'TEXINPUTS'} = './theme:' . ($ENV{'TEXINPUTS'} // '');
$pdf_mode = 4;   # 4 = lualatex
