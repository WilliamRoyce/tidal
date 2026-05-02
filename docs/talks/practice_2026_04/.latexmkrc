# Ensure the Cambridge house-style theme is on the search path for any
# latexmk invocation (make, LaTeX Workshop, CLI).
$ENV{'TEXINPUTS'} = './theme:' . ($ENV{'TEXINPUTS'} // '');
