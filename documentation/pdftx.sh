#!/bin/sh
# need to run dos2unix after editting
# run pdf latex with bibtex
pdflatex $1.tex &&
pdflatex $1.tex &&
cmd /c start SumatraPDF $1.pdf

rm *.log *.aux *.bbl *.blg *.ps *.out *.dvi