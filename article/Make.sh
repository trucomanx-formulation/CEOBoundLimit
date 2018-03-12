#!/bin/bash

pdflatex -interaction=nonstopmode ceoboundlimit.tex
pdflatex -interaction=nonstopmode ceoboundlimit.tex

./Clean.sh
