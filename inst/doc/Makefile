all: overview nimbleexpr affysnp clean

overview: 01Overview.tex
	cp -p ../scripts/01Overview.pdf .

nimbleexpr: 02NimbleExpressionPreprocessing.tex
	cp -p ../scripts/02NimbleExpressionPreprocessing.pdf .

affysnp: 03AffySnpGenotype.tex
	cp -p ../scripts/03AffySnpGenotype.pdf .

clean:
	$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	$(RM) -f Rplots.pdf 02NimbleExpressionPreprocessing-*
