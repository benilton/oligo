all: overview nimbleexpr affysnp experimental exongene clean

overview: V1Overview.tex
	cp -p ../inst/scripts/V1Overview.pdf .

nimbleexpr: V2NExpression.tex
	cp -p ../inst/scripts/V2NExpression.pdf .

affysnp: V3AffySnpGenotype.tex
	cp -p ../inst/scripts/V3AffySnpGenotype.pdf .

experimental: V4Experimental.tex
	cp -p ../inst/scripts/V4Experimental.pdf .

exongene: V4Experimental.tex
	cp -p ../inst/scripts/V5ExonGene.pdf .

clean:
	$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	$(RM) -f Rplots.pdf V2NExpression-*
