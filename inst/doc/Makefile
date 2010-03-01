all: overview nimbleexpr affysnp clean

overview: V1Overview.tex
	cp -p ../scripts/V1Overview.pdf .

nimbleexpr: V2NExpression.tex
	cp -p ../scripts/V2NExpression.pdf .

affysnp: V3AffySnpGenotype.tex
	cp -p ../scripts/V3AffySnpGenotype.pdf .

clean:
	$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	$(RM) -f Rplots.pdf V2NExpression-*
