.SUFFIXES: .tex .dvi .pdf .ps

.tex.dvi:
	platex-euc $<

.dvi.ps:
	dvips $<

.ps.pdf:
	ps2pdf14 $<

clean:
	rm *~ *.log