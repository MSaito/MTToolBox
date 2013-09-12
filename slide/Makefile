.SUFFIXES: .tex .dvi .pdf .ps

.tex.dvi:
	platex-utf8 $<

.dvi.ps:
	dvips $<

.ps.pdf:
	ps2pdf14 $<

clean:
	rm *~ *.log