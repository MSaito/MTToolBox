.SUFFIXES: .tex .dvi .pdf .ps

mttoolbox.pdf:mttoolbox.tex xorshift-11.cpp xorshift-12.cpp \
	xorshift-21.cpp xorshift-22.cpp

.tex.dvi:
	platex-utf8 $<

.dvi.ps:
	dvips $<

.ps.pdf:
	ps2pdf14 $<

clean:
	rm *~ *.log