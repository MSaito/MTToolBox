.SUFFIXES: .tex .dvi .pdf .ps

doc: mttoolbox.pdf

mttoolbox.pdf:mttoolbox.tex xorshift-11.cpp xorshift-12.cpp \
	xorshift-21.cpp xorshift-22.cpp

.tex.dvi:
	platex $<

.dvi.ps:
	dvips $<

.ps.pdf:
	ps2pdf14 $<
