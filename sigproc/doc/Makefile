all: sig_freq.pdf

sig_freq.dvi: sig_freq.eps

%.pdf: %.ps
	ps2pdf $<
%.ps: %.dvi
	dvips $<
%.dvi: %.tex
	latex $<
%.eps: %.png
	convert $< $@

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.dvi *.eps
