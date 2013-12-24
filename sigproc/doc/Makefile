all: sig_freq_cw.pdf

sig_freq_cw.dvi: sig_freq_cw.eps

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
