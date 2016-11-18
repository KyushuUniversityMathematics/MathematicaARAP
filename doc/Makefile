all: MARAP.pdf Appendix.pdf
MARAP.pdf: MARAP.texi BasicFunctions.texi GlobalInterpolations.texi ConstraintFunctions.texi LocalInterpolations.texi
	texi2pdf MARAP.texi
Appendix.pdf: Appendix.tex
	platex Appendix.tex
	dvipdfmx Appendix.dvi
MARAP.html: MARAP.texi
	texi2html MARAP.texi
clean:
	rm -rf *.dvi *.aux *.cp *.fn *.fns *.ky *.log *.pg *.toc *.tp *.vr *~
veryclean: clean
	rm -rf MARAP.html MARAP.pdf Appendix.pdf
