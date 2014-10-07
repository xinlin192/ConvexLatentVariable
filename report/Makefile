all:
	pdflatex Examplar-based_MAD-Bayes.tex

push:
	git add Examplar-based_MAD-Bayes.tex Examplar-based_MAD-Bayes.pdf Makefile
	git commit -m "report tex file and pdf file upgrade"
	git push

clean: 
	rm -rf *.aux *.log *.bak *.out
