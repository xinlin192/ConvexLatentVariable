# PROPOSAL FOR CONVEX LATENT VARIABLE PROJECT MAKEFILE
# Jimmy Lin <jimmylin@utexas.edu> 

#######################################################################
# DO NOT EDIT THIS MAKEFILE UNLESS YOU KNOW WHAT YOU ARE DOING. 
#######################################################################	

fname=nips2014
tex_file=$(fname).tex
pdf_file=$(fname).pdf

all:
	pdflatex $(tex_file)

push:
	git add $(tex_file) $(pdf_file)
	git commit -m "new update"
	git push

clean:
	rm -f $(fname).aux $(fname).log $(fname).snm $(fname).out texput.log
