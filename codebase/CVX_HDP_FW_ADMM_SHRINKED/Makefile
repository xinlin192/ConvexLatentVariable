cpps=$(shell ls *.cpp)

all: 
	g++ -O3 $(cpps) -o cvx_hdp_medoids 

push:
	git add cvx_hdp_medoids.cpp Makefile cvx_hdp_medoids.h
	git commit -m "upgrade"
	git push

clean:
	rm -rf cvx_hdp_medoids opt* dist_mat w_out plot* word2doc 
