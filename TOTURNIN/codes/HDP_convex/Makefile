cpps=$(shell ls *.cpp)

all: 
	g++ -fopenmp -O3 $(cpps) -o HDP_convex

push:
	git add HDP_convex.cpp Makefile HDP_convex.h
	git commit -m "upgrade"
	git push

clean:
	rm -rf HDP_convex opt* dist_mat w_out plot* word2doc 
