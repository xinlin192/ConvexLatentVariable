cpps=$(shell ls *.cpp)

all: 
	g++ -fopenmp -O3 $(cpps) -o DP_convex

push:
	git add DP_convex.cpp DP_convex.h Makefile 
	git commit -m "upgrade"
	git push

clean:
	rm -rf DP_convex opt* dist_mat w_out plot*
