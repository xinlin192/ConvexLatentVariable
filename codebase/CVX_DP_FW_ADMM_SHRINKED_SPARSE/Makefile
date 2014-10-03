cpps=$(shell ls *.cpp)

all: 
	g++ -O3 cvx_clustering.cpp -o cvx_clustering 

push:
	git add cvx_clustering.cpp Makefile  cvx_clustering.h
	git commit -m "upgrade"
	git push

test:
	g++ -O3 gs_test.cpp -o gs_test


clean:
	rm cvx_clustering opt* dist_mat w_out plot*
