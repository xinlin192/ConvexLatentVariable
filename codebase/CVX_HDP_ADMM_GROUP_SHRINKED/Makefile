cpps=$(shell ls *.cpp)

all: 
	g++ -O3 $(cpps) -o cvx_hdp_shrinked

push:
	git add cvx_hdp_shrinked.cpp cvx_hdp_shrinked.h Makefile 
	git commit -m "upgrade"
	git push

clean:
	rm cvx_hdp_shrinked opt* dist_mat w_out plot*
