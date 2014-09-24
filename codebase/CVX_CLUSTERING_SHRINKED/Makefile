cpps=$(shell ls *.cpp)

all: 
	g++ -O3 $(cpps) -o cvx_clustering 

push:
	git add cvx_clustering.cpp Makefile 
	git commit -m "upgrade"
	git push

clean:
	rm cvx_clustering opt*
