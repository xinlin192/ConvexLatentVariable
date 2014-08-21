cpp=$(shell ls *.cpp) 
headers=$(shell ls *.h)
executable=$(shell ls *.cpp) | tr -d '.cpp'
README=$(shell ls README*)

all: 
	g++ -O3 $(cpp) -o $(executable)

push:
	git add $(cpp) $(headers) $(README) Makefile
	git commit
	git push

clean:
	rm $(executable) 
