cpp=$(shell ls *.cpp) 
headers=$(shell ls *.h)
executable=$(shell echo $(cpp) | cut -d '.' -f 1)
README=$(shell ls README*)

all: 
	g++ -O3 $(cpp) -o $(executable)

push:
	git add $(cpp) $(headers) $(README) Makefile
	git commit 
	git push

clean:
	rm -rf $(executable) opt*
