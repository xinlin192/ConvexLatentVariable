cpp=$(shell ls *.cpp) 
headers=$(shell ls *.h)
executable=$(shell ls *.cpp) | tr -d '.cpp'

all: 
	g++ -O3 $(cpp) -o $(executable)

push:
	git add $(cpp) $(headers)
	git commit
	git push

clean:
	rm $(executable) 
