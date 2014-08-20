# APPLICATION MAKEFILE
# Jimmy Lin <JimmyLin@utexas.edu> 

#######################################################################
# DO NOT EDIT THIS MAKEFILE UNLESS YOU KNOW WHAT YOU ARE DOING. 
#######################################################################	

all: 
	g++ -O3 single.cpp -o single

HDP:
	g++ -O3 HDP.cpp -o HDP

test: testEsmat testHDP

testEsmat: testEsmat.cpp
	g++ -O3 testEsmat.cpp -o testEsmat
	./testEsmat

testHDP:
	g++ -O3 testHDP.cpp -o testHDP
	./testHDP 

clean:
	rm -rf single HDP testEsmat testHDP
