#include<iostream>
#include<fstream>
using namespace std;

int main(){
	
	ifstream fin("input");
	int tmp;
	fin >> tmp ;
	fin.get();
	
	char c;
	int i=0;
	while(!fin.eof()){

		c = fin.get();
		cerr << i << ":" << (int)c << endl;
		i++;
	}
	fin.close();

	cerr << (int)('\n') << endl;
	return 0;
}
