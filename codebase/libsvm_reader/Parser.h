#include"Instance.h"
#include<stdio.h>
#include<stdlib.h>
#include<fstream>
#include<vector>
#include<map>
#include<iostream>
#include<cstring>

using namespace std;

class Parser{
	
	public:
	static  vector< Instance* >*  parseSVM(char* fileName, int& numFea){
		
		ifstream fin(fileName);
		if( fin.fail() ){
			cerr << "cannot find file." << endl;
			exit(0);
		}
		
		vector< Instance* >* data = new vector< Instance* >();
		
		char line[10000];
		vector<string> tokens;
		numFea=0;
		int id=0;
		while( !fin.eof() ){
			
			Instance* ins = new Instance(id++);
			
			fin.getline(line,10000);
			string str = string(line);
			tokens = split(str," ");
			//if( str.length() < 2 )
			//	continue;

			//label
			ins->label = atoi(tokens[0].c_str());

			//feature
			for(int i=1;i<tokens.size();i++){
				
				vector<string> pv = split(tokens[i],":");
				pair<int,double> pair;
				pair.first = atoi(pv[0].c_str());
				pair.second = atof(pv[1].c_str());
				ins->feature.push_back(pair);
			}
			//cerr << "fea="<< ins->feature.back().second << endl;
			//cerr << data->size() << ", " << ins->feature.size() <<  endl;
			if( ins->feature.size()>0 && ins->feature.back().first > numFea )
				numFea = ins->feature.back().first;
			
			data->push_back( ins );
		}
		
		data->pop_back();
		return data;
	}

	static vector<string> split(string str, string pattern){

		vector<string> str_split;
		size_t i=0;
		size_t index=0;
		while( index != string::npos ){

			index = str.find(pattern,i);
			str_split.push_back(str.substr(i,index-i));

			i = index+1;
		}
		
		if( str_split.back()=="" )
			str_split.pop_back();
		
		return str_split;
	}

};
