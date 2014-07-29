#include <iostream>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <map>

using namespace std;

typedef vector<int> FreqList;

int sample(FreqList& freq, int N, int alpha){
	
	int r = rand()% (N+alpha);
	int sum = 0;
	for(int i=0;i<freq.size();i++){
		sum += freq[i];
		if( r < sum )
			return i;
	}
	
	return freq.size();
}

void addToFreq(int x, FreqList& freq, int& N){
	
	if( x < freq.size() )
		freq[x]++;
	else{
		int ind = freq.size();
		while( ind < x ){
			freq.push_back(0);
			ind++;
		}
		freq.push_back(1);
	}
	N++;
}

void addToDoc(int x, map<int,int>& doc){
	
	map<int,int>::iterator it = doc.find(x);
	if( it != doc.end() )
		it->second++;
	else
		doc.insert(make_pair(x,1));
}

int main(int argc, char** argv){
		
	if( argc < 6 ){
		cerr << "./simdata_DP [num_doc] [doc_len] [alpha_topic] [alpha_topic_word] [alpha_word]" << endl;
		exit(0);
	}
	
	int num_doc = atoi(argv[1]);
	int num_word_p_doc = atoi(argv[2]);
	int alpha_topic = atoi(argv[3]);
	int alpha_topic_word = atoi(argv[4]);
	int alpha_word = atoi(argv[5]);
	
	vector<FreqList> topic_word_freq;
	vector<int> topic_word_sum;
	
	FreqList topic_freq;
	int topic_sum = 0;

	FreqList word_freq;
	int word_sum = 0;
	
	vector<map<int,int> > docs;
	int topic, word;
	for(int i=0;i<num_doc;i++){
		docs.push_back(map<int,int>());
		
		//Sample topic according to Hierachical DP
		topic = sample(topic_freq, topic_sum, alpha_topic);
		if( topic >= topic_freq.size() ){
			topic_word_freq.push_back(FreqList());
			topic_word_sum.push_back(0);
		}
		addToFreq(topic, topic_freq, topic_sum);
		
		for(int j=0;j<num_word_p_doc;j++){
			
			//Sample word according to DP|topic
			word = sample( topic_word_freq[topic], topic_word_sum[topic], alpha_topic_word);
			if( word >= topic_word_freq[topic].size() ){
				word = sample(word_freq, word_sum, alpha_word);
			}	
			addToFreq( word, word_freq, word_sum );
			addToFreq( word, topic_word_freq[topic], topic_word_sum[topic] );
			
			//assign the word to the doc
			addToDoc( word, docs[i] );
		}
	}
	
	//write result to file
	
	////docs
	ofstream fout("docs");
	for(int i=0;i<docs.size();i++){
	
		map<int,int>::iterator it;
		for(it=docs[i].begin();it!=docs[i].end();it++){
			fout << it->first << ":" << it->second << " ";
		}
		fout << endl;
	}
	fout.close();

	////topic freq
	fout.open("topic_freq");
	for(int i=0;i<topic_freq.size();i++){
		fout << i << " " << topic_freq[i] << endl;
	}
	fout.close();
	
	////word freq
	fout.open("word_freq");
	for(int i=0;i<word_freq.size();i++){
		fout << i << " " << word_freq[i] << endl;;
	}
	fout.close();
	
	////topic word freq
	fout.open("topic_word_freq");
	for(int i=0;i<topic_word_freq.size();i++){
		fout << i << " ";
		for(int j=0;j<topic_word_freq[i].size();j++){
			if( topic_word_freq[i][j] != 0)
				fout << j << ":" << topic_word_freq[i][j] << " ";
		}
		fout << endl;
	}
	fout.close();
}
