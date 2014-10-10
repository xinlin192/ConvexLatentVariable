using namespace std;
#include<cassert>
#include<string>
#include<fstream>

#ifndef UTIL_H
#include "../util.h"
#endif

/*
#ifndef EXSPARSE_MAT_H
#include "../exSparseMat.h"
#endif
*/

typedef struct {
    int nWords;
    int nVocs;
    int nDocs;
    vector< pair<int,int> >* doc_lookup;
    vector< pair<int,int> >* word_lookup;
    vector< vector<int> >* voc_lookup;
} Lookups ;

/* Note that operations within this function do not destroy original input */
void split (string input, vector<string>* elements, string delimiter) {
    string str (input);
    int index = str.find_first_of(delimiter);
    while (!(index < 0)) {
        elements->push_back(str.substr(0, index));
        str = str.substr(index+1);
        index = str.find_first_of(delimiter);
    }
    if (str.size() > 0) 
        elements->push_back(str) ;  
}

/* TODO: consider the vocabulary is 0-based or 1-based */
void voc_list_read (string fname, vector<string>* vocList) {
   	ifstream fin(fname.c_str());

    string line;
	while (!fin.eof()) {
        getline(fin, line);
		if ( fin.eof() ) break;
        vocList->push_back (line);
	}
	fin.close(); 
}

void voc_list_print (vector<string>* vocList) {
    int nVoc = vocList->size();
    for (int v = 0; v < nVoc; v ++) {
        cout << (*vocList)[v] << endl;
    }
}

/* word_lookup table restore the index in voc_list of vocabulary to which a word coresponds */
void document_list_read (string fname, Lookups* tables) {
    vector< pair<int,int> >* doc_lookup = tables->doc_lookup;
    vector< pair<int,int> >* word_lookup = tables->word_lookup; 
    vector< vector<int> >* voc_lookup = tables->voc_lookup;

   	ifstream fin(fname.c_str());
	string line = "";
    int doc_index_begin = 0;
    int doc_index_end = 0;
    int d = 0, w = 0;
	while (!fin.eof()) {
        getline(fin, line);
		if ( fin.eof() ) break;
        // process new document
        doc_index_begin = w;
        // split the string to several field (delimiter: whitespace) 
        vector<string> fields;
        split (line, &fields, " ");
        for (int f = 0; f < fields.size(); f ++) {
            vector<string> voc_freq_pair;
            split (fields[f], &voc_freq_pair, ":");
            // for each word, split voc_index and frequency by ":"
            int voc_index = atoi(voc_freq_pair[0].c_str());
            int freq = atoi(voc_freq_pair[1].c_str());
            // push to word_lookup table
            word_lookup->push_back(make_pair(voc_index, freq));
            ++ w;
        }
        doc_index_end = w;
        doc_lookup->push_back(make_pair(doc_index_begin, doc_index_end));
        ++ d;
	}
	fin.close(); 
    int nWords = w;
    for (int i = 0; i < nWords; i ++) {
        (*voc_lookup)[ (*word_lookup)[i].first-1].push_back(i);
    }
} 
