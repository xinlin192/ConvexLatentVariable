#include "../util.h"
#include <time.h>
#include <queue>

typedef struct {
    int nWords;
    int nVocs;
    int nDocs;
    vector< pair<int,int> >* doc_lookup; // document start,end word id pair
    vector< pair<int,int> >* word_lookup;  // word id to voc id
    vector< vector<int> >* voc_lookup;  // voc id to word id set
    vector< pair<int,int> >* word_in_doc;
} Lookups ;

class Int_Double_Pair_Dec
{
    public:
        bool operator() (pair<int, double> obj1, pair<int, double> obj2)
        {
            return obj1.second > obj2.second;
        }
};

void output_model (double** W, int R, int C) {
    ofstream model_out ("opt_model");
    vector<int> centroids;
    get_all_centroids (W, &centroids, R, C); // contains index of all centroids
    int nCentroids = centroids.size();
    model_out << "nCentroids: " << nCentroids << endl;
    for (int i = 0; i < nCentroids; i ++) {
        model_out << "centroids[" << i <<"]: " << centroids[i] << endl;
    }
    model_out.close();
}
void output_assignment (double ** W, Lookups * tables) {
    int N = tables->nWords;
    int D = tables->nDocs;
    vector< pair<int,int> > doc_lookup = *(tables->doc_lookup); ofstream asgn_out ("opt_assignment");
    vector< pair<int,int> > word_lookup = *(tables->word_lookup); 
    vector< vector<int> > voc_lookup = *(tables->voc_lookup);
    // get all centroids
    for (int d = 0; d < D; d++) {
        asgn_out << "d = " << d << endl;
        for (int i = doc_lookup[d].first; i < doc_lookup[d].second; i ++) {
            // output identification and its belonging
            asgn_out << "  id=" << i+1 << ", voc_index=" << word_lookup[i].first << ", "; 
            for (int j = 0; j < D; j ++) {
                if( fabs(W[i][j]) > 1e-2 ) {
                    asgn_out << j << "(" << W[i][j] << "), ";
                }
            }
            asgn_out << endl;
        }
    }
    asgn_out.close();
}

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
		if ( fin.eof() ) { 
            break;
        }
        vocList->push_back (line);
        // cout << line << endl;
	}
	// fin.close(); 
}

void voc_list_print (vector<string>* vocList) {
    int nVoc = vocList->size();
    for (int v = 0; v < nVoc; v ++) {
        cout << (*vocList)[v] << endl;
    }
}

/* word_lookup table restore the index in voc_list of vocabulary to which a word coresponds */
void document_list_read (string fname, Lookups* tables) {
    // document was zero-based 
    // word was zero-based
    // vocabulary was zero-based
    vector< pair<int,int> >* doc_lookup = tables->doc_lookup;
    vector< pair<int,int> >* word_lookup = tables->word_lookup; 
    vector< vector<int> >* voc_lookup = tables->voc_lookup;
    vector< pair<int,int> >* word_in_doc = tables->word_in_doc; 

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
            word_in_doc->push_back(make_pair(d, freq));
            ++ w;
        }
        doc_index_end = w;
        doc_lookup->push_back(make_pair(doc_index_begin, doc_index_end));
        ++ d;
	}
	fin.close(); 
    int nWords = w;
    for (int i = 0; i < nWords; i ++)  // 1-based index
        (*voc_lookup)[ (*word_lookup)[i].first ].push_back(i);
}
