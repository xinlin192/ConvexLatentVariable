#include <cassert>
#include <cmath>
#include <omp.h>
#define INTEGER_MAX 3000000
#include "../util.h"

using namespace std;

typedef struct {
    int nWords;
    int nVocs;
    int nDocs;
    vector< pair<int,int> >* doc_lookup;
    vector< pair<int,int> >* word_lookup;
    vector< vector<int> >* voc_lookup;
} Lookups ;

void get_doc_lookup (vector< Instance* > & data, vector<pair<int,int> >& doc_lookup) {
    // validate the input
    doc_lookup.clear();
    int N = data.size();
    for (int i = 1; i < N; i ++) 
        assert(data[i-1]->label <= data[i]->label);
    // compute list of document info
    int doc_begin = 0;
    int doc_end = 0;
    int last_doc = data[0]->label;
    for (int i = 0; i < N ; i++) {
        int curr_doc = data[i]->label;
        if (curr_doc != last_doc) {
            doc_end = i;
            //cerr << "(" << doc_begin << ", " << doc_end << ")" <<endl;
            doc_lookup.push_back(make_pair(doc_begin, doc_end));
            doc_begin = i;
            last_doc = curr_doc;
        } 
    }
    // cerr << "(" << doc_begin << ", " << doc_end << ")" <<endl;
    doc_lookup.push_back(make_pair(doc_begin, N));
}
