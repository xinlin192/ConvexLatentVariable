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

/*
void compute_means (vector<Instance*>& means, vector< pair<int,int> >& words, vector<int>& assignment, int DIM) { 
    // compute number of existing clusters
    int N = words.size();
    int nClusters = means.size();  // for 0-index
    assert (N == assignment.size());
    // clear means
    for (int c = 0; c < nClusters; c++) delete means[c];
    // compute sums and counts
    vector< vector<double> > means_count (nClusters, vector<int>(DIM, 0.0));
    for (int i = 0; i < N; i++) {
        int voc_idx = words[i].first;
        int freq = words[i].second;
        means_count[assignment[i]][voc_idx] += freq;
    }
    // compute means
    for (int c = 0; c < nClusters; c ++) {
        Instance* tmp_cluster = new Instance(1000+c);
        for (int j = 0; j < DIM; j ++) 
            if (means_count[c][j] > 0)
                tmp_cluster->fea.push_back(make_pair(j+1, means_count[c][j]));
        means.push_back(tmp_cluster);
    }
}
*/

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
