#include <cassert>
#include <cmath>
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


