#include "HDP_MEANS.h"

/*
    computes distance of one word to a document 
 */
double word_doc_dist (pair<int, int> word_pair, Instance* doc_ins) {
    int voc_idx = word_pair.first;
    int count_w_d1 = word_pair.second;
    double dist;
    map<int, double> distribution ();
    // a. compute sum of word frequency
    int sumFreq = 0;
    int num_word_in_doc = doc_ins->fea.size();
    for (int f = 0; f < num_word_in_doc; f++) 
        sumFreq += doc_ins->fea[f].second;
    // b. compute distribution
    for (int f = 0; f < num_word_in_doc; f++) {
        int voc_index = doc_ins->fea[f].first - 1;
        double prob = 1.0 * doc_ins->fea[f].second / sumFreq;
        distribution.insert(pair<int, double> (voc_index, prob));
    }
    // c. compute weight of word within document
    map<int, double>::const_iterator iter;
    doc_ins->doc_distribution[voc_idx];
    count_w_di = freq;
    if (iter == distribution.end()) { // no count
        prob_w_d2 = 0.0;
        dist = INF;
    } else {
        prob_w_d2 = iter->second;
        dist = - count_w_d1 * log(prob_w_d2);
    }
    return dist;
}

/* 
    data: N data instance
    assignment: N by 1 vector indicating assignment of one atom
*/
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

void compute_assignment (vector<int>& assignment, vector< vector<int> >& k, vector<int>& z, Lookups* tables) {
    int N = z.size();
    int D = tables->nDocs; 
    assert (N == assignment.size());
    vector< pair<int,int> > doc_lookup = *(tables->doc_lookup);
    for (int j = 0; j < D; j ++)
        for (int i = doc_lookup[j].first; i < doc_lookup[j].second; i++) 
            assignment[i] = k[j][z[i]];
}

double HDP_MEANS (vector<Instance*>& data, Lookups* tables, vector<double> lambdas, dist_func df) {
    // STEP ZERO: validate input and initialization
    int N = tables->nWords;
    int D = tables->nDocs;
    int V = tables->nVocs;
    double lambda_global = lambdas[0];
    double lambda_local = lambdas[1];

    vector<Instance*> global_means (1, new Instance(0));
    vector< vector<Instance*> > local_means (D, vector<Instance*>());
    vector< vector<int> > k (D, vector<int>(1,0));  // global association
    vector<int> z (N, 0); // local assignment
    vector<int> assignment (N, 0); // global assignment

    // STEP ONE: a. set initial *global* medoid as global mean
    compute_assignment (assignment, k, z, tables);
    compute_means (global_means, word_lookup, vector<int>& assignment, V); 

    while (true) {
        // 4. for each point x_ij
        for (int j = 0; j < D; j ++) {
            for (int i = doc_lookup[j].first; i < doc_lookup[j].second; i++) {
                int num_global_means = global_means.size();
                vector<double> d_ij (num_global_means, 0.0);
                for (int p = 0; p < num_global_means; p ++) 
                    d_ij[p] = word_doc_dist(word_lookup[i], global_means[p]);
                set<int> temp;
                for (int p = 0; p < num_global_means; p ++) temp.insert(p);
                int num_local_means = k[j].size();
                for (int q = 0; q < num_local_means; q ++) temp.erase(k[q]);
                set<int>::iterator it; 
                for (it=temp.begin(); it!=temp.end();++it) d_ij[*it] += lambda_local;
                int min_p = -1; double min_dij = INF;
                for (int p = 0; p < num_global_means; p ++) 
                    if (d_ij[p] < min_dij) {
                        min_p = p;
                        min_dij = d_ij[p];
                    }
                if (min_dij > lambda_global + lambda_local) {
                    z[i] = num_local_means; 
                    k[j].push_back(num_global_means);
                    global_means.push_back();
                } else {

                }
            }
        }
    }
    /*
    vector< vector<double> > last_means = new_means; 
    double last_cost = INF, new_cost = INF; // compute last cost
    while (true) {
        // STEP TWO: compute dist square to new medoids d_ic
        for (int i = 0; i < N; i ++) {
            int nClusters = new_means.size();
            // cout << "nClusters: " << nClusters << endl;
            vector<double> dist_vec (nClusters, 0.0);
            for (int c = 0; c < nClusters; c ++) {
                Instance* mean_ins = new Instance (10000+c);
                for (int f = 0; f < D; f++) 
                    mean_ins->fea.push_back(make_pair(f+1, new_means[c][f]));
                dist_vec[c] = df(mean_ins, data[i], D);
                delete mean_ins;
            }
            int min_index = -1;
            double min_value = INF;
            for (int j = 0; j < nClusters; j++) 
                if (dist_vec[j] < min_value) {
                    min_index = j;
                    min_value = dist_vec[j];
                }
            // cout << "min_value: " << min_value << endl;
            if (min_value <= lambda) 
                assignment[i] = min_index;
            else {
                assignment[i] = nClusters;
                vector<double> tmp_mean (D, 0.0);
                int F = data[i]->fea.size();
                for (int f = 0; f < F; f ++) 
                    tmp_mean[data[i]->fea[f].first-1] = data[i]->fea[f].second;
                new_means.push_back(tmp_mean);
            } 
        }
        compute_means (data, assignment, D, new_means);
        // STEP THREE: compute cost
        int nClusters = new_means.size();
        vector<Instance*> means_ins (nClusters, NULL);
        for (int c = 0; c < nClusters; c++) {
            means_ins[c] = new Instance (5000+c);
            for (int j = 0; j < D; j ++) 
                means_ins[c]->fea.push_back(make_pair(j+1, new_means[c][j]));
        }
        new_cost = 0.0;
        for (int i = 0; i < N; i ++) 
            new_cost += 0.5 * df (means_ins[assignment[i]], data[i], D);
        new_cost += lambda * nClusters;
        for (int c = 0; c < nClusters; c++) delete means_ins[c];
        cout << "CLUSTERING COST: " << new_cost << endl;
        cout << "c: " << nClusters << endl;
        // STEP FOUR: convergence evaluation
        if (new_cost != last_cost) {
            last_cost = new_cost;
            last_means = new_means;
        } else break;
    }
    */
    return last_cost;
}

// entry main function
int main (int argc, char ** argv) {
    if (argc < 5) {
        cerr << "Usage: " << endl;
        cerr << "\thdp_medoids [voc_dataFile] [doc_dataFile] [lambda_global] [lambda_local]" << endl;
        exit(-1);
    }

    // PARSE arguments
    string voc_file (argv[1]);
    string doc_file (argv[2]);
    vector<double> LAMBDAs (2, 0.0);
    LAMBDAs[0] = atof(argv[3]); // lambda_document
    LAMBDAs[1] = atof(argv[4]); // lambda_topic

    // preprocess the input dataset
    vector<string> voc_list;
    voc_list_read (voc_file, &voc_list);
    int nVocs = voc_list.size();

    // init lookup_tables
    vector< pair<int,int> > doc_lookup;
    vector< pair<int,int> > word_lookup;
    vector< vector<int> > voc_lookup (nVocs, vector<int>());
    Lookups lookup_tables;
    lookup_tables.doc_lookup = &doc_lookup;
    lookup_tables.word_lookup = &word_lookup;
    lookup_tables.voc_lookup = &voc_lookup;

    document_list_read (doc_file, &lookup_tables);
    lookup_tables.nDocs = lookup_tables.doc_lookup->size();
    lookup_tables.nWords = lookup_tables.word_lookup->size();
    lookup_tables.nVocs = nVocs;

    int seed = time(NULL);
    srand (seed);
    cerr << "###########################################" << endl;
    cerr << "nVocs = " << lookup_tables.nVocs << endl; // # vocabularies
    cerr << "nDocs = " << lookup_tables.nDocs << endl; // # documents
    cerr << "nWords = " << lookup_tables.nWords << endl; // # words
    cerr << "lambda_global = " << LAMBDAs[0] << endl;
    cerr << "lambda_local = " << LAMBDAs[1] << endl;
    cerr << "TRIM_THRESHOLD = " << TRIM_THRESHOLD << endl;
    cerr << "seed = " << seed << endl;
    cerr << "###########################################" << endl;

    int N = lookup_tables.nWords;
    int D = lookup_tables.nDocs;
    double** dist_mat = mat_init (N, D);
    HDP_MEANS (data, N, D, lambda, df);

    /* Output objective, model and assignments */
    /*
    output_objective(clustering_objective (dist_mat, min_w, N));
    output_model (min_w, N);
    output_assignment (min_w, data, N);
    */
    /* Deallocation */

}
