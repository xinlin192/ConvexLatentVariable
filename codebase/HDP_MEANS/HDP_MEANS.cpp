#include "HDP_MEANS.h"

Instance* vec2ins (vector<double> vec) {
    Instance* ins = new Instance (1);
    int DIM = vec.size();
    for (int f = 0; f < DIM; f ++) 
        ins->fea.push_back(make_pair(f+1, vec[f]));
    return ins;
}
/* 
    data: N data instance
    assignment: N by 1 vector indicating assignment of one atom
*/
void compute_means (vector<Instance*>& data, vector<int>& assignment, int FIX_DIM, vector< vector<double> >& means) {
    int N = data.size();
    assert (assignment.size() == N);
    // compute number of existing clusters
    int nClusters = means.size();  // for 0-index
    // clear means
    for (int j = 0; j < nClusters; j++) 
        std::fill(means[j].begin(), means[j].end(), 0.0);
    // compute sums and counts
    vector<int> means_count (nClusters, 0);
    for (int i = 0; i < N; i++) {
        int belongsto = assignment[i];
        for (int f = 0; f < data[i]->fea.size(); f++) {
            int j = data[i]->fea[f].first - 1;
            means[belongsto][j] += data[i]->fea[f].second;
        }
        ++ means_count[belongsto];
    }
    // compute means
    for (int c = 0; c < nClusters; c ++) 
        for (int j = 0; j < FIX_DIM; j++) 
            means[c][j] = 1.0 * means[c][j] / means_count[c];
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

double HDP_MEANS (vector<Instance*>& data, vector<Instance*>& means, Lookups* tables, vector<double> lambdas, dist_func df, int FIX_DIM) {
    // STEP ZERO: validate input and initialization
    int N = tables->nWords;
    int D = tables->nDocs;
    vector< pair<int, int> > doc_lookup = *(tables->);
    double lambda_global = lambdas[0];
    double lambda_local = lambdas[1];

    vector< vector<double> > global_means (1, new Instance(0));
    vector< vector<int> > k (D, vector<int>(1,0));  // global association
    vector<int> z (N, 0); // local assignment
    vector<int> global_asgn (N, 0); // global assignment

    // STEP ONE: a. set initial *global* medoid as global mean
    compute_assignment (global_asgn, k, z, tables);
    compute_means (data, global_asgn, FIX_DIM, global_means);

    while (true) {
        // 4. for each point x_ij,
        for (int j = 0; j < D; j ++) {
            for (int i = doc_lookup[j].first; i < doc_lookup[j].second; i++) {
                int num_global_means = global_means.size();
                vector<double> d_ij (num_global_means, 0.0);
                for (int p = 0; p < num_global_means; p ++) {
                    Instance* temp_ins = vec2ins(global_means[p]);
                    double euc_dist = df(data[i], temp_ins);
                    d_ij[p] = euc_dist * euc_dist;
                    delete temp_ins;
                }
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
                    Instance* new_g = new Instance* (10000+num_global_means);
                    new_g->fea = data[i]->fea;
                    global_means.push_back(new_g);
                } else {
                    bool c_exist = false;
                    for (int c = 0; c < num_local_means; c ++) 
                        if (k[j][c] == min_p) {
                            z[i] = c;
                            c_exist = true;
                        }
                    if (!c_exist) {
                        z[i] = num_local_means;
                        k[j].push_back(min_p);
                    }
                }
            }
        }
        // 5. for all local clusters,
        for (int j = 0; j < D; j ++) {
            int doc_len = doc_lookup[j].second - doc_lookup[j].first;
            int num_local_means = k[j].size();
            vector<int> local_asgn;
            vector<Instance*> local_data (data[doc_lookup[j].first], data[doc_lookup[j].second]);
            vector< vector<double> > local_means ( vector<Instance*>());
            for (i = doc_lookup[j].first; i < doc_lookup[j].second; i ++) 
                local_asgn.push_back(z[i]);
            compute_means (local_data, local_asgn, FIX_DIM, local_means);
            num_local_means = local_means.size();
            for (int c = 0; c < num_local_means; c ++) {
                int num_global_means = global_means.size();
                vector<double> d_jcp (num_global_means, 0);
                for (int p = 0; p < num_global_means; p ++) {


                }
                int min_p = -1; double min_d_ijp = INF;
                for (int p = 0; p < num_global_means; p ++) {

                }
                if (min_d_ijp > lambda_global + ) {
                    global_means.push_back(); //  push mu_jc
                    k[j].push_back(num_global_means);
                } else {
                    k[j][c] = min_p;
                }
            }
        }
        // 6. for each global clusters,
        compute_assignment (global_asgn, k, z, tables);
        compute_means (data, global_asgn, FIX_DIM, global_means);
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
        cerr << "\thdp_medoids [word_dataFile] [lambda_global] [nRuns] [lambda_local]" << endl;
        exit(-1);
    }

    // PARSE arguments
    char* dataFile = argv[1];
    vector<double> LAMBDAs (2, 0.0);
    LAMBDAs[0] = atof(argv[2]); // lambda_global
    LAMBDAs[1] = atof(argv[3]); // lambda_local
    int FW_MAX_ITER = atoi(argv[4]);
    int ADMM_MAX_ITER = atoi(argv[5]);
    int SS_PERIOD = 100;

    // read in data
    int FIX_DIM;
    Parser parser;
    vector<Instance*>* pdata;
    vector<Instance*> data;
    pdata = parser.parseSVM(dataFile, FIX_DIM);
    data = *pdata;

    // init lookup_tables
    vector< pair<int,int> > doc_lookup;
    get_doc_lookup (data, doc_lookup);
    Lookups lookup_tables;
    lookup_tables.doc_lookup = &doc_lookup;
    lookup_tables.nWords = data.size();
    lookup_tables.nDocs = lookup_tables.doc_lookup->size();
    int seed = time(NULL);
    srand (seed);
    cerr << "###########################################" << endl;
    cerr << "nDocs = " << lookup_tables.nDocs << endl; // # documents
    cerr << "nWords = " << lookup_tables.nWords << endl; // # words
    cerr << "lambda_global = " << LAMBDAs[0] << endl;
    cerr << "lambda_local = " << LAMBDAs[1] << endl;
    cerr << "TRIM_THRESHOLD = " << TRIM_THRESHOLD << endl;
    cerr << "seed = " << seed << endl;
    cerr << "###########################################" << endl;

    // Run sparse convex clustering
    int N = lookup_tables.nWords;
    int D = lookup_tables.nDocs;
    double** W = mat_init (N, N);
    // dist_mat computation and output
    dist_func df = L2norm;
    vector< Instance* > means;
    double obj = HDP_MEANS (data, &means, &lookup_tables, lambdas, df, FIX_DIM);

    /* Output objective */ 
    output_objective (dist_mat, W, &lookup_tables, r, LAMBDAs);
    /* Output cluster centroids */
    output_model (W, &lookup_tables);
    /* Output assignment */
    output_assignment (W, &lookup_tables);

    /* reallocation */
    mat_free (W, N, N);
    mat_free (dist_mat, N, N);
}
