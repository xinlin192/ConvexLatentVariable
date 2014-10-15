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

double compute_cost (vector<Instance*> data, vector< vector<double> >& global_means, vector<vector<int> > k, vector<int> z, vector<double> lambdas, Lookups* tables, dist_func df, int FIX_DIM) {
    double lambda_global = lambdas[0];
    double lambda_local = lambdas[1];
    int num_global_means = global_means.size();
    int num_local_means = 0;
    int D = tables->nDocs;
    for (int d = 0; d < D; d++) 
        num_local_means += k[d].size();
    double global_penalty = lambda_global * num_global_means;
    double local_penalty = lambda_local * num_local_means;

    int N = tables->nWords;
    vector< pair<int, int> > doc_lookup = *(tables->doc_lookup);
    vector<Instance*> temp_global_means (num_global_means, NULL);
    for (int p = 0; p < num_global_means; p ++) 
        temp_global_means[p] = vec2ins (global_means[p]);
    vector<int> global_asgn (N, 0);
    compute_assignment (global_asgn, k, z, tables);
    double loss = 0;
    for (int d = 0; d < D; d ++) 
        for (int i = doc_lookup[d].first; i < doc_lookup[d].second; i++) {
            double dist = df(data[i], temp_global_means[global_asgn[i]], FIX_DIM);
            loss += dist * dist;
        }
    loss *= 0.5;
    for (int p = 0; p < num_global_means; p ++) delete temp_global_means[p];

    double total = loss + global_penalty + local_penalty;
    cout << "loss: " << loss 
        << ", global: " << global_penalty 
        << ", local: " << local_penalty 
        << ", total: " << total << endl;
    return total;
}

double HDP_MEANS (vector<Instance*>& data, vector< vector<double> >& means, Lookups* tables, vector<double> lambdas, dist_func df, int FIX_DIM) {
    // STEP ZERO: validate input and initialization
    int N = tables->nWords;
    int D = tables->nDocs;
    vector< pair<int, int> > doc_lookup = *(tables->doc_lookup);
    double lambda_global = lambdas[0];
    double lambda_local = lambdas[1];

    vector< vector<double> > global_means (1, vector<double>(FIX_DIM, 0.0));
    vector< vector<int> > k (D, vector<int>(1,0));  // global association
    vector<int> z (N, 0); // local assignment
    vector<int> global_asgn (N, 0); // global assignment

    // STEP ONE: a. set initial *global* medoid as global mean
    compute_assignment (global_asgn, k, z, tables);
    compute_means (data, global_asgn, FIX_DIM, global_means);

    double last_cost = compute_cost (data, global_means, k, z, lambdas, tables, df, FIX_DIM);
    double new_cost = last_cost;
    while (true) {
        // 4. for each point x_ij,
        for (int j = 0; j < D; j ++) {
            for (int i = doc_lookup[j].first; i < doc_lookup[j].second; i++) {
                int num_global_means = global_means.size();
                vector<double> d_ij (num_global_means, 0.0);
                for (int p = 0; p < num_global_means; p ++) {
                    Instance* temp_ins = vec2ins(global_means[p]);
                    double euc_dist = df(data[i], temp_ins, FIX_DIM);
                    d_ij[p] = euc_dist * euc_dist;
                    delete temp_ins;
                }
                set<int> temp;
                for (int p = 0; p < num_global_means; p ++) temp.insert(p);
                int num_local_means = k[j].size();
                for (int q = 0; q < num_local_means; q ++) temp.erase(k[j][q]);
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
                    vector<double> new_g(FIX_DIM, 0.0);
                    for (int f = 0; f < data[i]->fea.size(); f++)
                        new_g[data[i]->fea[f].first-1] = data[i]->fea[f].second;
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
            int begin_i = doc_lookup[j].first;
            int end_i = doc_lookup[j].second;
            int doc_len = doc_lookup[j].second - doc_lookup[j].first;
            int num_local_means = k[j].size();
            // compute means of local clusters
            vector<int> local_asgn;
            vector< vector<double> > local_means (num_local_means, vector<double>(FIX_DIM, 0.0) );
            for (int i = begin_i; i < end_i; i ++) 
                local_asgn.push_back(z[i]);
            vector<Instance*> local_data (data.begin()+begin_i,data.begin()+end_i);
            compute_means (local_data, local_asgn, FIX_DIM, local_means);
            num_local_means = local_means.size();
            // compute distance of local clusters to each global cluster
            int num_global_means = global_means.size();
            vector<Instance*> temp_global_means (num_global_means, NULL);
            for (int p = 0; p < num_global_means; p ++) 
                 temp_global_means[p] = vec2ins (global_means[p]);
            vector< vector<double> > d_jcp (num_local_means, vector<double>(num_global_means, 0.0));
            vector<Instance*> temp_local_means (num_local_means, NULL);
            for (int c = 0; c < num_local_means; c ++) {
                 temp_local_means[c] = vec2ins (local_means[c]);
            }
            vector<double> sum_d_ijc (num_local_means, 0.0); 
            for (int i = doc_lookup[j].first; i < doc_lookup[j].second; i ++) {
                double local_dist = df (temp_local_means[z[i]], data[i], FIX_DIM);
                sum_d_ijc[z[i]] += local_dist * local_dist;
                for (int p = 0; p < num_global_means; p ++) {
                    double dist = df (temp_global_means[p], data[i], FIX_DIM);
                    d_jcp[z[i]][p] += dist * dist;
                }
            }
            for (int c = 0; c < num_local_means; c ++) {
                temp_local_means[c]->fea.clear();
                delete temp_local_means[c];
            }
            for (int p = 0; p < num_global_means; p ++) delete temp_global_means[p];
            for (int c = 0; c < num_local_means; c ++) {
                int num_global_means = global_means.size();
                int min_p = -1; double min_d_jcp = INF;
                for (int p = 0; p < num_global_means; p ++) 
                    if (d_jcp[c][p] < min_d_jcp) {
                        min_p = p;
                        min_d_jcp = d_jcp[c][p];
                    }
                if (min_d_jcp > lambda_global + sum_d_ijc[c]) {
                    global_means.push_back(local_means[c]); //  push mu_jc
                    k[j].push_back(num_global_means);
                } else {
                    k[j][c] = min_p;
                }
            }
        }
        // 6. for each global clusters,
        compute_assignment (global_asgn, k, z, tables);
        compute_means (data, global_asgn, FIX_DIM, global_means);

        // 7. convergence?
        new_cost = compute_cost (data, global_means, k, z, lambdas, tables, df, FIX_DIM);
        if (new_cost != last_cost) {
            last_cost = new_cost;
        } else break;
    }
    return last_cost;
}

// entry main function
int main (int argc, char ** argv) {
    if (argc < 5) {
        cerr << "Usage: " << endl;
        cerr << "\thdp_medoids [dataFile] [nRuns] [lambda_global] [lambda_local]" << endl;
        exit(-1);
    }

    // PARSE arguments
    char* dataFile = argv[1];
    int nRuns = atoi(argv[2]);
    vector<double> LAMBDAs (2, 0.0);
    LAMBDAs[0] = atof(argv[3]); // lambda_global
    LAMBDAs[1] = atof(argv[4]); // lambda_local

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
    double** dist_mat = mat_init (N, N);
    compute_dist_mat (data, dist_mat, N, FIX_DIM, df, true);
    vector< vector<double> > means;
    double obj = HDP_MEANS (data, means, &lookup_tables, LAMBDAs, df, FIX_DIM);
     
    /* Output objective */ 
   // output_objective (dist_mat, W, &lookup_tables, r, LAMBDAs);
    /* Output cluster centroids */
   // output_model (W, &lookup_tables);
    /* Output assignment */
   // output_assignment (W, &lookup_tables);

    /* reallocation */
    mat_free (W, N, N);
    mat_free (dist_mat, N, N);
}
