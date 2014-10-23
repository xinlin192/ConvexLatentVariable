#include "HDP_MEANS.h"
ofstream objmin_trace("../obj_vs_time_dp/wholesale/HDP-MEAN");
double objmin = 1e300;
double start_time;

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
        std::fill(means[j].begin(), means[j].end(), -INF);
    vector< vector<double> > means_sum (nClusters, vector<double> (FIX_DIM,0.0));
    // compute sums and counts
    vector<int> means_count (nClusters, 0);
    for (int i = 0; i < N; i++) {
        int belongsto = assignment[i];
        for (int f = 0; f < data[i]->fea.size(); f++) {
            int j = data[i]->fea[f].first - 1;
            means_sum[belongsto][j] += data[i]->fea[f].second;
        }
        ++ means_count[belongsto];
    }
    // compute means
    for (int c = 0; c < nClusters; c ++) {
        // assert (means_count[c] != 0);
        if (means_count[c] > 0)
            for (int j = 0; j < FIX_DIM; j++) 
                means[c][j] = 1.0 * means_sum[c][j] / means_count[c];
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

int get_num_global_means (vector<vector<int> > k) {
    int nDocs = k.size();
    set<int> gmeans_index;
    for (int d = 0; d < nDocs; d ++) {
        int num_local_means = k[d].size();
        for (int c = 0; c < num_local_means; c++)
            gmeans_index.insert(k[d][c]);
    }
    return gmeans_index.size();
}
int get_num_local_means (vector<int> z, Lookups* tables) {
    int nDocs = tables->nDocs;
    vector< pair<int, int> > doc_lookup = *(tables->doc_lookup);
    vector< set<int> > lmeans_sets (nDocs, set<int> ());
    for (int j = 0; j < nDocs; j ++) 
        for (int i = doc_lookup[j].first; i < doc_lookup[j].second; i++) 
            lmeans_sets[j].insert(z[i]);
    int sum_locals = 0;
    for (int j = 0; j < nDocs; j ++)
        sum_locals += lmeans_sets[j].size();
    return sum_locals;
}
double compute_cost (vector<Instance*> data, vector< vector<double> >& global_means, vector<vector<int> > k, vector<int> z, vector<double> lambdas, Lookups* tables, dist_func df, int FIX_DIM) {
    double lambda_global = lambdas[0];
    double lambda_local = lambdas[1];
    int num_global_means = global_means.size();
    int num_local_means = 0;
    int D = tables->nDocs;
    for (int d = 0; d < D; d++) 
        num_local_means += k[d].size();
    double global_penalty = lambda_global * get_num_global_means(k);
    double local_penalty = lambda_local * get_num_local_means(z, tables);

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
    for (int p = 0; p < num_global_means; p ++) 
        delete temp_global_means[p];

    double total = loss + global_penalty + local_penalty;
    cerr << "loss: " << loss 
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
                    // cout << "global and local increment" << endl;
                } else {
                    bool c_exist = false;
                    for (int c = 0; c < num_local_means; c ++) 
                        if (k[j][c] == min_p) {
                            z[i] = c;
                            c_exist = true;
                            break;
                        }
                    if (!c_exist) {
                        z[i] = num_local_means;
                        k[j].push_back(min_p);
                       // cout << "local increment" << endl;
                    }
                }
            }
        }
        /*
        cout << "half..........." << endl;
        cout << "#global created: " << global_means.size() 
            << ", #global used: " << get_num_global_means(k);
            */
        new_cost = compute_cost (data, global_means, k, z, lambdas, tables, df, FIX_DIM);
        // 5. for all local clusters,
        for (int j = 0; j < D; j ++) {
            int begin_i = doc_lookup[j].first;
            int end_i = doc_lookup[j].second;
            int doc_len = doc_lookup[j].second - doc_lookup[j].first;
            int num_local_means = k[j].size();

            // all local clusters are distinct to each other
            /*
            set<int> temp;
            for (int y = 0; y < num_local_means; y++)
                temp.insert(k[j][y]);
            cout << temp.size() << " ==? " << num_local_means << endl;
            assert (temp.size() == num_local_means);
            */

            // compute means of local clusters
            vector< vector<double> > local_means (num_local_means, vector<double>(FIX_DIM, 0.0));
            vector<int> local_asgn (z.begin()+begin_i, z.begin()+end_i);
            vector<Instance*> local_data (data.begin()+begin_i,data.begin()+end_i);
            compute_means (local_data, local_asgn, FIX_DIM, local_means);
            assert (num_local_means == local_means.size());

            // pre-compute instances for global means 
            int num_global_means = global_means.size();
            vector<Instance*> temp_global_means (num_global_means, NULL);
            for (int p = 0; p < num_global_means; p ++) 
                temp_global_means[p] = vec2ins (global_means[p]);

            // pre-compute instances for local means 
            vector<Instance*> temp_local_means (num_local_means, NULL);
            for (int c = 0; c < num_local_means; c ++) 
                temp_local_means[c] = vec2ins (local_means[c]);

            for (int c = 0; c < num_local_means; c++) {
                // compute distance of local clusters to each global cluster
                num_global_means = global_means.size();
                vector<double> d_jcp (num_global_means, 0.0);
                double sum_d_ijc = 0.0; 
                for (int i = doc_lookup[j].first; i < doc_lookup[j].second; i ++) {
                    if (z[i] != c) continue;
                    double local_dist = df (data[i], temp_local_means[c], FIX_DIM);
                    sum_d_ijc += local_dist * local_dist;
                    for (int p = 0; p < num_global_means; p ++) {
                        double dist = df (data[i], temp_global_means[p], FIX_DIM);
                        d_jcp[p] += dist * dist;
                    }
                }
                int min_p = -1; double min_d_jcp = INF;
                for (int p = 0; p < num_global_means; p ++) 
                    if (d_jcp[p] < min_d_jcp) {
                        min_p = p;
                        min_d_jcp = d_jcp[p];
                    }
                assert (min_p >= 0);
                // cout << min_d_jcp << " " << lambda_global << " " << sum_d_ijc << endl;
                if (min_d_jcp > lambda_global + sum_d_ijc) {
                    global_means.push_back(local_means[c]); //  push mu_jc
                    temp_global_means.push_back(vec2ins (local_means[c]));
                    k[j][c] = num_global_means;
                    // cout << "global increment" << endl;
                } else {
                    k[j][c] = min_p;
                }
            }
            for (int c = 0; c < num_local_means; c ++) 
                delete temp_local_means[c];
            num_global_means = global_means.size();
            for (int p = 0; p < num_global_means; p ++) 
                delete temp_global_means[p];
        }
        // 6. for each global clusters,
        compute_assignment (global_asgn, k, z, tables);
        /*
        cout << "compute global means.." << endl;
        cout << "#global created: " << global_means.size() 
            << ", #global used: " << get_num_global_means(k);
            */
        compute_means (data, global_asgn, FIX_DIM, global_means);

        // 7. convergence?
        new_cost = compute_cost (data, global_means, k, z, lambdas, tables, df, FIX_DIM);

        if ( new_cost < objmin ) objmin = new_cost;
        objmin_trace << omp_get_wtime()-start_time << " " << objmin << endl;
        if (new_cost == last_cost)
            break;
        if (new_cost < last_cost) {
            last_cost = new_cost;
        } else {
            assert(false);    
        }
    }
    means = global_means;
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

    objmin_trace << "time objective" << endl;

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
    /*
    ofstream dist_mat_out ("dist_mat");
    dist_mat_out << mat_toString(dist_mat, N, N);
    dist_mat_out.close();
    */
    start_time = omp_get_wtime();
    double min_obj = INF;
    vector< vector<double> > min_means;
    for (int j = 0; j < nRuns; j ++) {
        vector< vector<double> > means;
        // inner-doc shuffle
        for (int d = 0; d < D; d++) {
            int begin_i = doc_lookup[d].first;
            int end_i = doc_lookup[d].second;
            random_shuffle(data.begin()+begin_i, data.begin()+end_i);
        }
        // between-doc shuffle
        vector<pair<int,int> > s_doc_lookup (doc_lookup);
        random_shuffle(s_doc_lookup.begin(), s_doc_lookup.end());
        vector<Instance*> s_data (N, NULL);
        int p = 0;
        for (int d = 0; d < D; d ++) {
            for (int i = s_doc_lookup[d].first; i < s_doc_lookup[d].second; i ++) {
                s_data[p] = data[i];
                p ++;
            }
        }
        lookup_tables.doc_lookup = &s_doc_lookup;
        double obj = HDP_MEANS (s_data, means, &lookup_tables, LAMBDAs, df, FIX_DIM);
        lookup_tables.doc_lookup = &doc_lookup;
        cerr << "###################################################" << endl;
        if (obj < min_obj) {
            min_obj = obj;
            min_means = means;
        }
    }
     
    /* Output objective */ 
    output_objective (min_obj);
    ofstream model_out ("opt_model");
    for (int i = 0; i < min_means.size(); i ++) {
        model_out << "mean[" << i << "] "; 
        for (int j = 0; j < min_means[i].size(); j ++) {
            model_out << min_means[i][j] << " " ;
        }
        model_out << endl;
    }
    model_out.close();
    /* Output cluster centroids */
   // output_model (W, &lookup_tables);
    /* Output assignment */
   // output_assignment (W, &lookup_tables);

    /* reallocation */
    mat_free (W, N, N);
    mat_free (dist_mat, N, N);
    objmin_trace.close();
}
