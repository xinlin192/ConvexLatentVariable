#include "HDP_convex.h"

#include <string>

ofstream ss_out ("../obj_vs_time_dp/wholesale/HDP-convex");
double objmin = INF;

/* algorithmic options */ 
#define EXACT_LINE_SEARCH  // comment this to use inexact search

/* dumping options */
// #define EXACT_LINE_SEARCH_DUMP
const double NOISE_EPS = 1e-2;
const double SPARSITY_TOL = 1e-5; //
const double FRANK_WOLFE_TOL = 1e-20; //
const double ADMM_EPS = 5e-3;  //
const double r = 1000000.0; // dummy penality rate
const double EPS = 1e-5;

/* Compute the mutual distance of input instances contained within "data" */

void frank_wolfe_solver (double ** dist_mat, double ** y, double ** z, double ** w, double rho, int R, int C, int FW_MAX_ITER, set<int>& col_active_set) {
    // cout << "within frank_wolfe_solver" << endl;
    // STEP ONE: compute gradient mat initially
    vector< set< pair<int, double> > > actives (R, set<pair<int,double> >());
    vector< priority_queue< pair<int,double>, vector< pair<int,double> >, Int_Double_Pair_Dec> > pqueues (R, priority_queue< pair<int,double>, vector< pair<int,double> >, Int_Double_Pair_Dec> ());
    //#pragma omp parallel for
    for (int i = 0; i < R; i++) {
        for (set<int>::iterator it=col_active_set.begin();it != col_active_set.end(); ++it) {
            int j = *it;
            // double grad=0.5*dist_mat[i][j]+y[i][j]+rho*(w[i][j]-z[i][j]); 
            double grad=dist_mat[i][j]+y[i][j]+rho*(w[i][j]-z[i][j]); 
            if (w[i][j] > 1e-10) 
                actives[i].insert(make_pair(j,grad));
            else 
                pqueues[i].push(make_pair(j, grad));
        }
    }
    // STEP TWO: iteration solve each row 
    int k = 0;  // iteration number
    vector<bool> is_fw_opt_reached (R, false);
    set<pair<int,double> >::iterator it;
    // cout << "within frank_wolfe_solver: start iteration" << endl;
    while (k < FW_MAX_ITER) { 
        // compute new active atom: can be in active set or not
        vector< pair<int, double> > s (R, pair<int,double>());
        vector<bool> isInActives (R, false);
        for (int i = 0; i < R; i++) {
            if (is_fw_opt_reached[i]) continue;
            if (pqueues[i].size() <= 0) continue;
            pair<int,double> tmp_pair = pqueues[i].top();
            s[i].first = tmp_pair.first;
            s[i].second = tmp_pair.second;
            // cout << s[i].first << ":" << s[i].second << endl;
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) {
                // take the minimal of each row
                if (it->second < s[i].second) {
                    isInActives[i] = true;
                    s[i].first = it->first;
                    s[i].second = it->second;
                }
            }
            // compute gamma: inexact or exact
            double gamma; // step size of line search
#ifdef EXACT_LINE_SEARCH
            double sum1=0.0, sum2=0.0, sum3=0.0, sum4=0.0;
            // gamma* = (sum1 + sum2 + sum3) / sum4, where
            // sum1 = 1/2 sum_n sum_k (w - s)_nk * || x_n - mu_k ||^2
            // sum2 = sum_n sum_k y_nk (w - s)_nk
            // sum3 = - rho * sum_n sum_k  (w - z) (w-s)
            // sum4 = sum_n sum_k rho * (s - w)^2
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) {
                double w_minus_s;
                double w_minus_z = w[i][it->first] - z[i][it->first];
                if (it->first == s[i].first) {
                    w_minus_s = w[i][it->first]-1.0;
                } else {
                    w_minus_s = w[i][it->first];
                }
                // sum1 += 0.5* w_minus_s * (dist_mat[i][it->first] -r);
                sum1 +=  w_minus_s * (dist_mat[i][it->first] -r);
                sum2 += y[i][it->first] * w_minus_s;
                sum3 += rho * w_minus_s * w_minus_z;
                sum4 += rho * w_minus_s * w_minus_s; 
            }
            if (!isInActives[i]) {
                // sum1 +=  0.5*(-1.0) * (dist_mat[i][s[i].first] - r);
                sum1 += (-1.0) * (dist_mat[i][s[i].first] - r);
                sum2 += y[i][it->first] * (-1.0);
                sum3 += rho * (-1.0) * (w[i][s[i].first]-z[i][s[i].first]);
                sum4 += rho;
            }

            if (fabs(sum4) > 0) {
                gamma = (sum1 + sum2 + sum3) / sum4;
#ifdef EXACT_LINE_SEARCH_DUMP
                cout << "[exact] i=" << i ;
                cout << ",k=" << k;
                cout << ",sum1="<< sum1;
                cout << ",sum2="<< sum2;
                cout << ",sum3="<< sum3;
                cout << ",sum4="<< sum4;
                cout << ",gamma="<< gamma;
                cout << endl;
#endif
                gamma = max(gamma, 0.0);
                gamma = min(gamma, 1.0);
            } else {
                gamma = 0.0;
                is_fw_opt_reached[i] = true;
            }
#else
            gamma = 2.0 / (k+2.0);
#endif
            // update w
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) 
                w[i][it->first] *= (1-gamma);
            w[i][s[i].first] += gamma;
            // update new actives 
            set< pair<int, double> > temp;
            if (!isInActives[i]) {
                actives[i].insert(pqueues[i].top());
                pqueues[i].pop();
            }
            double new_grad;
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) {
                int j = it->first;
                // new_grad=0.5*dist_mat[i][j]+y[i][j]+rho*(w[i][j]-z[i][j]); 
                new_grad=dist_mat[i][j]+y[i][j]+rho*(w[i][j]-z[i][j]); 
                temp.insert (make_pair(it->first, new_grad));
            }
            actives[i].swap(temp);
        }
        k ++;
    }
}
void skyline (double** wout, double**wbar, int R_start, int R_end, int C, double lambda, set<int> col_active_sets) {
    vector< vector< double > > alpha_vec (C, vector<double>());
    vector< int > num_alpha_elem (C, 0);
    for (int i = R_start; i < R_end; i ++) {
        set<int>::iterator it;
        for (it=col_active_sets.begin();it!=col_active_sets.end();++it) {
            int j = *it;
            if (wbar[i][j] > SPARSITY_TOL) {
                alpha_vec[j].push_back (abs(wbar[i][j]));
                ++ num_alpha_elem[j];
            }
        }
    }
    vector<double> max_term (C, -1e50);
    vector<double> separator (C, 0.0);
    int R = R_end - R_start;
    //#pragma omp parallel for
    for (int j = 0; j < C; j ++) {
        if (num_alpha_elem[j] == 0) continue;
        // 2. sorting
        std::sort (alpha_vec[j].begin(), alpha_vec[j].end(), double_dec_comp);
        // 3. find mstar
        int mstar = 0; // number of elements support the sky
        double new_term, sum_alpha = 0.0;
        for (int i = 0; i < num_alpha_elem[j]; i ++) {
            sum_alpha += alpha_vec[j][i];
            new_term = (sum_alpha - lambda) / (i + 1.0);
            if ( new_term > max_term[j] ) {
                separator[j] = alpha_vec[j][i];
                max_term[j] = new_term;
                mstar = i;
            }
        }
        if (max_term[j] < 0) 
            max_term[j] = (sum_alpha - lambda) / R;
    }
    for (int i = R_start; i < R_end; i ++) {
        // 4. assign closed-form solution to wout
        set<int>::iterator it;
        for (it=col_active_sets.begin();it!=col_active_sets.end();++it) {
            int j = *it;
            if ( max_term[j] < 0 ) {
                wout[i][j] = 0.0;
                continue;
            }
            double wbar_val = wbar[i][j];
            if ( abs(wbar_val) >= separator[j] ) 
                wout[i][j] = max_term[j];
            else 
                // its ranking is above m*, directly inherit the wbar
                wout[i][j] = max(wbar_val,0.0);
        }
    }
}
void group_lasso_solver (double** y, double** z, double** w, double rho, vector<double>& lambda, Lookups *tables, set<int> col_active_sets) {
    int R = tables->nWords;
    int C = tables->nWords;
    vector< pair<int,int> > doc_lookup = *(tables->doc_lookup);
    double global_lambda = lambda[0];
    double local_lambda = lambda[1];

    double** wlocal = mat_init (R, C);
    double** wbar = mat_init (R, C);
    for (int i = 0; i < R; i ++) 
        for (int j = 0; j < C; j ++) 
            wbar[i][j] = (rho * z[i][j] - y[i][j]) / rho;
    // extend the group lasso solver to both local and global
    for (int d = 0; d < tables->nDocs; d++) {
        int R_start = doc_lookup[d].first;
        int R_end = doc_lookup[d].second;
        skyline (wlocal, wbar, R_start, R_end, C, local_lambda, col_active_sets);
    }
    skyline (w, wlocal, 0, R, C, global_lambda, col_active_sets);
    
    mat_free (wlocal, R, C);
    mat_free (wbar, R, C);
}

double overall_objective (double ** dist_mat, vector<double>& lambda, int R, int C, double ** z, Lookups* tables) {
    int D = tables->nDocs;
    vector< pair<int,int> > doc_lookup = *(tables->doc_lookup);
    // STEP ONE: compute 
    //     loss = sum_i sum_j z[i][j] * dist_mat[i][j]
    double normSum = 0.0;
    for (int i = 0; i < R; i ++) 
        for (int j = 0; j < C; j ++) 
            normSum += z[i][j] * dist_mat[i][j];
    double loss =  normSum;
    // STEP TWO: compute dummy loss
    // sum4 = r dot (1 - sum_k w_nk) -> dummy
    double * temp_vec = new double [R];
    mat_sum_row (z, temp_vec, R, C);
    double dummy_penalty = 0.0;
    for (int i = 0; i < R; i ++) 
        dummy_penalty += r * max(1 - temp_vec[i], 0.0);
    delete[] temp_vec;
    // STEP THREE: compute group-lasso regularization
    double global_lasso = lasso_objective(z, lambda[0], 0, R, C);
    double sum_local_lasso = 0.0;
    vector<double> local_lasso (D, 0.0);
    for (int d = 0; d < D; d++) {
        local_lasso[d] = lasso_objective(z, lambda[1], doc_lookup[d].first, doc_lookup[d].second, C); 
        sum_local_lasso += local_lasso[d];
    }
    double overall = loss + global_lasso + sum_local_lasso + dummy_penalty;
    cerr << "loss: " << loss << ", dummy=" << dummy_penalty
        << ", global_lasso=" << global_lasso 
        << ", sum_local_lasso=" << sum_local_lasso 
        << ", overall=" << overall
        << endl;
    return loss + global_lasso + sum_local_lasso;
}


void cvx_hdp_medoids (double ** dist_mat, int fw_max_iter, vector<double>& lambda, double ** W, int ADMM_max_iter, int SS_PERIOD, Lookups * tables) {
    int N = tables->nWords;
    int D = tables->nDocs;
    //vector< pair<int,int> > doc_lookup = *(tables->doc_lookup);
    // parameters 
    double alpha = 0.1;
    double rho = 1;
    ss_out << "Time Objective" << endl;
    double cputime = 0;
    double prev = omp_get_wtime();
    // iterative optimization 
    double error = INF;
    double ** wone = mat_init (N, N);
    double ** wtwo = mat_init (N, N);
    double ** yone = mat_init (N, N);
    double ** ytwo = mat_init (N, N);
    double ** z = mat_init (N, N);
    double ** z_old = mat_init (N, N);

    // variables for shriking method
    set<int> col_active_sets;
    // set initial active_set as all elements
    for (int i = 0; i < N; i++) 
        col_active_sets.insert(i);

    cputime += omp_get_wtime() - prev;
    ss_out << cputime << " " << 0 << endl;
    prev = omp_get_wtime();
    int iter = 0; 
    bool no_active_element = false, admm_opt_reached = false;
    while ( iter < ADMM_max_iter ) { 

        // STEP ONE: resolve w_1 and w_2
        frank_wolfe_solver (dist_mat, yone, z, wone, rho, N, N, fw_max_iter, col_active_sets);
        group_lasso_solver (ytwo, z, wtwo, rho, lambda, tables, col_active_sets);

        // STEP TWO: update z by averaging w_1 and w_2
        // STEP THREE: update the y_1 and y_2 by w_1, w_2 and z
        set<int>::iterator it;
        for (it=col_active_sets.begin();it!=col_active_sets.end();++it) {
            int j = *it;
            for (int i = 0; i < N; i ++) {
                z[i][j] = (wone[i][j] + wtwo[i][j]) / 2.0;
                yone[i][j] += alpha * (wone[i][j] - z[i][j]) ;
                ytwo[i][j] += alpha * (wtwo[i][j] - z[i][j]) ;
            }
        }

        // Shrinking Method:
        // STEP ONE: reduce number of elements considered in next iteration
        vector<int> col_to_shrink;
        for (it=col_active_sets.begin();it!=col_active_sets.end();++it) {
            int j = *it;
            bool is_primal_shrink=false, is_dual_shrink=false;
            int i;
            for (i = 0; i < N; i++) {
                // (A) primal shrinking:
                if ( z[i][j] > ADMM_EPS || z_old[i][j] > ADMM_EPS ) break;
                // (B) dual shrinking:
                if ( wone[i][j] > ADMM_EPS || wtwo[i][j] > ADMM_EPS ) break;
            }
            // cache index of element to be removed
            if ( i == N ) 
                col_to_shrink.push_back(j);
        }
        
        // remove shrinked row/column from row/column active sets
        int num_col_to_shrink = col_to_shrink.size();
        for (int s = 0; s < num_col_to_shrink; s ++) {
            int j_shr = col_to_shrink[s];
            col_active_sets.erase(j_shr);
            for(int i=0;i<N;i++){
                wone[i][j_shr] = 0.0;
                wtwo[i][j_shr] = 0.0;
                z[i][j_shr] = 0.0;
                z_old[i][j_shr] = 0.0;
            }
        }
        // update z_old
        for (it=col_active_sets.begin();it!=col_active_sets.end();++it) {
            int j = *it;
            for (int i = 0; i < N; i++) 
                z_old[i][j] = z[i][j];
        }
        // count number of active elements
        int num_active_cols = col_active_sets.size();
        // STEP FOUR: trace the objective function
        if ( (iter+1) % SS_PERIOD == 0) {
            cputime += omp_get_wtime() - prev;
            error = overall_objective (dist_mat, lambda, N, N, wone, tables);
            cerr << "[Overall] iter = " << iter 
                << ", num_active_cols: " << num_active_cols
                << ", Loss Error: " << error << endl;
            if (error < objmin) objmin = error;
            ss_out << cputime << " " << objmin << endl;
            prev = omp_get_wtime();
        }
        int num_active_elements=N*num_active_cols;
        // STEP TWO: consider to open all elements to check optimality
        /*
        if (num_active_elements == 0 && !no_active_element) {
            no_active_element = true;
            // open all elements to verify the result
            cout << "open all elements for optimality checking!" << endl;
            col_active_sets.clear();
            for (int i = 0; i < N; i++) {
                col_active_sets.insert(i);
            }
        } else if (num_active_elements == 0 && no_active_element) 
            admm_opt_reached = true;
        else if (num_active_elements > 0 && no_active_element) {
            no_active_element = false;
            cout << "fail to reach global ADMM optima!" << endl;
        }
        */
        iter ++;
    }

    // STEP FIVE: memory recollection
    mat_free (wone, N, N);
    mat_free (wtwo, N, N);
    mat_free (yone, N, N);
    mat_free (ytwo, N, N);
    mat_free (z_old, N, N);
    // STEP SIX: put converged solution to destination W
    mat_copy (z, W, N, N);
    mat_free (z, N, N);
    ss_out.close();
}

// entry main function
int main (int argc, char ** argv) {
    if (argc < 6) {
        cerr << "Usage: " << endl;
        cerr << "\tcvx_hdp_medoids [word_dataFile] [lambda_global] [lambda_local] [FW_MAX_ITER] [ADMM_MAX_ITER]" << endl;
        exit(-1);
    }

    // PARSE arguments
    char* dataFile = argv[1];
    vector<double> LAMBDAs (2, 0.0);
    LAMBDAs[0] = atof(argv[2]); // lambda_global
    LAMBDAs[1] = atof(argv[3]); // lambda_local
    int FW_MAX_ITER = atoi(argv[4]);
    int ADMM_MAX_ITER = atoi(argv[5]);
    int SS_PERIOD = 20;

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
    double** noise_mat = mat_init(N, N);
    compute_dist_mat (data, dist_mat, N, FIX_DIM, df, true);
    //adding perturbation
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            noise_mat[i][j] = NOISE_EPS * ((double)rand()/RAND_MAX);
        }
    }

    ofstream dmat_out ("dist_mat");
    
    dmat_out << mat_toString (dist_mat, N, N);
    dmat_out.close();
    cerr << "dist_mat output finished.." << endl;
    
    mat_add( dist_mat, noise_mat, dist_mat, N, N);
    cvx_hdp_medoids (dist_mat, FW_MAX_ITER, LAMBDAs, W, ADMM_MAX_ITER, SS_PERIOD, &lookup_tables);
    mat_sub( dist_mat, noise_mat, dist_mat, N, N); 

   // compute the mean and update objective to compare with DP_MEANS
    cerr << "==================================================" << endl; 
    cerr << "Computing the mean of derived group members ..." << endl;
    vector<int> centroids;
    get_all_centroids(W, &centroids, N, N);
    int nCentroids = centroids.size();
    vector< vector<int> > members (nCentroids, vector<int>());
    for (int c = 0; c < nCentroids; c++) {
        int j = centroids[c];
        // cout << "centroid: " << j << endl;
        for (int i = 0; i < N; i++) {
            if (W[i][j] < 1.01 && W[i][j] > 0.99) {
                members[c].push_back(i);
               // cerr << "member: "<< i << endl;
            }
        }
    }
    // compute the means
    vector< vector<double> > means (nCentroids, vector<double>(FIX_DIM, 0.0));
    for (int c = 0; c < nCentroids; c++) {
        int numMembers = members[c].size();
        for (int x = 0; x < numMembers; x ++) {
            int i = members[c][x];
            Instance* ins = data[i];
            int fea_size = ins->fea.size();
            for (int f = 0; f < fea_size; f++) 
                means[c][ins->fea[f].first-1] += ins->fea[f].second;
        }
        for (int f = 0; f < D; f++) 
            means[c][f] = means[c][f] / numMembers;
    }
    // compute distance
    double global_reg = LAMBDAs[0] * nCentroids;
    double local_reg = 0.0;
    for (int d = 0; d < D; d++) {
        int begin_i = doc_lookup[d].first;
        int end_i = doc_lookup[d].second;
        vector<int> local_centroids;
        double ** temp_W = mat_init (end_i-begin_i,N);
        for (int i = begin_i; i < end_i; i++) 
            for (int j = 0; j < N; j++)
                temp_W[i-begin_i][j] = W[i][j];
        get_all_centroids(temp_W, &local_centroids, end_i-begin_i, N);
        int local_nCentroids = local_centroids.size();
        local_reg += LAMBDAs[1] * local_nCentroids;
        mat_free(temp_W, end_i-begin_i,N);
    }
    double sum_means_loss = 0.0;
    for (int c = 0; c < nCentroids; c++) {
        double mean_loss = 0.0;
        int numMembers = members[c].size();
        Instance* mean_ins = new Instance(c+100000);
        for (int f = 0; f < D; f++) 
            mean_ins->fea.push_back(make_pair(f+1, means[c][f]));
        for (int x = 0; x < numMembers; x ++) {
            int i = members[c][x];
            Instance* ins = data[i];
            double dist = df(mean_ins, ins, D);
            mean_loss +=  dist * dist;
        }
        delete mean_ins;
        cerr << "c=" << centroids[c] << ", mean_loss: " << mean_loss << endl;
        sum_means_loss += mean_loss;
    }
    // output distance
    cerr << "Total_means_loss: " << sum_means_loss
        << ", Global_reg : " << global_reg
        << ", Local_reg : " << local_reg << endl;
    cerr << "Total_Error: " << sum_means_loss + global_reg + local_reg << endl;

    /* Output objective */ 
    output_objective (dist_mat, W, &lookup_tables, r, LAMBDAs);
    /* Output cluster centroids */
    output_model (W, &lookup_tables, LAMBDAs);
    /* Output assignment */
    output_assignment (W, &lookup_tables, LAMBDAs);
    
    /* reallocation */
    mat_free (W, N, N);
    mat_free (dist_mat, N, N);
}
