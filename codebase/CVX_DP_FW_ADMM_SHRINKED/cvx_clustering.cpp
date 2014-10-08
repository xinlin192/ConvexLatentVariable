#include "cvx_clustering.h"

/* algorithmic options */ 
#define EXACT_LINE_SEARCH  // comment this to use inexact search

/* dumping options */
// #define EXACT_LINE_SEARCH_DUMP

const double FRANK_WOLFE_TOL = 1e-20;
const double ADMM_EPS = 0.02;
const double SPARSITY_TOL = 1e-5;
const double r = 1000000.0;

void frank_wolfe_solver (double ** dist_mat, double ** yone, double ** zone, double ** wone, double rho, int N, int K, set<int>& col_active_set) {
    // cout << "within frank_wolfe_solver" << endl;
    // STEP ONE: compute gradient mat initially
    vector< set< pair<int, double> > > actives (N, set<pair<int,double> >());
    vector< priority_queue< pair<int,double>, vector< pair<int,double> >, Compare> > pqueues (N, priority_queue< pair<int,double>, vector< pair<int,double> >, Compare> ());
    for (int i = 0; i < N; i++) {
        for (set<int>::iterator it=col_active_set.begin();it != col_active_set.end(); ++it) {
            int j = *it;
            double grad=0.5*dist_mat[i][j]+yone[i][j]+rho*(wone[i][j]-zone[i][j]); 
            if (wone[i][j] > 1e-10) 
                actives[i].insert(make_pair(j,grad));
            else 
                pqueues[i].push(make_pair(j, grad));
        }
    }
    // STEP TWO: iteration solve each row 
    int k = 0;  // iteration number
    vector<bool> is_fw_opt_reached (N, false);
    set<pair<int,double> >::iterator it;
    // cout << "within frank_wolfe_solver: start iteration" << endl;
    while (k < K) { // TODO: change to use portional criteria
        // compute new active atom: can be in active set or not
        vector< pair<int, double> > s (N, pair<int,double>());
        vector<bool> isInActives (N, false);
        for (int i = 0; i < N; i++) {
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
                double w_minus_z = wone[i][it->first] - zone[i][it->first];
                if (it->first == s[i].first) {
                    w_minus_s = wone[i][it->first]-1.0;
                } else {
                    w_minus_s = wone[i][it->first];
                }
                sum1 += 0.5 * w_minus_s * (dist_mat[i][it->first] -r);
                sum2 += yone[i][it->first] * w_minus_s;
                sum3 += rho * w_minus_s * w_minus_z;
                sum4 += rho * w_minus_s * w_minus_s; 
            }
            if (!isInActives[i]) {
                sum1 += 0.5 * (-1.0) * (dist_mat[i][s[i].first] - r);
                sum2 += yone[i][it->first] * (-1.0);
                sum3 += rho * (-1.0) * (wone[i][s[i].first]-zone[i][s[i].first]);
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
            // update wone
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) {
                wone[i][it->first] *= (1-gamma);
            }
            wone[i][s[i].first] += gamma;
            // update new actives 
            set<pair<int, double> > temp;
            if (!isInActives[i]) {
                actives[i].insert(pqueues[i].top());
                pqueues[i].pop();
            }
            double new_grad;
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) {
                int j = it->first;
                new_grad=0.5*dist_mat[i][j]+yone[i][j]+rho*(wone[i][j]-zone[i][j]); 
                temp.insert (make_pair(it->first, new_grad));
            }
            actives[i].swap(temp);
            // cout << "actives[" << i << "]: " << actives[i].size() << endl;
        }
        // cout << "within frank_wolfe_solver: next iteration" << endl;
        k ++;
    }
}


void skyline (double** wout, double** wbar, int R, int C, double lambda, set<int> col_active_sets) {
    vector< vector< double > > alpha_vec (C, vector<double>());
    vector< int > num_alpha_elem (C, 0);
    for (int i = 0; i < R; i ++) {
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
    for (int i = 0; i < R; i ++) {
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
void group_lasso_solver (double** y, double** z, double** w, int R, int C, double rho, double lambda, set<int> col_active_sets) {

    double** wbar = mat_init (R, C);
    mat_zeros (wbar, R, C);
    for (int i = 0; i < R; i ++) {
        for (int j = 0; j < C; j ++) {
            wbar[i][j] = (rho * z[i][j] - y[i][j]) / rho;
        }
    }
    skyline (w, wbar, R, C, lambda, col_active_sets);

    mat_free (wbar, R, C);
}

double overall_objective (double ** dist_mat, double lambda, int N, double ** z) {
    // N is number of entities in "data", and z is N by N.
    // z is current valid solution (average of w_1 and w_2)
    double sum = 0.0;
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            sum += z[i][j];
    // cerr << "sum=" << sum/N << endl;
    // STEP ONE: compute 
    //     loss = sum_i sum_j z[i][j] * dist_mat[i][j]
    double normSum = 0.0;
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            normSum += z[i][j] * dist_mat[i][j];
        }
    }
    double loss = 0.5 * (normSum);
    cout << "loss=" << loss;
    // STEP TWO: compute dummy loss
    // sum4 = r dot (1 - sum_k w_nk) -> dummy
    double * temp_vec = new double [N];
    mat_sum_row (z, temp_vec, N, N);
    double dummy_penalty=0.0;
    double avg=0.0;
    for (int i = 0; i < N; i ++) {
        avg += temp_vec[i];
        dummy_penalty += r * max(1 - temp_vec[i], 0.0) ;
    }
    cout << ", dummy= " << dummy_penalty;
    // STEP THREE: compute group-lasso regularization
    double * maxn = new double [N]; 
    for (int i = 0;i < N; i ++) { // Ian: need initial 
        maxn[i] = -INF;
    }
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            if ( fabs(z[i][j]) > maxn[j])
                maxn[j] = fabs(z[i][j]);
        }
    }
    double sumk = 0.0;
    for (int i = 0; i < N; i ++) {
        sumk += lambda*maxn[i];
    }
    double reg = sumk; 
    cout << ", reg=" << reg ;
    delete[] maxn;
    delete[] temp_vec;
    double overall = loss + reg + dummy_penalty;
    cout << ", overall=" <<  overall << endl;
    return loss + reg;
}

void cvx_clustering (double ** dist_mat, int fw_max_iter, int D, int N, double lambda, double ** W, int ADMM_max_iter, int SS_PERIOD) {
    // parameters 
    double alpha = 0.2;
    double rho = 1;
    ofstream ss_out ("plot_cputime_objective");
    ss_out << "Time Objective" << endl;
    double cputime = 0;
    clock_t prev = clock();
    // iterative optimization 
    double error = INF;
    double ** wone = mat_init (N, N);
    double ** wtwo = mat_init (N, N);
    double ** yone = mat_init (N, N);
    double ** ytwo = mat_init (N, N);
    double ** z = mat_init (N, N);
    double ** z_old = mat_init (N, N);
    double ** diffzero = mat_init (N, N);
    double ** diffone = mat_init (N, N);
    double ** difftwo = mat_init (N, N);
    mat_zeros (wone, N, N);
    mat_zeros (wtwo, N, N);
    mat_zeros (yone, N, N);
    mat_zeros (ytwo, N, N);
    mat_zeros (z, N, N);
    mat_zeros (z_old, N, N);
    mat_zeros (diffzero, N, N);
    mat_zeros (diffone, N, N);
    mat_zeros (difftwo, N, N);

    // variables for shriking method
    set<int> col_active_sets;
    // set initial active_set as all elements
    for (int i = 0; i < N; i++) 
        col_active_sets.insert(i);

    cputime += clock() - prev;
    ss_out << cputime << " " << 0 << endl;
    prev = clock();
    int iter = 0; // Ian: usually we count up (instead of count down)
    bool no_active_element = false, admm_opt_reached = false;
    while ( iter < ADMM_max_iter ) { // stopping criteria

        // STEP ONE: resolve w_1 and w_2
        frank_wolfe_solver (dist_mat, yone, z, wone, rho, N, fw_max_iter, col_active_sets);
        group_lasso_solver (ytwo, z, wtwo, N, N, rho, lambda, col_active_sets);

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

        // STEP FOUR: trace the objective function
        if (iter < 3 * SS_PERIOD || (iter+1) % SS_PERIOD == 0) {
            cputime += (double)(clock() - prev) / CLOCKS_PER_SEC;
            error = overall_objective (dist_mat, lambda, N, wone);
            cout << "[Overall] iter = " << iter 
                << ", Loss Error: " << error << endl;
            ss_out << cputime << " " << error << endl;
            prev = clock();
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
            for (int i = 0; i < N; i++) {
                z_old[i][j] = z[i][j];
            }
        }
        // count number of active elements
        int num_active_cols = col_active_sets.size();
        if ((iter+1) % SS_PERIOD == 0) {
            ;
            /*
            cout << "iter: " << iter;
            cout << ", num_active_cols: " << num_active_cols <<endl;
            */
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
    mat_free (diffone, N, N);
    mat_free (difftwo, N, N);
    mat_free (z_old, N, N);
    // STEP SIX: put converged solution to destination W
    mat_copy (z, W, N, N);
    mat_free (z, N, N);
    ss_out.close();
}

// entry main function
int main (int argc, char ** argv) {
    // exception control: illustrate the usage if get input of wrong format
    if (argc < 6) {
        cerr << "Usage: cvx_clustering [dataFile] [fw_max_iter] [ADMM_max_iter] [lambda] [sc_period]" << endl;
        cerr << "Note: dataFile must be scaled to [0,1] in advance." << endl;
        exit(-1);
    }

    // parse arguments
    char * dataFile = argv[1];
    int fw_max_iter = atoi(argv[2]);
    int ADMM_max_iter = atoi(argv[3]);
    double lambda_base = atof(argv[4]);
    int screenshot_period = atoi(argv[5]);

    // read in data
    int FIX_DIM;
    Parser parser;
    vector<Instance*>* pdata;
    vector<Instance*> data;
    pdata = parser.parseSVM(dataFile, FIX_DIM);
    data = *pdata;

    // explore the data 
    int dimensions = -1;
    int N = data.size(); // data size
    for (int i = 0; i < N; i++) {
        vector< pair<int,double> > * f = &(data[i]->fea);
        int last_index = f->size()-1;
        if (f->at(last_index).first > dimensions) {
            dimensions = f->at(last_index).first;
        }
    }
    assert (dimensions == FIX_DIM);

    int D = dimensions;
    cerr << "D = " << D << endl; // # features
    cerr << "N = " << N << endl; // # instances
    cerr << "lambda = " << lambda_base << endl;
    cerr << "r = " << r << endl;
    cerr << "Screenshot period = " << screenshot_period << endl;
    int seed = time(NULL);
    srand (seed);
    cerr << "seed = " << seed << endl;

    double lambda = lambda_base;

    // pre-compute distance matrix
    dist_func df = L2norm;
    double ** dist_mat = mat_init (N, N);
    //  double ** dist_mat = mat_read (dmatFile, N, N);
    mat_zeros (dist_mat, N, N);
    compute_dist_mat (data, dist_mat, N, D, df, true); 

    // Run sparse convex clustering
    double ** W = mat_init (N, N);
    mat_zeros (W, N, N);
    cvx_clustering (dist_mat, fw_max_iter, D, N, lambda, W, ADMM_max_iter, screenshot_period);

    // Output cluster
    output_objective(clustering_objective (dist_mat, W, N));
    /* Output cluster centroids */
    output_model (W, N);
    /* Output assignment */
    output_assignment (W, data, N);

    /*
    // output the examplar 
    // compute the mean and update objective to compare with DP_MEANS
    vector< int > centroids ();
    get_all_centroids(W, &centroids, N, N);
    int nCentroids = centroids.size();
    vector< vector<int> > members (nCentroids, vector<int> ());
    for (int j = 0; j < nCentroids; i++) {
        for (int i = 0; i < N; i++) {
            if (W[i][j] < 1.001 && W[i][j] > 0.999)
                members[j].push_back(i);
        }
    }
    // compute the means
    vector< vector<double> > means (nCentroids, vector<double>());
    for (int j = 0; j < nCentroids; j++) {
        int numMembers = members[j].size();
        for (int x = 0; x < numMembers; x ++) {
            int i = members[j][x];
            for (int )
        }
    }
    */

    // compute distance
    
    // output distance

    /* reallocation */
    mat_free (dist_mat, N, N);
    mat_free (W, N, N);
}
