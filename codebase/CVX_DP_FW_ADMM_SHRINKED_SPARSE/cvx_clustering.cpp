/*############################################################### 
## MODULE: cvx_clustering.cpp
## VERSION: 1.0 
## SINCE 2014-06-14
## AUTHOR:
##     Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
## DESCRIPTION: 
##
################################################################# 
## Edited by MacVim
## Class Info auto-generated by Snippet 
################################################################*/

#include "cvx_clustering.h"
#include <cassert>
#include <queue>
#include <time.h>

#include "../util.h"

/* algorithmic options */ 
#define EXACT_LINE_SEARCH  // comment this to use inexact search

/* dumping options */
// #define FRANK_WOLFE_DUMP
// #define EXACT_LINE_SEARCH_DUMP
// #define BLOCKWISE_DUMP
// #define SUBPROBLEM_DUMP 

const double FRANK_WOLFE_TOL = 1e-20;
const double ADMM_EPS = 1e-2;
typedef double (* dist_func) (Instance*, Instance*, int); 
const double r = 10000.0;

/* \lambda_g \sumk \maxn |\wnk| */
double compute_group_lasso (Esmat* w, double lambda) {
    if (w->val.size() == 0) 
        return 0.0;
    Esmat* maxn = esmat_init (1, w->nCols);
    Esmat* sumk = esmat_init (1, 1);
    esmat_max_over_col (w, maxn);
    esmat_sum_row (maxn, sumk);
    double lasso = -INF;
    if (sumk->val.size() > 0)
        lasso = lambda * sumk->val[0].second; 
    else 
        lasso = 0.0;
    esmat_free (sumk);
    esmat_free (maxn);
    return lasso;
}

void subproblem_objective (double** dist_mat, Esmat* y, Esmat* z, Esmat* w, double rho, int N, double lambda) {
    Esmat* diff = esmat_init (N, N);
    esmat_zeros (diff);
    // reg = 0.5 * sum_k max_n | w_nk |  -> group-lasso
    double group_lasso = compute_group_lasso(w, lambda); 
    // loss = 0.5 * sum_n sum_k (w_nk * d^2_nk) -> loss
    double loss = 0.5 * esmat_frob_prod (dist_mat, w);
    // linear = y_1^T dot (w_1 - z) -> linear
    esmat_sub (w, z, diff); // temp = w_1 - z_1
    double linear = esmat_frob_prod (y, diff);
    // quadratic = 0.5 * rho * || w_1 - z_1 ||^2 -> quadratic
    double quadratic = 0.5 * rho * esmat_frob_norm (diff);
    // dummy = r dot (1 - sum_k w_nk) -> dummy
    dummy_penalty = esmat_compute_dummy (w, r);
    // double total = loss+linear+quadratic+dummy_penalty;
    cout << "(loss, lasso, linear, quadratic, dummy, total) = (" 
        << loss << ", " << lasso << ", " << linear << ", " <<
        quadratic << ", " << dummy_penalty << ", " << total <<  ")" << endl;
    esmat_free (diff);
    // return total;
}
double overall_objective (double** dist_mat, double* lambda, int N, Esmat* z) {
    // N is number of entities in "data", and z is N by N.
    // z is current valid solution (average of w_1 and w_2)
    // STEP ONE: compute 
    //     loss = sum_i sum_j z[i][j] * dist_mat[i][j]
    double loss = 0.5 * esmat_frob_prod(dist_mat, z);
    cout << "loss=" << loss;
    // STEP TWO: compute dummy loss
    // sum4 = r dot (1 - sum_k w_nk) -> dummy
    double dummy_penalty = esmat_compute_dummy (z, r);
    cout << ", dummy= " << dummy_penalty;
    // STEP THREE: compute group-lasso regularization
    double reg = compute_group_lasso(z, lambda); 
    cout << ", reg=" << reg ;
    double overall = loss + reg + dummy_penalty;
    cout << ", overall=" <<  overall << endl;
    return loss + reg;
}

/* Compute the mutual distance of input instances contained within "data" */
void compute_dist_mat (vector<Instance*>& data, Esmat* dist_mat, int N, int D, dist_func df, bool isSym) {
    for (int j = 0; j < N; j ++) {
        Instance * muj = data[j];
        for (int i = 0; i < N; i ++) {
            Instance * xi = data[i];
            double dist_value = df (xi, muj, D);
            dist_mat->val.push_back(make_pair(j*N+i,dist_value));
        }
    }
}
void compute_dist_mat (vector<Instance*>& data, double** dist_mat, int N, int D, dist_func df, bool isSym) {
    for (int i = 0; i < N; i ++) {
        Instance * xi = data[i];
        for (int j = 0; j < N; j ++) {
            Instance * muj = data[j];
            double dist_value = df (xi, muj, D);
            dist_mat[i][j] = dist_value;
        }
    }
}
void frank_wolfe_solver (double** dist_mat, Esmat* yone, Esmat* zone, Esmat* wone, double rho, int N, int K, map<int,int>& col_active_map) {
    // cout << "within frank_wolfe_solver" << endl;
    // STEP ONE: compute gradient mat initially
    vector< set< pair<int, double> > > actives (N, set<pair<int,double> >());
    vector< priority_queue< pair<int,double>, vector< pair<int,double> >, Compare> > pqueues (N, priority_queue< pair<int,double>, vector< pair<int,double> >, Compare> ());
    // compute N by K active gradient matrix
    Esmat* assist = esmat_init (N, N);
    esmat_sub (wone, zone, assist);
    esmat_scalar_mult (rho, assist);
    esmat_add (assist, yone);
    int num_active_cols = col_active_map.size();
    double** grad = mat_init (N, num_active_cols);
    mat_zeros (grad, N, num_active_cols);
    // set up actives for wone > 1e-10
    int num_active_elements = wone->val.size();
    for (int i = 0; i < num_active_elements; i ++) {
        int esmat_index = wone->val[i].first;
        int row_index = esmat_index % wone->nCols;
        int col_index = esmat_index / wone->nCols;
        double value = grad[row_index][col_active_map[col_index]];
        actives[row_index].insert(make_pair(col_index, value));
    }
    // set up queue with ruling out inactive columns
    for (int i = 0; i < N; i++) {
        int j = *it;
        double grad=0.5*dist_mat[i][j]+yone[i][j]+rho*(wone[i][j]-zone[i][j]); 
        // pqueues[i].push(make_pair(j, grad));

    }
    // STEP TWO: iteration solve each row 
    int k = 0;  // iteration number
    vector<bool> is_fw_opt_reached (N, false);
    set<pair<int,double> >::iterator it;
    // cout << "within frank_wolfe_solver: start iteration" << endl;
    while (k < K) { 
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
                cout << ",sum1=" << sum1;
                cout << ",sum2=" << sum2;
                cout << ",sum3=" << sum3;
                cout << ",sum4=" << sum4;
                cout << ",gamma=" << gamma;
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
void group_lasso_solver (Esmat* Y, Esmat* Z, Esmat* w, double RHO, double lambda) {
    // STEP ONE: compute the optimal solution for truncated problem
    Esmat* wbar = esmat_init (Z);
    Esmat* temp = esmat_init (Z);
    esmat_scalar_mult (RHO, Z, temp); // wbar = RHO * z
    // cout << wbar->nRows << "," << wbar->nCols << endl;
    // cout << Y->nRows << "," << Y->nCols << endl;
    esmat_sub (temp, Y, wbar); // wbar = RHO * z - y
    esmat_scalar_mult (1.0/RHO, wbar); // wbar = (RHO * z - y) / RHO
    esmat_free (temp);
#ifdef GROUP_LASSO_DEBUG
    cout << "[wbar]" << endl;
    cout << esmat_toString(wbar);
    cout << "lambda: " << lambda << endl;
#endif
    // STEP TWO: find the closed-form solution for second subproblem
    int SIZE = wbar->val.size();
    int R = wbar->nRows; int C = wbar->nCols;

    if (wbar->val.size() == 0) {
        // no need to solve all-zero matrix w
        return ;
    }
    // i is index of element in esmat->val, j is column index
    int i = 0; int j = 0;
    int begin_idx, end_idx;
    int col_es_begin = 0;
    vector< pair<int,double> > alpha_vec;
    while ((i < SIZE && j < C)) {
        begin_idx = j * R;
        end_idx = (j+1) * R;
        int esWBAR_index = wbar->val[i].first;
        // cout << "i: " << i << " , j: " << j << endl;
        if (esWBAR_index >= end_idx) {
            int nValidAlpha = alpha_vec.size();
            // cout << "nValidAlpha: " << nValidAlpha << endl;
            if (nValidAlpha == 0) {
                j ++;
                continue;
            }
            // a) sort existing temp_vec
            std::sort (alpha_vec.begin(), alpha_vec.end(), pair_Second_Elem_Comparator);
            // b) find mstar
            int mstar = 0; // number of elements supporting the sky
            double separator; 
            double max_term = -INF, new_term;
            double sum_alpha = 0.0;
            for (int v = 0; v < nValidAlpha; v ++) {
                sum_alpha += alpha_vec[v].second;
                new_term = (sum_alpha - lambda) / (v + 1.0);
                // cout << "new_term: " << new_term << endl;
                if ( new_term > max_term ) {
                    separator = alpha_vec[v].second;
                    max_term = new_term;
                    ++ mstar;
                }
            }
            // cout << "mstar: " << mstar << ", max_term: " << max_term << endl;
            // c) assign closed-form solution of current column to w
            if (mstar <= 0) {
                ; // this column of w is all-zero, hence we do nothing for that 
            } else {
                for (int esi = col_es_begin; wbar->val[esi].first < end_idx; esi ++) {
                    double pos = wbar->val[esi].first;
                    double value = wbar->val[esi].second;
                    if (fabs(value) >= separator) 
                        w->val.push_back(make_pair(pos, max(max_term, 0.0)));
                    else 
                        // w->val.push_back(make_pair(pos, max(value, 0.0)));
                        w->val.push_back(make_pair(pos, value));
                }
            }
            // d) clear all elements in alpha_vec 
            alpha_vec.clear();
            // e) push current element to the cleared alpha_vec
            double value = wbar->val[i].second;
            alpha_vec.push_back (make_pair(esWBAR_index % R, fabs(value)));
            // f) go to operate next element and next column
            col_es_begin = i;
            ++ i; ++ j;
        } else if (esWBAR_index >= begin_idx) {
            // a) push current element to the cleared temp_vec
            double value = wbar->val[i].second;
            alpha_vec.push_back (make_pair(esWBAR_index % R, fabs(value)));
            // b) go to operate next element with fixed column index (j)
            ++ i; 
        } else { // impossible to occur
            assert (false);
        }
    }
    if (alpha_vec.size() > 0) {
        // a) sort existing temp_vec
        std::sort (alpha_vec.begin(), alpha_vec.end(), pair_Second_Elem_Comparator);
        // b) find mstar
        int mstar = 0; // number of elements supporting the sky
        double separator;
        double max_term = -INF, new_term;
        double sum_alpha = 0.0;
        int nValidAlpha = alpha_vec.size();
        for (int v = 0; v < nValidAlpha; v ++) {
            sum_alpha += alpha_vec[v].second;
            new_term = (sum_alpha - lambda) / (v + 1.0);
            if ( new_term > max_term ) {
                separator = alpha_vec[v].second;
                max_term = new_term;
                ++ mstar;
            }
        }
        // cout << "mstar: " << mstar << ", max_term: " << max_term << endl;
        // c) assign closed-form solution of current column to w
        if (nValidAlpha == 0 || mstar <= 0) {
            ; // this column of w is all-zero, hence we do nothing for that 
        } else {
            for (int esi = col_es_begin; esi < SIZE; esi ++) {
                double pos = wbar->val[esi].first;
                double value = wbar->val[esi].second;
                if (fabs(value) >= separator) 
                    w->val.push_back(make_pair(pos, max(max_term, 0.0)));
                else 
                    // w->val.push_back(make_pair(pos, max(value, 0.0)));
                    w->val.push_back(make_pair(pos, value));
            }
        }
    }
#ifdef GROUP_LASSO_DEBUG
    esmat_print (w, "[solved w]");
#endif
    esmat_trim (w);
#ifdef GROUP_LASSO_DEBUG
    esmat_print (w, "[w after trimming]");
#endif
    // STEP THREE: recollect temporary variable - wbar
    esmat_free (wbar);
}

void cvx_clustering (Esmat* dist_mat, int fw_max_iter, int D, int N, double lambda, double** W, int ADMM_max_iter, int SS_PERIOD) {
    // parameters 
    double alpha = 0.1;
    double rho = 1;
    ofstream ss_out ("plot_cputime_objective");
    ss_out << "Time Objective" << endl;
    clock_t cputime = 0;
    clock_t prev = clock();
    // iterative optimization 
    double error = INF;
    Esmat* wone = esmat_init (N, N);
    Esmat* wtwo = esmat_init (N, N);
    Esmat* yone = esmat_init (N, N);
    Esmat* ytwo = esmat_init (N, N);
    Esmat* z = esmat_init (N, N);
    Esmat* z_old = esmat_init (N, N);
    Esmat* diffzero = esmat_init (N, N);
    Esmat* diffone = esmat_init (N, N);
    Esmat* difftwo = esmat_init (N, N);
    esmat_zeros (wone);
    esmat_zeros (wtwo);
    esmat_zeros (yone);
    esmat_zeros (ytwo);
    esmat_zeros (z);
    esmat_zeros (z_old);
    esmat_zeros (diffzero);
    esmat_zeros (diffone);
    esmat_zeros (difftwo);

    // variables for shriking method
    set<int> col_active_sets;
    // set initial active_set as all elements
    for (int i = 0; i < N; i++) {
        col_active_sets.insert(i);
    }

    cputime += clock() - prev;
    ss_out << cputime << " " << 0 << endl;
    prev = clock();
    int iter = 0; // Ian: usually we count up (instead of count down)
    bool no_active_element = false, admm_opt_reached = false;
    while ( iter < ADMM_max_iter ) { // stopping criteria

        // STEP ONE: resolve w_1 and w_2
        frank_wolfe_solver (dist_mat, yone, z, wone, rho, N, fw_max_iter, col_active_sets);
        group_lasso_solver (ytwo, z, wtwo, rho, lambda);

#ifdef SUBPROBLEM_DUMP
     cout << "[Frank_wolfe]";
     subproblem_objective (dist_mat, yone, zone, wone, rho, N, lambda);
     cout << "[Blockwise]";
     subproblem_objective (dist_mat, ytwo, ztwo, wtwo, rho, N, lambda);
#endif

        // STEP TWO: update z by averaging w_1 and w_2
        // STEP THREE: update the y_1 and y_2 by w_1, w_2 and z
        set<int>::iterator it;
        Esmat* temp = esmat_init (N, N);
        Esmat* diff = esmat_init (N, N);
        esmat_add (wone, wtwo, temp);
        esmat_scalar_mult (0.5, temp, z);

        esmat_copy (yone, temp)
        esmat_sub (wone, z, diff);
        esmat_scalar_mult (alpha, diff);
        esmat_add (temp, diff, yone);

        esmat_copy (ytwo, temp)
        esmat_sub (wtwo, z, diff);
        esmat_scalar_mult (alpha, diff);
        esmat_add (temp, diff, ytwo);
        esmat_free (diff);
        esmat_free (temp);

        // STEP FOUR: trace the objective function
        if (iter < 3 * SS_PERIOD || (iter+1) % SS_PERIOD == 0) {
            cputime += clock() - prev;
            error = overall_objective (dist_mat, lambda, N, z);
            cout << "[Overall] iter = " << iter 
                << ", Loss Error: " << error << endl;
            ss_out << cputime << " " << error << endl;
            prev = clock();
        }
        // Shrinking Method:
        // STEP ONE: reduce number of elements considered in next iteration
        /*
        esmat_trim (z, ADMM_EPS);
        vector<int> col_to_shrink;
        // TODO: detect active set O(N*K*eps)
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
        // TODO: modify
        int num_col_to_shrink = col_to_shrink.size();
        for (int s = 0; s < num_col_to_shrink; s ++) {
            int j_shr = col_to_shrink[s];
            col_active_sets.erase(j_shr);
            for(int i=0;i<N;i++){
                wone[i][j_shr] = 0.0;
                z[i][j_shr] = 0.0;
                z_old[i][j_shr] = 0.0;
            }
        }
        // update z_old
        esmat_copy (z, z_old);
        // count number of active elements
        int num_active_cols = col_active_sets.size();
        if ((iter+1) % SS_PERIOD == 0) {
            cout << "iter: " << iter;
            cout << ", num_active_cols: " << num_active_cols <<endl;
        }
        int num_active_elements=N*num_active_cols;
        */
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
    esmat_free (wone);
    esmat_free (wtwo);
    esmat_free (yone);
    esmat_free (ytwo);
    esmat_free (diffone);
    esmat_free (difftwo);
    esmat_free (z_old);
    // STEP SIX: put converged solution to destination W
    esmat_copy (z, W);
    esmat_free (z);
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
    double lambda = atof(argv[4]);
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

    // pre-compute distance matrix
    dist_func df = L2norm;
    double** dist_mat = esmat_init (N, N);
    //  Esmat* dist_mat = esmat_read (dmatFile, N, N);
    mat_zeros (dist_mat, N, N);
    compute_dist_mat (data, dist_mat, N, D, df, true); 

    // Run sparse convex clustering
    Esmat* esmatW = esmat_init (N, N);  esmat_zeros (esmatW);
    cvx_clustering (dist_mat, fw_max_iter, D, N, lambda, esmatW, ADMM_max_iter, screenshot_period);

    double** W = esmat2mat (esmatW);
    /* Output cluster */
    output_objective(clustering_objective (dist_mat, W, N));
    /* Output cluster centroids */
    output_model (W, N);
    /* Output assignment */
    output_assignment (W, data, N);

    /* reallocation */
    esmat_free (dist_mat, N, N);
    esmat_free (esmatW);
    mat_free (W, N, N);
}
