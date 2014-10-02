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
#define INF_INT 60000
/* algorithmic options */ 
 #define EXACT_LINE_SEARCH  // comment this to use inexact search

/* dumping options */
// #define FRANK_WOLFE_DUMP
// #define EXACT_LINE_SEARCH_DUMP
// #define BLOCKWISE_DUMP
// #define SUBPROBLEM_DUMP 

void frank_wolfe_solver (double** dist_mat, Esmat* yone, Esmat* zone, Esmat* wone, double rho, int N, int K, vector<int> col_active_map) {
    // cout << "within frank_wolfe_solver" << endl;
    // STEP ONE: compute gradient mat initially
    vector< vector<cell> > actives (N, vector<cell>());
    vector< priority_queue<cell, vector<cell>, CellCompare> > pqueues (N, priority_queue<cell,vector<cell>, CellCompare> ());

    int count = 0;
    for (int i = 0; i < N; i ++) {
        // cout << col_active_map[i] << endl;
        if (col_active_map[i] < 0) continue;
        else count ++;
    }
    // cout << "active: " << count << endl;
    
    int num_active_cols = col_active_map.size();
    double** grad = mat_init (N, num_active_cols);
    mat_zeros (grad, N, num_active_cols);
    // compute N by K active gradient matrix TODO: put this copy to overall ADMM
    for (int j = 0; j < N; j++) {
        int gj = col_active_map[j];
        if (gj < 0) continue;
        for (int i = 0; i < N; i ++)
            grad[i][gj] = 0.5*dist_mat[i][j];
    }
    // cout << "[grad]" << endl;
    // cout << mat_toString(grad, N, num_active_cols);
    Esmat* assist = esmat_init (N, N);
    Esmat* temp = esmat_init (N, N);
    esmat_sub (wone, zone, temp);
    esmat_scalar_mult (rho, temp);
    esmat_add (temp, yone, assist);
    esmat_free (temp);
    int assist_size = assist->val.size();
    for (int i = 0; i < assist_size; i++) {
        int assist_esmat_index = assist->val[i].first;
        int col_index = assist_esmat_index / assist->nRows;
        int grad_col_index = col_active_map[col_index];
        if (grad_col_index < 0) continue;
        int row_index = assist_esmat_index % assist->nRows;
        double value = assist->val[i].second;
        grad[row_index][grad_col_index] += value;
    }
    esmat_free (assist);
    // set up actives and pqueues
    int grad_esmat_index = 0;
    int yone_index = -1, yone_esmat_index = -1;
    int zone_index = -1, zone_esmat_index = -1;
    int wone_index = -1, wone_esmat_index = -1;
    int w_size = wone->val.size();
    int y_size = yone->val.size();
    int z_size = zone->val.size();
    int grad_size = N* num_active_cols;
    if (w_size == 0) wone_esmat_index = INF_INT;
    else {
        while (true) {
            ++ wone_index;
            if (wone_index == w_size) wone_esmat_index = INF_INT;
            else {
                wone_esmat_index = wone->val[wone_index].first;
                int wone_col_index = col_active_map[wone_esmat_index/wone->nRows];
                if (wone_col_index < 0) continue;
                wone_esmat_index = wone_col_index*wone->nRows + wone_esmat_index % wone->nRows;
            }
            break;
        }
    }
    if (y_size == 0) yone_esmat_index = INF_INT;
    if (z_size == 0) zone_esmat_index = INF_INT;
    // TODO: conversion from grad_esmat_index to NN esmat index
    while (grad_esmat_index < grad_size) {
        int row_index = grad_esmat_index % N;
        int col_index = grad_esmat_index / N;
        double g = grad[row_index][col_index];
        double w,y,z;
        // get y
        while (yone_esmat_index < grad_esmat_index) {
            ++ yone_index;
            if (yone_index == y_size) yone_esmat_index = INF_INT;
            else {
                yone_esmat_index = yone->val[yone_index].first;
                int yone_col_index = col_active_map[yone_esmat_index/yone->nRows];
                if (yone_col_index < 0) continue;
                yone_esmat_index = yone_col_index*yone->nRows + yone_esmat_index % yone->nRows;
            }
        }
        if (yone_esmat_index > grad_esmat_index) y = 0.0;
        else if (yone_esmat_index == grad_esmat_index) 
            y = yone->val[yone_index].second;
        // get z 
        while (zone_esmat_index < grad_esmat_index) {
            ++ zone_index;
            if (zone_index == z_size) zone_esmat_index = INF_INT;
            else {
                zone_esmat_index = zone->val[zone_index].first;
                int zone_col_index = col_active_map[zone_esmat_index/zone->nRows];
                if (zone_col_index < 0) continue;
                zone_esmat_index = zone_col_index*zone->nRows + zone_esmat_index % zone->nRows;
            }
        }
        if (zone_esmat_index > grad_esmat_index) z = 0.0;
        else if (zone_esmat_index == grad_esmat_index) 
            z = zone->val[zone_index].second;
        // get w and insert to set or push to stack
        if (grad_esmat_index == wone_esmat_index) {
            w = wone->val[wone_index].second; 
            actives[row_index].push_back(cell(col_index,w,z,y,g)); // active
            while (true) {
            ++ wone_index;
            if (wone_index == w_size) wone_esmat_index = INF_INT;
            else {
                wone_esmat_index = wone->val[wone_index].first;
                int wone_col_index = col_active_map[wone_esmat_index/wone->nRows];
                if (wone_col_index < 0) continue;
                wone_esmat_index = wone_col_index*wone->nRows + wone_esmat_index % wone->nRows;
            }
            break;
            }
        } else if (grad_esmat_index < wone_esmat_index) { // TODO: this is the problem
            w = 0.0;
            pqueues[row_index].push(cell(col_index,w,z,y,g)); // potentially active
        }
        ++ grad_esmat_index;
    }
    /*
    for (int i = 0;i < N; i ++) {
        cout << "size: " << pqueues[i].size() << endl;
    }
    */
    // STEP TWO: iteration solve each row 
    int k = 0;  // iteration number
    vector<bool> is_fw_opt_reached (N, false);
    // cout << "within frank_wolfe_solver: start iteration" << endl;
    while (k < K) { 
        // cout << "k: " << k << endl;
        // compute new active atom: can be in active set or not
        vector<cell *> s (N, NULL);
        vector<bool> isInActives (N, false);
        for (int i = 0; i < N; i++) {
            if (is_fw_opt_reached[i]) continue;
            if (pqueues[i].size() <= 0) continue;
            cell tmp_cell;
            tmp_cell.copy_from(pqueues[i].top());
            s[i] = &tmp_cell;
            // cout << s[i].first << ":" << s[i].second << endl;
            vector<cell>::iterator it;
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) {
                // take the minimal of each row
                if (it->grad < s[i]->grad) {
                    isInActives[i] = true;
                    s[i] = &(*it);
                    // s[i].copy_from(*it);
                }
            }
            /*
            cout << "s[" << i <<  "]" 
                 << "j="<< s[i]->index
                << ",w=" << s[i]->w 
                << ",y=" << s[i]->y 
                << ",z=" << s[i]->z 
                << ",grad=" << s[i]->grad << endl;
                */
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
                double w_minus_z = it->w - it->z;
                if (it->index == s[i]->index) {
                    w_minus_s = it->w -1.0;
                } else {
                    w_minus_s = it->w;
                }
                sum1 += 0.5 * w_minus_s * (dist_mat[i][it->index] -r);
                sum2 += it->y * w_minus_s;
                sum3 += rho * w_minus_s * w_minus_z;
                sum4 += rho * w_minus_s * w_minus_s; 
            }
            if (!isInActives[i]) {
                sum1 += 0.5 * (-1.0) * (dist_mat[i][s[i]->index] - r);
                sum2 += s[i]->y * (-1.0);
                sum3 += rho * (-1.0) * (s[i]->w - s[i]->z);
                sum4 += rho;
            }

            if (fabs(sum4) > 0) {
                gamma = (sum1 + sum2 + sum3) / sum4;
                gamma = max(gamma, 0.0);
                gamma = min(gamma, 1.0);
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
            } else {
                gamma = 0.0;
                is_fw_opt_reached[i] = true;
            }
#else
            gamma = 2.0 / (k+2.0);
#endif
            // update wone
            int active_size = actives[i].size();
            for (int actj = 0; actj < active_size; actj++) 
                actives[i][actj].w *= (1-gamma);
            s[i]->w += gamma; // gamma * (1.0)
            // update new actives 
            if (!isInActives[i]) {
                actives[i].push_back(*(s[i]));
                pqueues[i].pop();
                ++ active_size;
            }
            // iterate through all active elements
            for (int actj = 0; actj < active_size; actj++) {
                int j = actives[i][actj].index;
                cell* tc = &(actives[i][actj]);
                tc->grad = 0.5*dist_mat[i][j] + tc->y + rho*(tc->w - tc->z); 
            }
            // cout << "actives[" << i << "]: " << actives[i].size() << endl;
        }
        // cout << "within frank_wolfe_solver: next iteration" << endl;
        k ++;
    }
    // reconstruct esmat wone
    esmat_zeros(wone);
    for (int i = 0; i < N; i++) {
        int num_active_elem = actives[i].size();
        for (int j = 0; j < num_active_elem; j++) {
            int esmat_index = i + actives[i][j].index * N;
            /*
            cout << "esmat_index: " << esmat_index
                << ", r="<< esmat_index % N
                << ", c="<< esmat_index / N  << endl;
                */
            double value = actives[i][j].w;
            if (value > ADMM_EPS)
                wone->val.push_back(make_pair(esmat_index, value));
        }
    }
    esmat_align (wone); // Elog(E)
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

void cvx_clustering (double** dist_mat, int fw_max_iter, int D, int N, double lambda, Esmat* W, int ADMM_max_iter, int SS_PERIOD) {
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
    vector<int> col_active_map (N, -1);
    for (int i = 0; i < N; i ++) 
        col_active_map[i] = i;

    cputime += clock() - prev;
    ss_out << cputime << " " << 0 << endl;
    prev = clock();
    int iter = 0; // Ian: usually we count up (instead of count down)
    bool no_active_element = false, admm_opt_reached = false;
    while ( iter < ADMM_max_iter ) { // stopping criteria

        // STEP ONE: resolve w_1 and w_2
        frank_wolfe_solver (dist_mat, yone, z, wone, rho, N, fw_max_iter, col_active_map);
         cout << "[wone]" << endl;
         cout << esmat_toString(wone);
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
        
        cout << "[z]" << endl;
        cout << esmat_toString (z) << endl;

        esmat_copy (yone, temp);
        esmat_sub (wone, z, diff);
        esmat_scalar_mult (alpha, diff);
        esmat_add (temp, diff, yone);

        esmat_copy (ytwo, temp);
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
        esmat_trim (z, ADMM_EPS);
        esmat_trim (wtwo, ADMM_EPS);
        Esmat* temp1 = esmat_init (wone);
        Esmat* temp2 = esmat_init (z);
        esmat_add (wone, wtwo, temp1); 
        esmat_add (z, z_old, temp2); 
        Esmat* temp3 = esmat_init (temp1);
        esmat_add (temp1,temp2,temp3);
        esmat_free (temp2);
        esmat_count_over_col (temp3, temp1);
        int temp_size = temp1->val.size();
        int map_index = 0;
        assert (temp1->nRows == 1);
        // TODO: can improve
        for (int i = 0; i < N; i ++) 
            col_active_map[i] = -1;
        for (int i = 0; i < temp_size; i ++) {
            int temp_index = temp1->val[i].first;
            col_active_map[temp_index] = map_index;
            ++ map_index;
        }
        esmat_free (temp1);
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
    cerr << "lambda = " << lambda << endl;
    cerr << "r = " << r << endl;
    cerr << "Screenshot period = " << screenshot_period << endl;
    int seed = time(NULL);
    srand (seed);
    cerr << "seed = " << seed << endl;

    // pre-compute distance matrix
    dist_func df = L2norm;
    double** dist_mat = mat_init (N, N);
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
    mat_free (dist_mat, N, N);
    esmat_free (esmatW);
    mat_free (W, N, N);
}
