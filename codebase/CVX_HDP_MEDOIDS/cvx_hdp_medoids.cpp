/*###############################################################
## MODULE: cvx_hdp_medoids.cpp
## VERSION: 2.0 
## SINCE 2014-07-21
## AUTHOR:
##     Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
##
## DESCRIPTION: 
##     Convex relaxation for cvx_hdp_medoids resolution
################################################################# 
## Edited by MacVim
## Class Info auto-generated by Snippet 
################################################################*/

/* TODO list:
  [DONE] 1. compute_dist_mat generates matrix with N by D -- checked
  [DONE] 2. fixup bugs in frank_wolfe_solver -- checked
  [DONE] 3. fixup bugs in exact line search  -- checked
  [DONE] 4. fixup problem in group_lasso_solver
  5. fixup local problem separation and -1 problem
  [DONE] 6. develop early stopping detection for frank_wolfe_solver
  [DONE] 7. refactor codes and refine the exSparseMat API
*/

/*
Shrinking methods for frank_wolfe_solver:
  1. inherently inactive for n,k that dist_mat[n][k] == inf
  2. inactive 
  3. gamma < TRIM_THRESHOLD  
 */

#include "cvx_hdp_medoids.h"
using namespace std;

/* algorithmic options */ 
#define EXACT_LINE_SEARCH  // comment this to use inexact search

/* dumping options */
// #define FRANK_WOLFE_DEBUG
// #define EXACT_LINE_SEARCH_DUMP
// #define COVERAGE_SUBPROBLEM_DUMP
// #define LOCAL_SUBPROBLEM_DUMP
// #define GROUP_LASSO_DEBUG
// #define BLOCKWISE_DUMP
// #define DIST_MAT_DUMP
#define SUBPROBLEM_DUMP
#define EPSILON 1e-3
bool iszero (double value) { return (value > -EPSILON && value < EPSILON); }
bool isone (double value) { return (value > 1-EPSILON && value < 1+EPSILON); }

bool stopping (Esmat* z) {
    bool tostop = true;
    int SIZE = z->val.size();
    for (int i = 0; i < SIZE; i ++) {
        double value = z->val[i].second;
        if (!iszero(value) && !isone(value)) {
            tostop = false;
        } 
    }
    return tostop;
}
/* Compute the mutual distance of input instances contained within "data" */
void compute_dist_mat (Esmat* dist_mat, Lookups* tables, int N, int D) {
    // STEP ZERO: parse input
    vector< pair<int,int> > doc_lookup = *(tables->doc_lookup);
    vector< pair<int,int> > word_lookup = *(tables->word_lookup); 

    // STEP ONE: compute distribution for each document
    vector< map<int, double> > distributions (D, map<int, double>());
    for (int d = 0; d < D; d ++) {
        // a. compute sum of word frequency
        int sumFreq = 0;
        for (int w = doc_lookup[d].first; w < doc_lookup[d].second; w++) {
            sumFreq += word_lookup[w].second;
        }
        // b. compute distribution
        for (int w = doc_lookup[d].first; w < doc_lookup[d].second; w++) {
            int voc_index = word_lookup[w].first;
            double prob = 1.0 * word_lookup[w].second / sumFreq;
            distributions[d].insert(pair<int, double> (voc_index, prob));
        }
    }
    // STEP TWO: compute weight of word within one document
    esmat_zeros(dist_mat);
    for (int j = 0; j < D; j ++) {
        for (int d = 0; d < D; d ++) {
            for (int w = doc_lookup[d].first; w < doc_lookup[d].second; w++) {
                int voc_index = word_lookup[w].first;
                int count_w_d1 = word_lookup[w].second;    
                double prob_w_d2;
                double dist;
                map<int, double>::const_iterator iter;
                iter = distributions[j].find(voc_index);
                if (iter == distributions[j].end()) {
                    prob_w_d2 = 0.0;
                    dist = INF;
                } else {
                    prob_w_d2 = iter->second;
                //   dist(w, d2) =  - count_w_d1 * log( prob_d2(w) )
                    dist = - count_w_d1 * log(prob_w_d2);
                }
                int esmat_index = w + N * j;
#ifdef DIST_MAT_DUMP
                cout << "index: " << esmat_index << ", dist: " << dist << ", ";
                cout << "count_w_d1: " << count_w_d1 << ", prob_w_d2: " << prob_w_d2 << endl;
#endif
                dist_mat->val.push_back(make_pair(esmat_index, dist));
            }
        }
    }
}
void get_all_centroids(Esmat* W, vector<int>* centroids) {
    Esmat* max_belonging = esmat_init();
    esmat_max_over_col (W, max_belonging);
    int size = max_belonging->val.size();
    for (int i = 0; i < size; i ++ ) {
        if (fabs(max_belonging->val[i].second) > 0.3) {
            centroids->push_back(max_belonging->val[i].first +1);
        }
    }
    esmat_free(max_belonging);
}
void output_model (Esmat* W) {
    vector<int> centroids;
    get_all_centroids (W, &centroids); // contains index of all centroids
    int nCentroids = centroids.size();
    ofstream model_out ("opt_model");
    model_out << "nCentroids: " << nCentroids << endl;
    for (int i = 0; i < nCentroids; i ++) {
        model_out << "centroids[" << i <<"]: " << centroids[i] << endl;
    }
    model_out.close();
}
void output_assignment (Esmat* W, Lookups * tables) {
    int D = tables->nDocs;
    int N = tables->nWords;
    vector< pair<int,int> >* word_lookup = tables->word_lookup;
    vector< pair<int,int> > doc_lookup = *(tables->doc_lookup);
    vector< vector< pair<int, double> > > words_asgn (N, vector< pair<int, double> > ());
    int nEntries = W->val.size();
    for (int i = 0; i < nEntries; i ++) {
        int esmat_index = W->val[i].first;
        int row_index = esmat_index % W->nRows;
        int col_index = esmat_index / W->nRows;
        double value = W->val[i].second;
        words_asgn[row_index].push_back(make_pair(col_index, value));
    }
    ofstream asgn_out ("opt_assignment");
    // get all centroids
    for (int d = 0; d < D; d++) {
        asgn_out << "d = " << d << endl;
        for (int i = doc_lookup[d].first; i < doc_lookup[d].second; i ++) {
            // output identification and its belonging
            asgn_out << "  id=" << i+1 << ", voc_index=" << (*word_lookup)[i].first << ", "; 
            for (int j = 0; j < words_asgn[i].size(); j ++) {
                int index = words_asgn[i][j].first;
                double value = words_asgn[i][j].second;
                if( fabs(value) > 0.0003 ) {
                    asgn_out << index+1 << "(" << value << ")";
                }
            }
            asgn_out << endl;
        }
    }
    asgn_out.close();
}
/* dummy_penalty = r dot (1 - sum_k w_nk) */
double get_dummy_loss (Esmat* Z) {
    Esmat* temp_vec = esmat_init (Z->nRows, 1);
    esmat_sum_row (Z, temp_vec);
    double dummy_loss = esmat_compute_dummy (temp_vec);
    esmat_free (temp_vec);
    return dummy_loss;
}
/* \lambda_g \sumk \maxn |\wnk| */
double get_global_topic_reg (Esmat* absZ, double lambda) {
    if (absZ->val.size() == 0) {
        return 0.0;
    }
    Esmat* maxn = esmat_init (1, absZ->nCols);
    Esmat* sumk = esmat_init (1, 1);
    esmat_max_over_col (absZ, maxn);
    esmat_sum_row (maxn, sumk);
    double global_topic_reg = -INF;
    if (sumk->val.size() > 0)
        global_topic_reg = lambda * sumk->val[0].second; 
    else 
        global_topic_reg = 0.0;
    esmat_free (sumk);
    esmat_free (maxn);
    return global_topic_reg;
}
/* \lambdal \sum_d \sum_k \underset{n \in d}{\text{max}} |\wnk| */
double get_local_topic_reg (Esmat* absZ, double lambda, vector< pair<int,int> >* doc_lookup) {
    // STEP ONE: initialize sub matrix for each document
    int nDocs = doc_lookup->size();
    vector<Esmat*> sub_absZ (nDocs);
    for (int d = 0; d < nDocs; d ++) {
        sub_absZ[d] = esmat_init (0,0);
    }

    // STEP TWO: separate entire Z to submat Z
    esmat_submat_row (absZ, sub_absZ, doc_lookup);

    // STEP THREE: compute global topic regularizer for each localized doc
    double local_topic_reg = 0.0;
    for (int d = 0; d < nDocs; d ++) {
        // cout << "size[" << d << "]: " << sub_absZ[d]->val.size() << endl;
        local_topic_reg += get_global_topic_reg (sub_absZ[d], lambda); 
    }

    // Final: free resource
    esmat_free_all (sub_absZ);
    return local_topic_reg;
}
double subproblem_objective (int prob_index, Esmat* Y, Esmat* Z, Esmat* W, double RHO, double lambda, Lookups* tables, Esmat* dist_mat) {
    string title = "";
    vector< pair<int,int> >* doc_lookup = tables->doc_lookup;
    vector< pair<int,int> >* word_lookup = tables->word_lookup; 
    vector< vector<int> >* voc_lookup = tables->voc_lookup;
    /*
       esmat_print (W, "W["+ to_string(prob_index) + "]")
       esmat_print (Y, "Y["+ to_string(prob_index) + "]")
       esmat_print (Z, "Z["+ to_string(prob_index) + "]")
       */

    // STEP ONE: compute main term
    double total = 0.0;
    double main = -1.0;
    if (prob_index == 1) {
        double loss = esmat_frob_prod (dist_mat, Z);
        total += loss;
#ifdef SUBPROBLEM_DUMP
        cout << "loss: " << loss << ", ";
#endif
        title = "dummy";
        // cout << "begin to compute dummy loss" << endl;
        // dummy_penalty = r dot (1 - sum_k w_nk)
        main = get_dummy_loss (W);
    } else {
        Esmat* absW = esmat_init (W);
        esmat_abs (W, absW);
        if (prob_index == 2) {
            title = "Global_Reg";
            main = get_global_topic_reg (absW, lambda);
        } else if (prob_index == 3) {
            title = "Local_Reg";
            main = get_local_topic_reg (absW, lambda, doc_lookup);
        }
        esmat_free (absW);
    }
    Esmat* w_minus_z = esmat_init ();
    // STEP TWO: compute linear term: linear = y^T dot (w - z) 
    esmat_sub (W, Z, w_minus_z); // temp = w - z
    double linear = esmat_frob_prod (Y, w_minus_z);
    // STEP THREE: compute quadratic term: quadratic = 0.5 * RHO * || w - z ||^2 
    double quadratic = 0.5 * RHO * esmat_frob_norm (w_minus_z);
    total += main + linear + quadratic;
    esmat_free (w_minus_z);
#ifdef SUBPROBLEM_DUMP
    cout << title << ": " << main << ", ";
    cout << "linear: " << linear << ", ";
    cout << "quadratic: " << quadratic << ", ";
    cout << "total: " << total << endl;
#endif
    return total;
}
double original_objective (Esmat* Z, vector<double> LAMBDAs, Lookups* tables) {
    vector< pair<int,int> >* doc_lookup = tables->doc_lookup;
    vector< pair<int,int> >* word_lookup = tables->word_lookup; 
    vector< vector<int> >* voc_lookup = tables->voc_lookup;

    Esmat* absZ = esmat_init (Z);
    esmat_abs (Z, absZ);
    // STEP ONE: compute dummy loss
    double dummy = get_dummy_loss (Z);
    // cout << "dummy =" << dummy << endl;

    // STEP TWO: compute "GLOBAL TOPIC" group-lasso regularization
    double global_topic_reg = get_global_topic_reg (absZ, LAMBDAs[0]);
    esmat_free (absZ); 
    return dummy + global_topic_reg;
}
void frank_wolfe_solver (Esmat* dist_mat, Esmat * Y_1, Esmat * Z_1, Esmat * w_1, double RHO) {

    bool is_global_optimal_reached = false;

    // cout << "within frank_wolfe_solver" << endl;
    // This can be computed by using corner point. 
    Esmat * gradient = esmat_init (w_1);
    Esmat * s = esmat_init (w_1);
    Esmat * old_s = esmat_init (w_1);

#ifndef EXACT_LINE_SEARCH
    int K = 300;
#else
    int K = 50;
    Esmat * w_minus_s = esmat_init ();
    Esmat * w_minus_z = esmat_init ();
#endif

    int k = 0; // iteration number
    double gamma; // step size
    double penalty;
    Esmat * tempS = esmat_init (w_1);

    //esmat_zeros (w_1);
    double sum1, sum2, sum3, sum4;
    // cout << "within frank_wolfe_solver: start iteration" << endl;
    while (k < K && !is_global_optimal_reached) {
        // STEP ONE: find s that minimizes <s, grad f>
        // gradient[i][j] = dist_mat[i][j] +  Y_1[i][j] + RHO * (w_1[i][j] - Z_1[i][j]) ;  //- r[i];
        esmat_sub (w_1, Z_1, tempS);
        esmat_scalar_mult (RHO, tempS);
        esmat_add (Y_1, tempS, gradient);
        esmat_copy (gradient, tempS);
        esmat_add (tempS, dist_mat, gradient);
        esmat_min_row (gradient, s);

#ifdef FRANK_WOLFE_DEBUG
        cout << "=============================" << endl;
        esmat_print(w_1, "[w_1]");
        esmat_print(gradient, "[gradient]");
        esmat_print(s, "[s]");
        esmat_print(old_s, "[old_s]");
        cout << "esmat_equal(s, old_s): " << esmat_equal(s, old_s) << endl;
#endif 

#ifdef EXACT_LINE_SEARCH
        // early stopping strategy
        if (esmat_equal(s, old_s)) {
            is_global_optimal_reached = true;
            break;
        } else {
            esmat_copy (s, old_s);
        }
#endif

        // STEP TWO: apply exact or inexact line search to find solution
#ifndef EXACT_LINE_SEARCH
        // Here we use inexact line search
        gamma = 2.0 / (k + 2.0);
#else
        // Here we use exact line search 
        if (k == 0) {
            gamma = 1.0;
        } else {
            // gamma* = (sum1 + sum2 + sum3) / sum4, where
            // sum1 = 1/2 sum_n sum_k (w - s)_nk * || x_n - mu_k ||^2
            // sum2 = sum_n sum_k (w - s)_nk
            // sum3 = - RHO * sum_n sum_k  (w - z) 
            // sum4 = sum_n sum_k RHO * (s - w)
            esmat_sub (w_1, s, w_minus_s);
            // cout << esmat_toString(w_minus_s);
            esmat_sub (w_1, Z_1, w_minus_z);

#ifdef EXACT_LINE_SEARCH_DUMP
            cout << "[w_minus_s]" << endl;
            cout << esmat_toString(w_minus_s);
            cout << "[w_minus_z]" << endl;
            cout << esmat_toString(w_minus_z);
#endif
            // NOTE: in case of ||w_1 - s||^2 = 0, not need to optimize anymore
            // since incremental term = w + gamma (s - w), and whatever gamma is,
            // w^(k+1) = w^(k), this would be equivalent to gamma = 0
            if (esmat_frob_norm (w_minus_s) == 0) {
                gamma = 0;
                is_global_optimal_reached = true;
                // reach the exit condition, do not make more iteration
            } else {
                sum1 = 0.5 * esmat_frob_prod (w_minus_s, dist_mat);
                sum2 = esmat_frob_prod (Y_1, w_minus_s);
                sum3 = RHO * esmat_frob_prod (w_minus_z, w_minus_s);
                sum4 = RHO * esmat_frob_norm (w_minus_s);
                // gamma should be within interval [0,1]
                // gamma = (sum1 + sum2 + sum3) / sum4;
                gamma = (sum1 + sum2 + sum3) / sum4;

#ifdef EXACT_LINE_SEARCH_DUMP
                cout << "[exact line search] (sum1, sum2, sum3, sum4, gamma) = ("
                    << sum1 << ", " << sum2 << ", " << sum3 << ", " << sum4 << ", " << gamma
                    << ")"
                    << endl;
#endif
            }
            /* postcondition for early stopping */
            if (gamma < TRIM_THRESHOLD * 1e-3) {
                is_global_optimal_reached = true;
            }
        }
#endif
        // update the w^(k+1)
        esmat_scalar_mult (gamma, s); 
        // esmat_print (s, "[gamma*s]");
        esmat_scalar_mult (1.0-gamma, w_1);
        esmat_copy (w_1, tempS);
        esmat_add (tempS, s, w_1);
        // esmat_print (w_1, "[post w_1]");

        // cout << "within frank_wolfe_solver: step three finished" << endl;
        // report the #iter and objective function
        /*
           cout << "[Frank-Wolfe] iteration: " << k << ", first_subpro_obj: " << penalty << endl;
           */

        k ++;
        // cout << "within frank_wolfe_solver: next iteration" << endl;
    }
    // cout << "within frank_wolfe_solver: to free gradient" << endl;
    esmat_free (gradient);
    // cout << "within frank_wolfe_solver: to free temps" << endl;
    esmat_free (tempS);
    // cout << "within frank_wolfe_solver: to free s " << endl;
    esmat_free (s);
    // cout << "end frank_wolfe_solver: finished! " << endl;
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

/* three subproblems that employed group_lasso_solver in different ways */
void global_topic_subproblem (Esmat* Y, Esmat* Z, Esmat* w, double RHO, double lambda) {
    group_lasso_solver (Y, Z, w, RHO, lambda);
}
void local_topic_subproblem (Esmat* Y, Esmat* Z, Esmat* w, double RHO, double lambda, Lookups* tables) {
    vector< pair<int,int> >* doc_lookup = tables->doc_lookup;
    vector< pair<int,int> >* word_lookup = tables->word_lookup; 
    vector< vector<int> >* voc_lookup = tables->voc_lookup;

#ifdef LOCAL_SUBPROBLEM_DUMP
    cout << "[local w input]" << w->nRows << "," << w->nCols << "," << w->val.size() << endl;
    esmat_toString (w);
#endif

    int nDocs = tables->nDocs;
    Esmat* tempW = esmat_init (w);
    vector<Esmat*> subW (nDocs); 
    vector<Esmat*> subY (nDocs);
    vector<Esmat*> subZ (nDocs);
    // STEP ZERO: initialize all submats
    for (int d = 0; d < nDocs; d ++) {
        subW[d] = esmat_init (0,0);
        subY[d] = esmat_init (0,0);
        subZ[d] = esmat_init (0,0);
    }
    // STEP ONE: separate esmat Y, Z, w to multiple small-sized esmat
    // NOTE: all esmat position has to be recomputed
    esmat_submat_row (w, subW, doc_lookup);
    esmat_submat_row (Y, subY, doc_lookup);
    esmat_submat_row (Z, subZ, doc_lookup);

    for (int d = 0; d < nDocs; d ++) {
#ifdef LOCAL_SUBPROBLEM_DUMP
        cout << "subW[d" << d << "]" << esmat_toInfo(subW[d]);
        cout << "subY[d" << d << "]" << esmat_toInfo(subY[d]);
        cout << esmat_toString(subY[d]);
        cout << "subZ[d" << d << "]" << esmat_toInfo(subZ[d]);
        cout << esmat_toString(subZ[d]);
#endif
        // STEP TWO: invoke group_lasso_solver to each individual group
        group_lasso_solver (subY[d], subZ[d], subW[d], RHO, lambda);

#ifdef LOCAL_SUBPROBLEM_DUMP
        cout << "res_subW[d" << d << "]" << esmat_toInfo(subW[d]);
        cout << esmat_toString(subW[d]);
        cout << endl;
#endif

        // STEP TREE: merge solution of each individual group (place back)
        // NOTE: all esmat position has to be recomputed
        int start_row = (*doc_lookup)[d].first;
        int end_row = (*doc_lookup)[d].second;
        esmat_merge_row (subW[d], start_row, end_row, tempW);
    }
    // realign the mat->val with index-increasing order
    esmat_align (tempW);
    // STEP FIVE: free auxiliary resource
    for (int d = 0; d < nDocs; d ++) {
        esmat_free (subW[d]);
        esmat_free (subY[d]);
        esmat_free (subZ[d]);
    }
#ifdef LOCAL_SUBPROBLEM_DUMP
    cout << "[local w before output]" << w->nRows << "," << w->nCols << "," << w->val.size() << endl;
    cout << "[tempW] " << esmat_toInfo(tempW);
#endif

    // FINAL: update merged solution to w
    esmat_copy (tempW, w);

#ifdef LOCAL_SUBPROBLEM_DUMP
    cout << "[local w output]" << w->nRows << "," << w->nCols << "," << w->val.size() << endl;
#endif
}
void cvx_hdp_medoids (Esmat* dist_mat, vector<double> LAMBDAs, Esmat* W, Lookups* tables) {
    // SET MODEL-RELEVANT PARAMETERS 
    assert (LAMBDAs.size() == 2);
    double ALPHA = 0.1;
    double RHO = 1.0;
    int N = tables->nWords;
    int D = tables->nDocs;
    /* DECLARE AND INITIALIZE INVOLVED VARIABLES AND MATRICES */
    Esmat* w_1 = esmat_init (N, D);
    Esmat* w_2 = esmat_init (N, D);
    Esmat* w_3 = esmat_init (N, D);
    Esmat* y_1 = esmat_init (N, D);
    Esmat* y_2 = esmat_init (N, D);
    Esmat* y_3 = esmat_init (N, D);
    Esmat* z = esmat_init (N, D);
    Esmat* diff_1 = esmat_init (N, D);
    Esmat* diff_2 = esmat_init (N, D);
    Esmat* diff_3 = esmat_init (N, D);

    /* SET ITERATION-RELEVANT VARIABLES */
    double error = INF;
    int iter = 0; 
    // int max_iter = 5000;
    /* ITERATIVE OPTIMIZATION */
    // while ( iter < max_iter ) { // STOPPING CRITERIA
    while ( true ) { // STOPPING CRITERIA
        cout << "###################[iter:"<<iter<<"]#####################" << endl;
        // STEP ZERO: RESET ALL SUBPROBLEM SOLUTIONS (OPTIONAL) 
        esmat_zeros (w_1);
        esmat_zeros (w_2);
        esmat_zeros (w_3);

        // STEP ONE: RESOLVE W_1, W_2, W_3, W_4
        // resolve w_1
        // cout << "[w_1 before frank_wolfe_solver] " << esmat_toInfo(w_1);
        // cout << esmat_toString (w_1);
        frank_wolfe_solver (dist_mat, y_1, z, w_1, RHO); 
        double sub1_obj = subproblem_objective (1, y_1, z, w_1, RHO, 0.0, tables, dist_mat);
        // cout << "[w_1 after frank_wolfe_solver]" << esmat_toInfo(w_1);
        // out << esmat_toString (w_1);

        // resolve w_2
        global_topic_subproblem (y_2, z, w_2, RHO, LAMBDAs[0]);
        // compute value of objective function
        double sub2_obj = subproblem_objective (2, y_2, z, w_2, RHO, LAMBDAs[0],tables, dist_mat);
        // cout << "sub2_objective: " << sub2_obj << endl;
        // esmat_print(w_2, "[w_2] ");

        // cout << "[w_3]" << w_3->nRows << "," << w_3->nCols << "," << w_3->val.size() << endl;
        local_topic_subproblem (y_3, z, w_3, RHO, LAMBDAs[1], tables);
        // cout << "[w_3]" << w_3->nRows << "," << w_3->nCols << "," << w_3->val.size() << endl;
        double sub3_obj = subproblem_objective (3, y_3, z, w_3, RHO, LAMBDAs[1], tables, dist_mat);
        // esmat_print(w_3, "[w_3] ");

        // STEP TWO: update z by averaging w_1, w_2 and w_4
        Esmat* temp = esmat_init (N, D);
        esmat_add (w_1, w_2, z);
        esmat_copy (z, temp);
        esmat_add (temp, w_3, z);
        esmat_scalar_mult (1.0/3.0, z);

        // STEP THREE: update the y_1 and y_2 by w_1, w_2 and z
        esmat_sub (w_1, z, diff_1);
        esmat_scalar_mult (ALPHA, diff_1);
        esmat_copy (y_1, temp);
        esmat_add (temp, diff_1, y_1);

        esmat_sub (w_2, z, diff_2);
        esmat_scalar_mult (ALPHA, diff_2);
        esmat_copy (y_2, temp);
        esmat_add (temp, diff_2, y_2);

        esmat_sub (w_3, z, diff_3);
        esmat_scalar_mult (ALPHA, diff_3);
        esmat_copy (y_3, temp);
        esmat_add (temp, diff_3, y_3);

        // double trace_wone_minus_z = esmat_frob_norm (diff_1); 
        // double trace_wtwo_minus_z = esmat_frob_norm (diff_2); 
        // double trace_wtwo_minus_z = esmat_frob_norm (diff_3); 
        // double trace_wfour_minus_z = esmat_frob_norm (diff_4); 

        /*
        if (iter % 10 == 0) {
            esmat_print(z, "[z]");
            esmat_print(y_1, "[y_1]");
            esmat_print(y_2, "[y_2]");
            esmat_print(y_3, "[y_3]");
        }
        */

        // STEP FOUR: trace the objective function
        iter ++;
        cout << endl;
        esmat_print (z, "[z]");
        if (stopping(z)) break;
        // esmat_print(z, "[z] ");
        // if (iter == 2) return;
    }

    // STEP FIVE: memory recollection
    esmat_free (w_1);
    esmat_free (w_2);
    esmat_free (w_3);

    esmat_free (y_1);
    esmat_free (y_2);
    esmat_free (y_3);

    esmat_free (diff_1);
    esmat_free (diff_2);
    esmat_free (diff_3);

    // STEP SIX: put converged solution to destinated W
    esmat_copy (z, W);
    esmat_free (z);
}

// entry main function
int main (int argc, char ** argv) {

    // EXCEPTION control: illustrate the usage if get input of wrong format
    if (argc < 5) {
        cerr << "Usage: " << endl;
        cerr << "\tcvx_hdp_medoids [voc_dataFile] [doc_dataFile] [lambda_global] [lambda_local]" << endl;
        exit(-1);
    }

    // PARSE arguments
    string voc_file (argv[1]);
    string doc_file (argv[2]);
    vector<double> LAMBDAs (2, 0.0);
    LAMBDAs[0] = atof(argv[3]); // lambda_document
    LAMBDAs[1] = atof(argv[4]); // lambda_topic
    // int max_iter = atoi(argv[5]); // max_iter

    // preprocess the input dataset
    vector<string> voc_list;
    voc_list_read (voc_file, &voc_list);
    cerr << "vocs read done! " << endl;
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
    cerr << "docs read done" << endl;

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
    cerr << "lambda_coverage = " << LAMBDAs[1] << endl;
    cerr << "TRIM_THRESHOLD = " << TRIM_THRESHOLD << endl;
    cerr << "seed = " << seed << endl;
    cerr << "###########################################" << endl;

    // Run sparse convex clustering
    int N = lookup_tables.nWords;
    int D = lookup_tables.nDocs;
    Esmat* W = esmat_init (lookup_tables.nWords, lookup_tables.nDocs);
    Esmat* dist_mat = esmat_init (N, D);
    compute_dist_mat (dist_mat, &lookup_tables, N, D);

    ofstream dmat_out ("dist_mat");
    dmat_out << esmat_toInfo(dist_mat);
    dmat_out << esmat_toString(dist_mat);
    cerr << "dist_mat output finished.." << endl;
    // cvx_hdp_medoids (dist_mat, LAMBDAs, W, &lookup_tables, max_iter);
    cvx_hdp_medoids (dist_mat, LAMBDAs, W, &lookup_tables);

    /* Output objective */
    output_objective(clustering_objective (dist_mat, W));

    /* Output cluster centroids */
    output_model (W);

    /* Output assignment */
    output_assignment (W, &lookup_tables);

    /* reallocation */
    esmat_free (W);
    esmat_free (dist_mat);
}
