/*###############################################################
## MODULE: HDP.cpp
## VERSION: 1.0 
## SINCE 2014-07-21
## AUTHOR:
##     Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
##
## DESCRIPTION: 
##     Convex relaxation for HDP resolution
################################################################# 
## Edited by MacVim
## Class Info auto-generated by Snippet 
################################################################*/

#include "HDP.h"
#include "exSparseMat.h"
#include <cassert>

/* algorithmic options */ 
#define EXACT_LINE_SEARCH  // comment this to use inexact search

/* dumping options */
// #define FRANK_WOLFE_DUMP
// #define EXACT_LINE_SEARCH_DUMP
// #define BLOCKWISE_DUMP
// #define NTOPIC_DUMP


double sign (int input) {
    if (input > 0) return 1.0;
    else if ( input < 0 ) return -1.0;
    else return 0.0;
}

bool pairComparator (const std::pair<int, double>& firstElem, const std::pair<int, double>& secondElem) {
    // sort pairs by second element with *decreasing order*
    return firstElem.second > secondElem.second;
}

typedef double (* dist_func) (Instance*, Instance*, int); 

/*{{{*/
/*
int get_nCentroids (double ** W, int nRows, int nCols) {

    int nTopics = 0;

    double * sum_belonging = new double [nCols];

    mat_max_col (W, sum_belonging, nRows, nCols);

    for (int i = 0; i < nCols; i ++ ) {
        if (sum_belonging[i] > 3e-1) { 
            nCentroids ++;
        }
    }
    
    delete[] sum_belonging;

    return nTopics;
}

vector<int> get_all_centroids(double ** W, int nRows, int nCols) {

    std::vector<int> topics;

    double * sum_belonging = new double [nCols];
    for (int i = 0; i < nCols; i ++) {
        sum_belonging[i] = 0.0;
    }

    mat_max_col (W, sum_belonging, nRows, nCols);

    for (int i = 0; i < nCols; i ++ ) {
        if (sum_belonging[i] > 3e-1) {
            centroids.push_back(i);
        }
    }
    
    delete[] sum_belonging;

    return topics;
}
*/
/*}}}*/

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
    Esmat* maxn = esmat_init (1, absZ->nCols);
    Esmat* sumk = esmat_init (1, 1);
    esmat_max_col (absZ, maxn);
    esmat_sum_row (maxn, sumk);
    double global_topic_reg = lambda * sumk->val[0].second; 
    esmat_free (sumk);
    esmat_free (maxn);
    return global_topic_reg;
}

/* \lambda_l \sumn \maxk |\wnk| */
double get_local_topic_reg (Esmat* absZ, double lambda) {
    Esmat* maxk = esmat_init (Z->nRows, 1);
    Esmat* sumn = esmat_init (1, 1);
    esmat_max_row (absZ, maxk);
    esmat_sum_col (maxk, sumn);
    double local_topic_reg = LAMBDAs[1] * sumn->val[0].second; 
    esmat_free (sumn);
    esmat_free (maxk);
    return local_topic_reg;
}
/* \lambdab \sumk \sum_v \max |\wnk^{(v)}| */
double get_coverage_reg (Esmat* absZ, double lambda) {
    // TODO
    
}

double sub_objective (int prob_index, Esmat* Y, Esmat* Z, Esmat* W, double RHO) {
    // STEP ONE: compute main term
    double main = -1.0;
    if (prob_index == 1) {
        // dummy_penalty = r dot (1 - sum_k w_nk)
        main = get_dummy_loss (W);
    } else {
        esmat_abs (W, absW);
        if (prob_index == 2) {
            main = get_global_topic_reg (absW);
        } else if (prob_index == 3) {
            main = get_local_topic_reg (absW);
        } else if (prob_index == 4) {
            // TODO:

        }
    }
    Esmat* w_minus_z = esmat_init (W);
    // STEP TWO: compute linear term: linear = y^T dot (w - z) 
    esmat_sub (W, Z, w_minus_z); // temp = w - z
    double linear = esmat_fdot (Y, w_minus_z);
    // STEP THREE: compute quadratic term: quadratic = 0.5 * RHO * || w - z ||^2 
    double quadratic = 0.5 * RHO * esmat_fnorm (w_minus_z);
    esmat_free (w_minus_z);
    /*
#ifdef FRANK_WOLFE_DUMP
    cout << "[Frank_wolfe] (loss, linear, quadratic, dummy, total) = (" 
         << sum1 << ", " << sum2 << ", " << sum3 << ", " << dummy_penalty << ", " << total
         <<  ")" << endl;
#endif
*/
    return main + linear + quadratic;
}
double original_objective (Esmat* Z, vector<double> LAMBDAs) {

    Esmat* absZ = esmat_init (Z);
    esmat_abs (Z, absz);

    // STEP ONE: compute dummy loss
    double dummy = get_dummy_loss (Z);
    // cout << "dummy =" << dummy << endl;

    // STEP TWO: compute "GLOBAL TOPIC" group-lasso regularization
    double global_topic_reg = get_global_topic_reg (absZ, LAMBDAs[0]);
    // STEP THREE: compute "LOCAL TOPIC" group-lasso regularization
    double local_topic_reg = get_local_topic_reg (absZ, LAMBDAs[1]);
    // STEP FOUR: TODO compute "TOPIC COVERAGE" group-lasso regularization
    double coverage_reg = 0.0;

    esmat_free (absZ);    
    return dummy + global_topic_reg + local_topic_reg + coverage_reg;
}
void frank_wolfe_solver (Esmat * Y_1, Esmat * Z_1, Esmat * w_1, double RHO, int N) {

    bool is_global_optimal_reached = false;

    // cout << "within frank_wolfe_solver" << endl;
    // This can be computed by using corner point. 
    Esmat * gradient = esmat_init (w_1);
    Esmat * s = esmat_init (w_1);

#ifndef EXACT_LINE_SEARCH
    int K = 300;
#else
    int K = 2;
    Esmat * w_minus_s = esmat_init (w_1);
    Esmat * w_minus_z = esmat_init (w_1);
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
        // gradient[i][j] = Y_1[i][j] + RHO * (w_1[i][j] - Z_1[i][j]) ;  //- r[i];
        esmat_sub (w_1, Z_1, tempS);
        esmat_scalar_mult (RHO, tempS);
        esmat_add (Y_1, tempS, gradient);
        esmat_min_row (gradient, s);

        /*
#ifdef FRANK_WOLFE_DUMP
        cout << "esmat_norm2 (w_1): " <<  esmat_norm2 (w_1) << endl;
        cout << "esmat_norm2 (Y_1): " <<  esmat_norm2 (Y_1) << endl;
        cout << "esmat_norm2 (gradient): " <<  esmat_norm2 (gradient) << endl;
        cout << "esmat_sum (s): " <<  esmat_sum (s) << endl;
#endif
*/
        // cout << "within frank_wolfe_solver: step one finished" << endl;

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
        esmat_sub (w_1, Z_1, w_minus_z);

        // NOTE: in case of ||w_1 - s||^2 = 0, not need to optimize anymore
        // since incremental term = w + gamma (s - w), and whatever gamma is,
        // w^(k+1) = w^(k), this would be equivalent to gamma = 0
        if (esmat_fnorm(w_minus_s) == 0) {
            gamma = 0;
            is_global_optimal_reached = true;
            // reach the exit condition, do not make more iteration
        } else {
            /*
            esmat_fdot (w_minus_s, dist_mat, tempS);
            sum1 = 0.5 * esmat_sum (tempS);
            */
            sum2 = esmat_fdot (Y_1, w_minus_s);
            sum3 = RHO * esmat_fdot (w_minus_z, w_minus_s);
            sum4 = RHO * esmat_fnorm (w_minus_s);
            // gamma should be within interval [0,1]
            gamma = (sum1 + sum2 + sum3) / sum4;

#ifdef FRANK_WOLFE_DUMP
            cout << "esmat_norm2 (w_minus_s)" << esmat_fnorm (w_minus_s) << endl;
            cout << "esmat_norm2 (w_minus_z)" << esmat_fnorm (w_minus_z) << endl;
#endif
            // cout << "within frank_wolfe_solver: step two finished" << endl;

#ifdef EXACT_LINE_SEARCH_DUMP
            cout << "[exact line search] (sum1, sum2, sum3, sum4, gamma) = ("
                << sum1 << ", " << sum2 << ", " << sum3 << ", " << sum4 << ", " << gamma
                << ")"
                << endl;
#endif
        }

        }
#endif
        // update the w^(k+1)
        esmat_scalar_mult (gamma, s); 
        esmat_scalar_mult (1.0-gamma, w_1);
        esmat_copy (w_1, tempS);
        esmat_add (tempS, s, w_1);

        // compute value of objective function
        penalty = sub1_objective (Y_1, Z_1, w_1, RHO, N);
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



void group_lasso_solver (Esmat* Y_2, Esmat* Z, Esmat* w_2, double RHO, double lambda) {

    // STEP ONE: compute the optimal solution for truncated problem
    Esmat* wbar = esmat_init (w_2);
    esmat_zeros (wbar);
    esmat_scalar_mult (RHO, Z, wbar); // wbar = RHO * z_2
    esmat_sub (wbar, Y_2, wbar); // wbar = RHO * z_2 - y_2
    esmat_scalar_mult (1.0/RHO, wbar); // wbar = (RHO * z_2 - y_2) / RHO

    // STEP TWO: find the closed-form solution for second subproblem
    for (int j = 0; j < w_2->nCols; j ++) {
        // 1. bifurcate the set of values
        vector< pair<int,double> > alpha_vec;
        for (int i = 0; i < w_2->nRows; i ++) {
            double value = wbar[i][j];
            alpha_vec.push_back (make_pair(i, abs(value)));
        }

        // 2. sorting
        std::sort (alpha_vec.begin(), alpha_vec.end(), pairComparator);
        /* 
          // check the decreasing order 
           for (int i = 0; i < N; i ++) {
           if (alpha_vec[i].second != 0)
           cout << alpha_vec[i].second << endl;
           }
           */

        // 3. find mstar
        int mstar = 0; // number of elements support the sky
        double separator;
        double max_term = -INF, new_term;
        double sum_alpha = 0.0;
        for (int i = 0; i < N; i ++) {
            sum_alpha += alpha_vec[i].second;
            new_term = (sum_alpha - lambda) / (i + 1.0);
            if ( new_term > max_term ) {
                separator = alpha_vec[i].second;
                max_term = new_term;
                mstar = i;
            }
        }

        // 4. assign closed-form solution to w_2
        if( max_term < 0 ){
            for(int i = 0; i < N; i++)
                w_2[i][j] = 0.0;
            continue;
        }

        for (int i = 0; i < N; i ++) {
            // harness vector of pair
            double value = wbar[i][j];
            if ( abs(value) >= separator ) {
                w_2[i][j] = max_term;
            } else {
                // its ranking is above m*, directly inherit the wbar
                w_2[i][j] = max(wbar[i][j],0.0);
            }
        }
    }

    // compute value of objective function
    double penalty = sub2_objective (Y_2, Z, w_2, RHO, lambda);
    // report the #iter and objective function
    /*cout << "[Blockwise] sub2_objective: " << penalty << endl;
      cout << endl;*/

    // STEP THREE: recollect temporary variable - wbar
    esmat_free (wbar);
}


/* Compute the mutual distance of input instances contained within "data" */
/*
void compute_match_mat (vector<Instance*>& data, double ** dist_mat, int N, int D, dist_func df, bool isSym) {

    // STEP ONE: convert input sparse data to desired format
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            ;
        }
    }
    // STEP TWO: check correctness by assertion
}
*/

void HDP (int D, int N, vector<double> LAMBDAs, Esmat* W) {
    // SET MODEL-RELEVANT PARAMETERS 
    assert (LAMBDAs.size() == 3);
    double LABMDA_DOC = LAMBDAs[0];
    double LAMBDA_TOPIC = LAMBDAs[1];
    double LAMBDA_BLOCK = LAMBDAs[2];
    double ALPHA = 1.0;
    double RHO = 1.0;

    /* DECLARE AND INITIALIZE INVOLVED VARIABLES AND MATRICES */
    Esmat* w_1 = esmat_init (N, N);
    Esmat* w_2 = esmat_init (N, N);
    Esmat* w_3 = esmat_init (N, N);
    Esmat* w_4 = esmat_init (N, N);

    Esmat* y_1 = esmat_init (N, N);
    Esmat* y_2 = esmat_init (N, N);
    Esmat* y_3 = esmat_init (N, N);
    Esmat* y_4 = esmat_init (N, N);

    Esmat* z = esmat_init (N, N);

    Esmat* diff_1 = esmat_init (N, N);
    Esmat* diff_2 = esmat_init (N, N);
    Esmat* diff_3 = esmat_init (N, N);
    Esmat* diff_4 = esmat_init (N, N);

    /* SET ITERATION-RELEVANT VARIABLES */
    double error = INF;
    int iter = 0; 
    int max_iter = 2000;

    /* ITERATIVE OPTIMIZATION */
    while ( iter < max_iter ) { // STOPPING CRITERIA
        // STEP ZERO: RESET ALL SUBPROBLEM SOLUTIONS (OPTIONAL) 
        esmat_zeros (w_1);
        esmat_zeros (w_2);
        esmat_zeros (w_3);
        esmat_zeros (w_4);
        /*
#ifdef ITERATION_TRACE_DUMP
        cout << "it is place 0 iteration #" << iter << ", going to get into frank_wolfe_solvere"  << endl;
#endif
*/
        // STEP ONE: RESOLVE W_1, W_2, W_3, W_4
        // resolve w_1
        frank_wolfe_solver ( y_1, z, w_1, RHO, N); 

        /*
#ifdef ITERATION_TRACE_DUMP
        cout << "frank_wolfe_solver done. norm2(w_1) = " << mat_norm2 (wone, N, N) << endl;
        cout << "it is place 1 iteration #" << iter << ", going to get into group_lasso_solver"<< endl;
#endif
*/
        // resolve w_2
        group_lasso_solver (y_2, z, w_2, RHO, LAMBDAs[0], N);
        /*
#ifdef ITERATION_TRACE_DUMP
        cout << "it is place 3 iteration #" << iter << endl;
        cout << "norm2(w_2) = " << mat_norm2 (wtwo, N, N) << endl;
#endif
*/
        // resolve w_3
        group_lasso_solver (y_3, z, w_3, RHO, LAMBDAs[1], N);
        // resolve w_4
        group_lasso_solver (y_4, z, w_4, RHO, LAMBDAs[2], N);
        
        // STEP TWO: update z by averaging w_1, w_2, w_3 and w_4
        Esmat* temp = esmat_init (0,0);
        esmat_add (w_1, w_2, z);
        esmat_copy (z, temp);
        esmat_add (temp, w_3, z);
        esmat_copy (z, temp);
        esmat_add (temp, w_4, z);
        esmat_scalar_mult (0.25, z);

        /*
#ifdef ITERATION_TRACE_DUMP
        cout << "it is place 4 iteration #" << iter << endl;
        cout << "norm2(z) = " << mat_norm2 (z, N, N) << endl;
#endif
*/
        // STEP THREE: update the y_1 and y_2 by w_1, w_2 and z
        esmat_sub (w_1, z, diff_1);
        // double trace_wone_minus_z = esmat_norm2 (diff_1); 
        esmat_scalar_mult (ALPHA, diff_1);
        esmat_copy (y_1, temp);
        esmat_add (temp, diff_1, y_1);

        esmat_sub (w_2, z, diff_2);
        //double trace_wtwo_minus_z = esmat_norm2 (diff_2); 
        esmat_scalar_mult (ALPHA, diff_2);
        esmat_copy (y_2, temp);
        esmat_add (temp, diff_2, y_2);

        esmat_sub (w_3, z, diff_3);
        //double trace_wthree_minus_z = esmat_norm2 (diff_3); 
        esmat_scalar_mult (ALPHA, diff_3);
        esmat_copy (y_3, temp);
        esmat_add (temp, diff_3, y_3);

        esmat_sub (w_4, z, diff_4);
        //double trace_wfour_minus_z = esmat_norm2 (diff_4); 
        esmat_scalar_mult (ALPHA, diff_4);
        esmat_copy (y_4, temp);
        esmat_add (temp, diff_4, y_4);

        // STEP FOUR: trace the objective function
        /*
        if (iter % 1 == 0) {
            // 1. trace the error
            error = original_objective (dist_mat, lambda, N, z);
            cout << "[Overall] iter = " << iter 
                 << ", Overall Error: " << error;
            /
#ifdef NTOPIC_DUMP
            // 2. get number of topic
            int nTopics = get_nTopics(z, N, N);
            cout << ", nTopics: " << nTopics;
#endif
            cout << endl;
        }
        */

        iter ++;
    }
    
    // STEP FIVE: memory recollection
    esmat_free (w_1);
    esmat_free (w_2);
    esmat_free (w_3);
    esmat_free (w_4);

    esmat_free (y_1);
    esmat_free (y_2);
    esmat_free (y_3);
    esmat_free (y_4);

    esmat_free (diff_1);
    esmat_free (diff_2);
    esmat_free (diff_3);
    esmat_free (diff_4);
    // STEP SIX: put converged solution to destinated W
    esmat_copy (z, W);
    esmat_free (z);
}

// entry main function
int main (int argc, char ** argv) {

    // EXCEPTION control: illustrate the usage if get input of wrong format
    if (argc < 3) {
        cerr << "Usage: HDP [dataFile] [lambda_document] [lambda_topic] [lambda_block]" << endl;
        cerr << "Note: dataFile must be scaled to [0,1] in advance." << endl;
        exit(-1);
    }

    // PARSE arguments
    char * dataFile = argv[1];
    vector<double> LAMBDAs (3, 0.0);
    LAMBDAs[0] = atof(argv[2]); // lambda_document
    LAMBDAs[1] = atof(argv[3]); // lambda_topic
    LAMBDAs[2] = atof(argv[4]); // lambda_block

    // READ in data
    vector<Instance*> data;
    //read2D (dataFile, data);  
    // NOTE: need to change the number once switch to a new dataset
    readFixDim (dataFile, data, 13);

    // EDA: explore the data 
    int dimensions = -1;
    int N = data.size(); // data size
    for (int i = 0; i < N; i++) {
        vector< pair<int,double> > * f = &(data[i]->fea);
        int last_index = f->size() - 1;
        if (f->at(last_index).first > dimensions) {
            dimensions = f->at(last_index).first;
        }
    }
    int D = dimensions;
    cerr << "D = " << D << endl; // # features
    cerr << "N = " << N << endl; // # instances
    cerr << "lambda_g = " << LAMBDAs[0] << endl;
    cerr << "lambda_l = " << LAMBDAs[1] << endl;
    cerr << "lambda_b = " << LAMBDAs[2] << endl;
    cerr << "r = " << TRIM_THRESHOLD << endl;

    int seed = time(NULL);
    srand (seed);
    cerr << "seed = " << seed << endl;

    // restore matchness matrix in sparse representation
    /* here we consider non-noise version of topic model
    double ** match_mat = mat_init (N, N);
    mat_zeros (match_mat, N, N);
    */

    // Run sparse convex clustering
    Esmat* W = esmat_init (N, 1);
    HDP (D, N, LAMBDAs, W);

    // Output results
    ofstream fout("result");

    // get all topics
    /*{{{*/
    /*
    vector<int> centroids = get_all_centroids(W, N, N); // contains index of all centroids
    
    int nCentroids = centroids.size();
    for (int i = 0; i < N; i ++) {
        // output identification and its belonging
        fout << "id=" << i+1 << ", fea[0]=" << data[i]->fea[0].second << ", ";  // sample id
        for (int j = 0; j < N; j ++) {
            if( fabs(W[i][j]) > 3e-1 ) {
                fout << j+1 << "(" << W[i][j] << "),\t";
            }
        }
	fout << endl;

        // output distance of one sample to each centroid 
        fout << "dist_centroids: (";
        for (int j = 0; j < nCentroids - 1; j ++) {
            fout << dist_mat[i][ centroids[j] ] << ", ";
        }
        fout << dist_mat[i][ centroids[nCentroids-1] ] << ")";
        fout << endl;
    }
    */
/*}}}*/
}
