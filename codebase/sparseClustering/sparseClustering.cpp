/*###############################################################
## MODULE: sparseKmeans.cpp
## VERSION: 1.0 
## SINCE 2014-06-14
## AUTHOR:
##     Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
## DESCRIPTION: 
##      ################################################################# 
## Edited by MacVim
## Class Info auto-generated by Snippet 
################################################################*/

#include "sparseClustering.h"
#include <cassert>

/* algorithmic options */ 
#define EXACT_LINE_SEARCH  // comment this to use inexact search

/* dumping options */
// #define FRANK_WOLFE_DUMP
// #define EXACT_LINE_SEARCH_DUMP
// #define BLOCKWISE_DUMP
#define NCENTROID_DUMP


typedef double (* dist_func) (Instance*, Instance*, int); 
const double r = 1000.0;

int get_nCentroids (double ** W, int nRows, int nCols) {

    int nCentroids = 0;
    double * sum_belonging = new double [nCols];
    for (int i = 0; i < nCols; i ++) {
        sum_belonging[i] = 0.0;
    }

    mat_sum_col (W, sum_belonging, nRows, nCols);

    for (int i = 0; i < nCols; i ++ ) {
        if (sum_belonging[i] > 1) { 
            nCentroids ++;
        }
    }

    return nCentroids;
}

vector<int> get_all_centroids(double ** W, int nRows, int nCols) {

    std::vector<int> centroids;

    double * sum_belonging = new double [nCols];
    for (int i = 0; i < nCols; i ++) {
        sum_belonging[i] = 0.0;
    }

    mat_sum_col (W, sum_belonging, nRows, nCols);

    for (int i = 0; i < nCols; i ++ ) {
        if (sum_belonging[i] > 1) { 
            centroids.push_back(i);
        }
    }
    
    return centroids;
}

double first_subproblm_obj (double ** dist_mat, double ** yone, double ** zone, double ** wone, double rho, int N) {

    double ** temp = mat_init (N, N);
    double ** diffone =mat_init (N, N);
    mat_zeros (diffone, N, N);

    // sum1 = 0.5 * sum_n sum_k (w_nk * d^2_nk) -> loss
    mat_zeros (temp, N, N);
    mat_times (wone, dist_mat, temp, N, N);
    double sum1 = 0.5 * mat_sum (temp, N, N);

    // sum2 = y_1^T dot (w_1 - z) -> linear
    mat_zeros (temp, N, N);
    mat_sub (wone, zone, diffone, N, N); // temp = w_1 - z_1
    mat_tdot (yone, diffone, temp, N, N);
    double sum2 = mat_sum (temp, N, N);

    // sum3 = 0.5 * rho * || w_1 - z_1 ||^2 -> quadratic
    mat_zeros (temp, N, N);
    mat_sub (wone, zone, temp, N, N);
    double sum3 = 0.5 * rho * mat_norm2 (temp, N, N);

    // sum4 = r dot (1 - sum_k w_nk) -> dummy
    double * temp_vec = new double [N];
    mat_sum_row (wone, temp_vec, N, N);
    double dummy_penalty = 0.0;
    for (int i = 0; i < N; i ++) {
        dummy_penalty += r*(1 - temp_vec[i]);
    }
    double total = sum1+sum2+sum3+dummy_penalty;
#ifdef FRANK_WOLFE_DUMP
    cout << "[Frank_wolfe] (loss, linear, quadratic, dummy, total) = (" 
         << sum1 << ", " << sum2 << ", " << sum3 << ", " << dummy_penalty << ", " << total
         <<  ")" << endl;
#endif

    mat_free (temp, N, N);
    mat_free (diffone, N, N);
    delete [] temp_vec;

    return total;
}

void frank_wolf (double ** dist_mat, double ** yone, double ** zone, double ** wone, double rho, int N) {

    bool is_global_optimal_reached = false;

    // cout << "within frank_wolf" << endl;
    // This can be computed by using corner point. 
    double ** gradient = mat_init (N, N);
    double ** s = mat_init (N, N);
    mat_zeros (gradient, N, N);
    mat_zeros (s, N, N);

#ifndef EXACT_LINE_SEARCH
    int K = 300;
#else
    int K = 50;
    double ** w_minus_s = mat_init (N, N);
    double ** w_minus_z = mat_init (N, N);
#endif

    int k = 0; // iteration number
    double gamma; // step size
    double penalty;
    double ** tempS = mat_init (N, N);
    mat_zeros (wone, N, N);
    double sum1, sum2, sum3, sum4;
    // cout << "within frank_wolf: start iteration" << endl;
    while (k < K && !is_global_optimal_reached) {
        // STEP ONE: find s minimize <s, grad f>
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                gradient[i][j] = 0.5 * dist_mat[i][j] + yone[i][j] + rho * (wone[i][j] - zone[i][j]) ;//- r[i];
            }
        }
        mat_min_row (gradient, s, N, N);

#ifdef FRANK_WOLFE_DUMP
        cout << "mat_norm2 (wone, N, N): " <<  mat_norm2 (wone, N, N) << endl;
        cout << "mat_norm2 (yone, N, N): " <<  mat_norm2 (yone, N, N) << endl;
        cout << "mat_norm2 (gradient, N, N): " <<  mat_norm2 (gradient, N, N) << endl;
        cout << "mat_sum (s, N, N): " <<  mat_sum (s, N, N) << endl;
#endif
        // cout << "within frank_wolf: step one finished" << endl;

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
        // sum3 = - rho * sum_n sum_k  (w - z) 
        // sum4 = sum_n sum_k rho * (s - w)
        mat_zeros (tempS, N, N);
        mat_zeros (w_minus_s, N, N);
        mat_zeros (w_minus_z, N, N);
        mat_sub (wone, s, w_minus_s, N, N);
        mat_sub (wone, zone, w_minus_z, N, N);

        // NOTE: in case of ||w_1 - s||^2 = 0, not need to optimize anymore
        // since incremental term = w + gamma (s - w), and whatever gamma is,
        // w^(k+1) = w^(k), this would be equivalent to gamma = 0
        if (mat_norm2(w_minus_s, N, N) == 0) {
            gamma = 0;
            is_global_optimal_reached = true;
        } else {
            mat_times (w_minus_s, dist_mat, tempS, N, N);
            sum1 = 0.5 * mat_sum (tempS, N, N);

            mat_zeros (tempS, N, N);
            mat_tdot (yone, w_minus_s, tempS, N, N);
            sum2 = mat_sum (tempS, N, N);

            mat_zeros (tempS, N, N);
            mat_tdot (w_minus_z, w_minus_s, tempS, N, N);
            sum3 = rho * mat_sum (tempS, N, N);

            mat_zeros (tempS, N, N);
            sum4 = rho * mat_norm2 (w_minus_s, N, N);

            // gamma should be within interval [0,1]
            gamma = (sum1 + sum2 + sum3) / sum4;

#ifdef FRANK_WOLFE_DUMP
            cout << "mat_norm2 (w_minus_s, N, N)" << mat_norm2 (w_minus_s, N, N) << endl;
            cout << "mat_norm2 (w_minus_z, N, N)" << mat_norm2 (w_minus_z, N, N) << endl;
#endif
            // cout << "within frank_wolf: step two finished" << endl;

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
        mat_zeros (tempS, N, N);
        mat_dot (gamma, s, tempS, N, N);
        mat_dot (1.0-gamma, wone, wone, N, N);
        mat_add (wone, tempS, wone, N, N);

        // compute value of objective function
        penalty = first_subproblm_obj (dist_mat, yone, zone, wone, rho, N);
    // cout << "within frank_wolf: step three finished" << endl;
        // report the #iter and objective function
        /*
        cout << "[Frank-Wolfe] iteration: " << k << ", first_subpro_obj: " << penalty << endl;
        */

        k ++;
    // cout << "within frank_wolf: next iteration" << endl;
    }
    // cout << "within frank_wolf: to free gradient" << endl;
    mat_free (gradient, N, N);
    // cout << "within frank_wolf: to free temps" << endl;
    mat_free (tempS, N, N);
    // cout << "within frank_wolf: to free s " << endl;
    mat_free (s, N, N);
    // cout << "end frank_wolf: finished! " << endl;
}

double sign (int input) {

    if (input > 0) return 1.0;
    else if ( input < 0 ) return -1.0;
    else return 0.0;

}

bool pairComparator (const std::pair<int, double>& firstElem, const std::pair<int, double>& secondElem) {
    // sort pairs by second element with decreasing order
    return firstElem.second > secondElem.second;
}

double second_subproblem_obj (double ** ytwo, double ** z, double ** wtwo, double rho, int N, double lambda) {

    double ** temp = mat_init (N, N);
    double ** difftwo = mat_init (N, N);
    mat_zeros (difftwo, N, N);

    // reg = 0.5 * sum_k max_n | w_nk |  -> group-lasso
    mat_zeros (temp, N, N);
    double * maxn = new double [N]; 
    for (int i = 0; i < N; i ++) { // Ian: need initial 
        maxn[i] = -INF;
    }

    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            if (wtwo[i][j] > maxn[j])
                maxn[j] = wtwo[i][j];
        }
    }
    double sumk = 0.0;
    for (int i = 0; i < N; i ++) {
        sumk += maxn[i];
    }
    double group_lasso = lambda * sumk; 

    // sum2 = y_2^T dot (w_2 - z) -> linear
    mat_zeros (temp, N, N);
    mat_sub (ytwo, z, difftwo, N, N);
    mat_tdot (ytwo, difftwo, temp, N, N);
    double sum2 = mat_sum (temp, N, N);

    // sum3 = 0.5 * rho * || w_2 - z_2 ||^2 -> quadratic mat_zeros (temp, N, N);
    mat_sub (wtwo, z, temp, N, N);
    double sum3 = 0.5 * rho * mat_norm2 (temp, N, N);

    mat_free (temp, N, N);

    // ouput values of each components
#ifdef BLOCKWISE_DUMP
    cout << "[Blockwise] (group_lasso, linear, quadratic) = ("
        << group_lasso << ", " << sum2 << ", " << sum3
        << ")" << endl;
#endif

    //cerr << group_lasso << ", " << sum2 << ", " << sum3 << endl;
    return group_lasso + sum2 + sum3;
}

void blockwise_closed_form (double ** ytwo, double ** ztwo, double ** wtwo, double rho, double lambda, int N) {

    // STEP ONE: compute the optimal solution for truncated problem
    double ** wbar = mat_init (N, N);
    mat_zeros (wbar, N, N);
    mat_dot (rho, ztwo, wbar, N, N); // wbar = rho * z_2
    mat_sub (wbar, ytwo, wbar, N, N); // wbar = rho * z_2 - y_2
    mat_dot (1.0/rho, wbar, wbar, N, N); // wbar = (rho * z_2 - y_2) / rho

    // STEP TWO: find the closed-form solution for second subproblem
    for (int j = 0; j < N; j ++) {
        // 1. bifurcate the set of values
        vector< pair<int,double> > alpha_vec;
        for (int i = 0; i < N; i ++) {
            double value = wbar[i][j];
            /*if( wbar[i][j] < 0 ){
              cerr << "wbar[" << i << "][" << j << "]" << endl;
              exit(0);
              }*/
            alpha_vec.push_back (make_pair(i, abs(value)));
        }

        // 2. sorting
        std::sort (alpha_vec.begin(), alpha_vec.end(), pairComparator);
        /*
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

        // 4. assign closed-form solution to wtwo
        if( max_term < 0 ){
            for(int i=0;i<N;i++)
                wtwo[i][j] = 0.0;
            continue;
        }

        for (int i = 0; i < N; i ++) {
            // harness vector of pair
            double value = wbar[i][j];
            if ( abs(value) >= separator ) {
                wtwo[i][j] = max_term;
            } else {
                // its ranking is above m*, directly inherit the wbar
                wtwo[i][j] = max(wbar[i][j],0.0);
            }
        }
    }

    // compute value of objective function
    double penalty = second_subproblem_obj (ytwo, ztwo, wtwo, rho, N, lambda);
    // report the #iter and objective function
    /*cout << "[Blockwise] second_subproblem_obj: " << penalty << endl;
      cout << endl;*/

    // STEP THREE: recollect temporary variable - wbar
    mat_free (wbar, N, N);
}

double L2norm (Instance * ins1, Instance * ins2, int D) {
    // TODO: 
    //   1. refine by using hash table to restore each instance
    //   2. avoid the apply memory for vec1, vec2, make it direct computation
    assert (ins1->fea.size() == D);
    assert (ins2->fea.size() == D);

    double * diff = new double [D];

    for (int i = 0; i < D; i ++) {
        diff[ ins1->fea[i].first ] = ins1->fea[i].second;
    }
    for (int i = 0; i < D; i ++) {
        diff[ ins2->fea[i].first ] -= ins2->fea[i].second;
    }

    double norm = 0.0;
    for (int i = 0; i < D; i ++) {
        norm += diff[i] * diff[i];
    }

    return norm;
}

double opt_objective (double ** dist_mat, double lambda, int N, double ** z) {
    // N is number of entities in "data", and z is N by N.
    // z is current valid solution (average of w_1 and w_2)

    // STEP ONE: compute loss function
    double normSum = 0.0;
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            normSum += z[i][j] * dist_mat[i][j];
        }
    }

    // sum4 = r dot (1 - sum_k w_nk) -> dummy
    double * temp_vec = new double [N];
    mat_sum_row (z, temp_vec, N, N);
    double dummy_penalty=0.0;
    double avg=0.0;
    for (int i = 0; i < N; i ++) {
        avg += temp_vec[i];
        dummy_penalty += r * max(1 - temp_vec[i], 0.0) ;
    }

    double loss = 0.5 * (normSum+dummy_penalty);
    // cout << "loss=" << loss << endl;

    // STEP TWO: compute group-lasso regularization
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
        sumk += maxn[i];
    }
    double reg = lambda * sumk; 
    //cout << "reg=" << reg << endl;

    delete[] temp_vec;
    return loss + reg;
}

/* Compute the mutual distance of input instances contained within "data" */
void compute_dist_mat (vector<Instance*>& data, double ** dist_mat, int N, int D, dist_func df, bool isSym) {

    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {

            if (j >= i || !isSym) { // invoke dist_func
                if (j == i) {
                    dist_mat[i][j] = 0;
                    continue;
                }
                Instance * xi = data[i];
                Instance * muj = data[j];
                dist_mat[i][j] = df (xi, muj, D);
            } else { // by symmetry 
                dist_mat[i][j] = dist_mat[j][i];
            }

        }
    }

}

void sparseClustering ( double ** dist_mat, int D, int N, double lambda, double ** W) {

    // parameters 
    double alpha = 1.0;
    double rho = 1.0;

    // iterative optimization 
    double error = INF;
    double ** wone = mat_init (N, N);
    double ** wtwo = mat_init (N, N);
    double ** yone = mat_init (N, N);
    double ** ytwo = mat_init (N, N);
    double ** z = mat_init (N, N);
    double ** diffzero = mat_init (N, N);
    double ** diffone = mat_init (N, N);
    double ** difftwo = mat_init (N, N);
    mat_zeros (wone, N, N);
    mat_zeros (wtwo, N, N);
    mat_zeros (yone, N, N);
    mat_zeros (ytwo, N, N);
    mat_zeros (z, N, N);
    mat_zeros (diffzero, N, N);
    mat_zeros (diffone, N, N);
    mat_zeros (difftwo, N, N);


    int iter = 0; // Ian: usually we count up (instead of count down)
    int max_iter = 300;

    while ( iter < max_iter ) { // stopping criteria

#ifdef SPARSE_CLUSTERING_DUMP
        cout << "it is place 0 iteration #" << iter << ", going to get into frank_wolfe"  << endl;
#endif

        // STEP ONE: resolve w_1 and w_2
        frank_wolf (dist_mat, yone, z, wone, rho, N); // for w_1
#ifdef SPARSE_CLUSTERING_DUMP
        cout << "frank_wolfe done. norm2(w_1) = " << mat_norm2 (wone, N, N) << endl;
        cout << "it is place 1 iteration #" << iter << ", going to get into blockwise_closed_form"<< endl;
#endif

        blockwise_closed_form (ytwo, z, wtwo, rho, lambda, N);  // for w_2
#ifdef SPARSE_CLUSTERING_DUMP
        cout << "it is place 3 iteration #" << iter << endl;
        cout << "norm2(w_2) = " << mat_norm2 (wtwo, N, N) << endl;
#endif
        /* // test
        mat_sub (wone, wtwo, diffzero, N, N);
        double trace_wone_minus_wtwo = mat_norm2 (diffzero, N, N);
        */
        
        // STEP TWO: update z by averaging w_1 and w_2
        mat_add (wone, wtwo, z, N, N);
        mat_dot (0.5, z, z, N, N);
#ifdef SPARSE_CLUSTERING_DUMP
        cout << "it is place 4 iteration #" << iter << endl;
        cout << "norm2(z) = " << mat_norm2 (z, N, N) << endl;
#endif

        // STEP THREE: update the y_1 and y_2 by w_1, w_2 and z
        double * temp_vec = new double [N];
        mat_sum_row (z, temp_vec, N, N);
        double dummy_penalty=0.0;
        for (int i = 0; i < N; i ++) {
            /*if( temp_vec[i] > 1 ){
              cerr << "tmp_vec[i]=" << temp_vec[i] << endl;
              mat_sum_row(ytwo, temp_vec, N, N);
              cerr << "y1_sum[i]=" << temp_vec[i] << endl;
              mat_sum_row(wtwo, temp_vec, N, N);
              cerr << "w2_sum[i]=" << temp_vec[i] << endl;
              exit(0);
              }*/
            dummy_penalty += r*( 1 - temp_vec[i] );
        }
        delete[] temp_vec;

        mat_sub (wone, z, diffone, N, N);
        //double trace_wone_minus_z = mat_norm2 (diffone, N, N); 
        mat_dot (alpha, diffone, diffone, N, N);
        mat_add (yone, diffone, yone, N, N);

        mat_sub (wtwo, z, difftwo, N, N);
        //double trace_wtwo_minus_z = mat_norm2 (difftwo, N, N); 
        mat_dot (alpha, difftwo, difftwo, N, N);
        mat_add (ytwo, difftwo, ytwo, N, N);

        // STEP FOUR: trace the objective function
        if (iter % 1 == 0) {
            error = opt_objective (dist_mat, lambda, N, z);
            // get current number of employed centroid
#ifdef NCENTROID_DUMP
            int nCentroids = get_nCentroids (z, N, N);
#endif
            cout << "[Overall] iter = " << iter 
                 << ", Overall Error: " << error 
#ifdef NCENTROID_DUMP
                 << ", nCentroids: " << nCentroids
#endif
                 << endl;
        }
        /*cout << "w1" << endl;
          error = opt_objective (dist_mat, lambda, N, wone);
          cout << "w2" << endl;
          error = opt_objective (dist_mat, lambda, N, wtwo);*/
        /*cout << "[Overall] || w_1 - w_2 || ^2 = " << trace_wone_minus_wtwo << endl;
          cout << "[Overall] || w_1 - z || ^2 = " << trace_wone_minus_z << endl;
          cout << "[Overall] || w_2 - z || ^2 = " << trace_wtwo_minus_z << endl;
          cout << endl;*/

        iter ++;
    }
    
    // STEP FIVE: memory recollection
    mat_free (wone, N, N);
    mat_free (wtwo, N, N);
    mat_free (yone, N, N);
    mat_free (ytwo, N, N);
    mat_free (diffone, N, N);
    mat_free (difftwo, N, N);
    // STEP SIX: put converged solution to destination W
    mat_copy (z, W, N, N);
}

// entry main function
int main (int argc, char ** argv) {

    // exception control: illustrate the usage if get input of wrong format
    if (argc < 3) {
        cerr << "Usage: sparseClustering [dataFile] [lambda]" << endl;
        cerr << "Note: dataFile must be scaled to [0,1] in advance." << endl;
        exit(-1);
    }

    // parse arguments
    char * dataFile = argv[1];
    double lambda = atof(argv[2]);

    // read in data
    vector<Instance*> data;
    //read2D (dataFile, data);
    readFixDim (dataFile, data, 13);

    // explore the data 
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
    cerr << "lambda = " << lambda << endl;
    cerr << "r = " << r << endl;
    int seed = time(NULL);
    srand (seed);
    cerr << "seed = " << seed << endl;

    // pre-compute distance matrix
    dist_func df = L2norm;
    double ** dist_mat = mat_init (N, N);
    mat_zeros (dist_mat, N, N);
    compute_dist_mat (data, dist_mat, N, D, df, true); 

    // Run sparse convex clustering
    map<int, Cluster*> clusters; 
    double ** W = mat_init (N, N);
    mat_zeros (W, N, N);
    sparseClustering (dist_mat, D, N, lambda, W);
    // Output results
    ofstream fout("result");

    // get all centroids
    vector<int> centroids = get_all_centroids(W, N, N); // contains index of all centroids
    
    int nCentroids = centroids.size();
    for (int i = 0; i < N; i ++) {
        // output identification and its belonging
        fout << "id=" << i+1 << ", fea[0]=" << data[i]->fea[0].second << ", ";  // sample id
        for (int j = 0; j < N; j ++) {
            if( fabs(W[i][j]) > 3e-1 ) {
                fout << "centroid: " << j+1 << "(" << W[i][j] << "),\t";
            }
        }
        // output distance of one sample to each centroid 
        fout << "dist_centroids: (";
        for (int j = 0; j < nCentroids - 1; j ++) {
            fout << dist_mat[i][ centroids[j] ] << ", ";
        }
        fout << dist_mat[i][ centroids[nCentroids-1] ] << ")";

        fout << endl;
    }
}
