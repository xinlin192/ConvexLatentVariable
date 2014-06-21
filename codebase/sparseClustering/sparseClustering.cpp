/*###############################################################
## MODULE: sparseKmeans.cpp
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

#include "sparseClustering.h"

typedef double (* dist_func) (Instance*, Instance*, int);

double first_subproblm_obj (double ** dist_mat, double ** yone, double ** zone, double ** wone, double rho, int N, double * r) {

    double ** temp = mat_init (N, N);
    
    // sum1 = 0.5 * sum_n sum_k (w_nk * d^2_nk)
    mat_times (wone, dist_mat, temp, N, N);
    double sum1 = 0.5 * mat_sum (temp, N, N);

    // sum2 = y_1^T dot w_1
    mat_tdot (yone, wone, temp, N, N);
    double sum2 = mat_sum (temp, N, N);

    // sum3 = 0.5 * rho * || w_1 - z_1 ||^2
    mat_sub (wone, zone, temp, N, N);
    double sum3 = 0.5 * rho* mat_norm2 (temp, N, N);
    
    // sum4 = r dot (1 - sum_k w_nk)
    double * temp_vec = new double [N];
    mat_sum_row (wone, temp_vec, N, N);
    for (int i = 0; i < N; i ++) {
        temp_vec[i] = 1 - temp_vec[i];
    }

    double sum4 = mat_dot (temp_vec, r, N);
    cout << "sum1: " << sum1 << endl;
    cout << "sum2: " << sum2 << endl;
    cout << "sum3: " << sum3 << endl;
    cout << "dummy term: " << sum4 << endl;

    mat_free (temp, N, N);
    delete temp_vec;

    return sum1 + sum2 + sum3 + sum4;
}

void frank_wolf (double ** dist_mat, double ** yone, double ** zone, double ** wone, double rho, int N) {

    // This can be computed by using corner point. 
    double ** gradient = mat_init (N, N);
    double ** s = mat_init (N, N);
    double * r = new double [N];
    for (int i = 0; i < N; i ++) {
        r[i] = 100;
    }
    
    int K = 100, k = 0; // iteration number
    double gamma; // step size
    double penalty;
    double ** tempS = mat_init(N, N);
    while (k < K) {
        // STEP ONE: find s minimize <s, grad f>
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                gradient[i][j] = 0.5 * dist_mat[i][j] + yone[i][j] + rho * (wone[i][j] - zone[i][j]) ;//- r[i];
            }
        }
        mat_min_row (gradient, s, N, N);
	
        // STEP TWO: apply exact or inexact line search to find solution
        // TODO: refine it by using exact line search algorithm
        // Here we use inexact line search
        gamma = 2.0 / (k + 2.0);
        mat_dot (gamma, s, tempS, N, N);
        mat_dot (1.0-gamma, wone, wone, N, N);
        mat_add (wone, tempS, wone, N, N);
	
        // compute value of objective function
        penalty = first_subproblm_obj (dist_mat, yone, zone, wone, rho, N, r);
        // report the #iter and objective function
        cout << "[Frank-Wolfe] iteration: " << k << ", first_subproblm_obj: " << penalty << endl;

        k ++;
    }
    mat_free (gradient, N, N);
    mat_free (tempS, N, N);
    mat_free (s, N, N);
}

double sign (int input) {

    if (input >= 0) return 1.0;
    else return -1.0;

}

bool pairComparator (const std::pair<int, double>& firstElem, const std::pair<int, double>& secondElem) {
    // sort pairs by second element with decreasing order
    return firstElem.second > secondElem.second;
}

double second_subproblm_obj (double ** ytwo, double ** z, double ** wtwo, double rho, int N, double lambda) {

    double ** temp = mat_init (N, N);

    // reg = 0.5 * sum_k max_n | w_nk | 
    double * maxn = new double [N]; 
    for (int i = 0;i < N; i ++) { // Ian: need initial 
    	maxn[i] = -INF;
    }

    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            if (z[j][i] > maxn[i])
                maxn[i] = z[j][i];
        }
    }
    double sumk = 0.0;
    for (int i = 0; i < N; i ++) {
        sumk = maxn[i];
    }
    double reg = lambda * sumk; 

    // sum2 = y_2^T dot w_2
    mat_tdot (ytwo, wtwo, temp, N, N);
    double sum2 = mat_sum (temp, N, N);

    // sum3 = 0.5 * rho * || w_2 - z_2 ||^2
    mat_sub (wtwo, z, temp, N, N);
    double sum3 = 0.5 * rho * mat_norm2 (temp, N, N);

    mat_free (temp, N, N);

    return reg + sum2 + sum3;
}

void blockwise_closed_form (double ** ytwo, double ** ztwo, double ** wtwo, double rho, double lambda, int N) {

    // STEP ONE: compute the optimal solution for truncated problem
    double ** wbar = mat_init (N, N);
    mat_dot (rho, ztwo, wbar, N, N); // wbar = rho * z_2
    mat_sub (wbar, ytwo, wbar, N, N); // wbar = rho * z_2 - y_2
    mat_dot (1.0/rho, wbar, wbar, N, N); // wbar = (rho * z_2 - y_2) / rho

    // STEP TWO: find the closed-form solution for second subproblem
    for (int j = 0; j < N; j ++) {
        // 1. bifurcate the set of values
        vector< pair<int,double> > alpha_vec;
        for (int i = 0; i < N; i ++) {
            double value = wbar[i][j];
            alpha_vec.push_back (make_pair(i, abs(value)));
        }

        // 2. sorting
        std::sort (alpha_vec.begin(), alpha_vec.end(), pairComparator);

        // 3. find mstar
        int mstar = 0; // number of elements support the sky
        double separator;
        double old_term = -INF, new_term;
        double sum_alpha = 0.0;
        for (int i = 0; i < N; i ++) {
            sum_alpha += alpha_vec[i].second;
            new_term = (sum_alpha - lambda) / (i + 1.0);
            if ( new_term < old_term ) {
                separator = alpha_vec[i].second;
                break;
            }
            mstar ++;
            old_term = new_term;
        }
        double max_term = old_term; 

        // 4. assign closed-form solution to wtwo
        for (int i = 0; i < N; i ++) {
            // harness vector of pair
            double value = wbar[i][j];
            if ( abs(value) > separator ) { 
                wtwo[i][j] = sign(wbar[i][j]) * max_term;
            } else {
                // its ranking is above m*, directly inherit the wbar
                wtwo[i][j] = wbar[i][j];
            }
        }
    }

    // compute value of objective function
    double penalty = second_subproblm_obj (ytwo, ztwo, wtwo, rho, N, lambda);
    // report the #iter and objective function
    cout << "[Blockwise]  second_subproblem_obj: " << penalty << endl;

    // STEP THREE: recollect temporary variable - wbar
    mat_free (wbar, N, N);

}

double L2norm (Instance * ins1, Instance * ins2, int N) {
    // TODO: 
    //   1. refine by using hash table to restore each instance
    //   2. avoid the apply memory for vec1, vec2, make it direct computation

    double * vec1 = new double [N];
    double * vec2 = new double [N];

    for (int i = 0; i < ins1->fea.size(); i ++) {
        vec1[ ins1->fea[i].first ] = ins1->fea[i].second;
    }
    for (int i = 0; i < ins2->fea.size(); i ++) {
        vec2[ ins2->fea[i].first ] = ins2->fea[i].second;
    }
    
    double norm = 0.0;
    for (int i = 0; i < N; i ++) {
        norm += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }

    delete vec1, vec2;

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
    double loss = 0.5 * normSum;

    // STEP TWO: compute group-lasso regularization
    double * maxn = new double [N]; 
    for (int i = 0;i < N; i ++) { // Ian: need initial 
    	maxn[i] = -INF;
    }

    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            if (z[j][i] > maxn[i])
                maxn[i] = z[j][i];
        }
    }
    double sumk = 0.0;
    for (int i = 0; i < N; i ++) {
        sumk = maxn[i];
    }
    double reg = lambda * sumk; 

    return loss + reg;
}

/* Compute the mutual distance of input instances contained within "data" */
void compute_dist_mat (vector<Instance*>& data, double ** dist_mat, int N, dist_func df, bool isSym) {

    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {

            if (j >= i || !isSym) { // invoke dist_func
                Instance * xi = data[i];
                Instance * muj = data[j];
                dist_mat[i][j] = df (xi, muj, N);
            } else { // by symmetry 
                dist_mat[i][j] = dist_mat[j][i];
            }

        }
    }

}

void sparseClustering ( vector<Instance*>& data, int D, int N, double lambda, double ** W) {

    // parameters 
    double alpha = 0.1;
    double rho = 1.0;
    dist_func df = L2norm;

    // iterative optimization 
    double error = INF;
    double ** wone = mat_init (N, N);
    double ** wtwo = mat_init (N, N);
    double ** yone = mat_init (N, N);
    double ** ytwo = mat_init (N, N);
    double ** z = mat_init (N, N);
    double ** diffone = mat_init (N, N);
    double ** difftwo = mat_init (N, N);
    
    double ** dist_mat = mat_init (N, N);
    compute_dist_mat (data, dist_mat, N, df, true); 

    int iter = 0; // Ian: usually we count up (instead of count down)
    int max_iter = 1000;

    while ( iter < max_iter ) { // stopping criteria
        // STEP ONE: resolve w_1 and w_2
        frank_wolf (dist_mat, yone, z, wone, rho, N); // for w_1
        blockwise_closed_form (ytwo, z, wtwo, rho, lambda, N);  // for w_2

        // STEP TWO: update z by averaging w_1 and w_2
        mat_add (wone, wtwo, z, N, N);
        mat_dot (0.5, z, z, N, N);

        // STEP THREE: update the y_1 and y_2 by w_1, w_2 and z
        mat_sub (wone, z, diffone, N, N);
        double trace_wone_minus_z = mat_norm2 (diffone, N, N); 
        mat_dot (alpha, diffone, diffone, N, N);
        mat_sub (yone, diffone, yone, N, N);

        mat_sub (wtwo, z, difftwo, N, N);
        double trace_wtwo_minus_z = mat_norm2 (diffone, N, N); 
        mat_dot (alpha, difftwo, difftwo, N, N);
        mat_sub (ytwo, difftwo, ytwo, N, N);

        // STEP FOUR: trace the objective function
        error = opt_objective (dist_mat, lambda, N, z);
        cout << "[Overall] iter = " << iter << ", Overall Error: " << error << endl;
        cout << "[Overall] || w_1 - z || ^2 = " << trace_wone_minus_z << endl;
        cout << "[Overall] || w_2 - z || ^2 = " << trace_wtwo_minus_z << endl;

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
    read2D (dataFile, data);

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
    cout << "D = " << D << endl; // # features
    cout << "N = " << N << endl; // # instances
    cout << "lambda = " << lambda << endl;
    int seed = time(NULL);
    srand(seed);
    cout << "seed = " << seed << endl;

    // Run sparse convex clustering
    map<int, Cluster*> clusters; 
    double ** W = mat_init (N, N);
    sparseClustering (data, D, N, lambda, W);
    // Output results
    
}
