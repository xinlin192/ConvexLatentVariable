#include "cvx_clustering.h"
ofstream ss_out ("../obj_vs_time_dp/iris/CVX-DP-MEDOID");

/* algorithmic options */ 
#define EXACT_LINE_SEARCH  // comment this to use inexact search

/* dumping options */
// #define EXACT_LINE_SEARCH_DUMP
//#define LAMBDA_K_PLOT

const double NOISE_EPS = 0.0;
const double FRANK_WOLFE_TOL = 1e-20;
const double ADMM_EPS = 0.000001;
const double SPARSITY_TOL = 1e-5;
const double r = 10000000.0;

// single column skyline 
void skyline (double* wout, double* wbar, int R, double lambda) {
    vector< double > alpha_vec (R, 0.0);
    int num_alpha_elem = 0;
    for (int i = 0; i < R; i ++) {
        if (wbar[i] > SPARSITY_TOL) {
            alpha_vec.push_back (abs(wbar[i]));
            ++ num_alpha_elem;
        }
    }
    double max_term = -1e50;
    double separator = 0.0;

    if (num_alpha_elem == 0) return;
    // 2. sorting
    std::sort (alpha_vec.begin(), alpha_vec.end(), double_dec_comp);
    // 3. find mstar
    int mstar = 0; // number of elements support the sky
    double new_term, sum_alpha = 0.0;
    for (int i = 0; i < num_alpha_elem; i ++) {
        sum_alpha += alpha_vec[i];
        new_term = (sum_alpha - lambda) / (i + 1.0);
        if ( new_term > max_term ) {
            separator = alpha_vec[i];
            max_term = new_term;
            mstar = i;
        }
    }
    if (max_term < 0) 
        max_term = (sum_alpha - lambda) / R;

    for (int i = 0; i < R; i ++) {
        // 4. assign closed-form solution to wout
        set<int>::iterator it;
        if ( max_term < 0 ) {
            wout[i] = 0.0;
            continue;
        }
        double wbar_val = wbar[i];
        if ( abs(wbar_val) >= separator ) {
            wout[i] = max_term;
        } else { 
            // its ranking is above m*, directly inherit the wbar
            wout[i] = max(wbar_val,0.0);
        }
    }
}

void col_BCD_solver (double* d, double ** w, double * y, double* s, int N, double rho, int active_col, double lambda, double & PG) {
    double * g = new double [N];
    double * w_half = new double [N];
    double * w_new = new double [N];

    // compute: g_j
    for (int i = 0; i < N; i ++) {
        g[i] = rho * (s[i] + y[i]/rho - 1) + d[i];
    }
    // update: W_j = W_j - g_j /rho
    for (int i = 0; i < N; i++) {
        w_half[i] =  w[active_col][i] - g[i] / rho;
    }
    // compute: w_{t+1} through skyline TODO
    skyline(w_new, w_half, N, lambda);

    // compute: PG
    PG = 0.0;
    for (int i = 0; i < N; i++) {
        double tmp = w_new[i] - w[active_col][i];
        PG += tmp * tmp;
    }

    // update the s_{t+1}
    for (int i = 0; i < N; i++) {
        s[i] += w_new[i] - w[active_col][i];
    }
    // update back to matrix w
    for (int i = 0; i < N; i ++) {
        w[active_col][i] = w_new[i];
    }

    delete[] w_new;
    delete[] w_half;
    delete[] g;
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
    double loss =  (normSum);
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
    for (int i = 0;i < N; i ++) 
        maxn[i] = -INF;
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

bool is_col_zero (double ** mat, int col_index, int R) {
    bool isallzero = true;
    for (int i = 0; i < R; i ++) {
        if (abs(mat[col_index][i]) > SPARSITY_TOL) {
            isallzero = false;
            break;
        }
    }
    return isallzero;
}

void cvx_clustering (double ** dist_mat, int D, int N, double lambda, double ** W) {
    // parameters 
    double alpha = 0.2;
    double rho = 1;
    ss_out << "Time Objective" << endl;
    double cputime = 0;
    double last_cost = INF;
    double prev = omp_get_wtime();
    // iterative optimization 
    double error = INF;
    double ** w = mat_init (N, N);
    double * s = new double [N];
    double * y = new double [N];
    for (int i = 0; i < N; i++) {
        s[i] = 0.0; // summation of each row
        y[i] = 0.0;
    }
    cerr << "cvx_clustering init" <<endl;

    // variables for shriking method
    vector<int> active_set(N, 0);
    int active_size = active_set.size();
    // set initial active_set as all elements
    for (int i = 0; i < N; i++) 
        active_set[i] = i;

    cputime += omp_get_wtime() - prev;
    ss_out << cputime << " " << 0 << endl;
    prev = omp_get_wtime();
    int iter = 0; 
    bool no_active_element = false, admm_opt_reached = false;
    double PG, PGmax; // TODO
    double INNER_EPS = 0.5;
    int max_inner_iter;
    int active_size_last = N;

    bool greedy = false;
    while ( true ) { // stopping criteria

        // STEP ONE:  greedily or randomly process active column
        active_size = N;
        random_shuffle (active_set.begin(),active_set.end());
        cerr << "shuffled .." << endl;

        // STEP TWO: blockwise coordinate descent (BCD)
        max_inner_iter = N/active_size_last;
        int inner_iter= 0;
        while( inner_iter < max_inner_iter ){
            PGmax = -1e300;
            for(int j=0;j<active_size;j++){
                int k = active_set[j]; // column to be solved
                col_BCD_solver (dist_mat[k], w, y, s, N, rho, k, lambda, PG); 
                if( is_col_zero(w, k, N) && PG <= 0.0 ){
                    active_size--;
                    // swap(active_set[j], active_set[active_size]);
                    int tmp = active_set[j];
                    active_set[j] = active_set[active_size];
                    active_set[active_size] = tmp;
                    continue;
                }
                if ( PG > PGmax ) PGmax = PG;
            }
            if( PGmax < INNER_EPS )//stopping condition
                break;
            inner_iter++;
        }
        active_size_last = active_size;
        // TODO: y update
        for (int i = 0; i < N; i ++) {
            y[i] += rho * (s[i] - 1 );
        }
        // STEP THREE: trace the objective function
        // if (iter < 3 * SS_PERIOD || (iter+1) % SS_PERIOD == 0) {
        if (true) {
            cputime += (omp_get_wtime() - prev);
            error = overall_objective (dist_mat, lambda, N, w);
            cout << "[Overall] iter = " << iter 
                << ", Loss Error: " << error
                << ", cputime: " << cputime << endl;
            // NOTE: here is just for plotting
            if (error >= last_cost) error = last_cost;
            else last_cost = error;
            ss_out << cputime << " " << error << endl;
            prev = omp_get_wtime();
        }

        // STEP FOUR: stopping criteria 

        iter ++;
    }

    // STEP FIVE: memory recollection
    delete[] y;
    delete[] s;
    mat_free (w, N, N);
    // STEP SIX: put converged solution to destination W
    mat_copy (w, W, N, N);
    ss_out.close();
}

// entry main function
int main (int argc, char ** argv) {
    // exception control: illustrate the usage if get input of wrong format
    if (argc < 3) {
        cerr << "Usage: cvx_clustering [dataFile] [lambda] " << endl;
        cerr << "Note: dataFile must be scaled to [0,1] in advance." << endl;
        exit(-1);
    }

    // parse arguments
    char * dataFile = argv[1];
    double lambda_base = atof(argv[2]);

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
    double lambda = lambda_base;

    // pre-compute distance matrix
    dist_func df = L2norm;
    double ** dist_mat = mat_init (N, N);
    compute_dist_mat (data, dist_mat, N, D, df, true); 
    // add noise
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            dist_mat[i][j] += NOISE_EPS* ((double)rand() / RAND_MAX);
        }
    }

    ofstream dmat_out ("dist_mat");
    dmat_out << mat_toString(dist_mat, N, N);
    dmat_out.close();
    //  double ** dist_mat = mat_read (dmatFile, N, N);

    // Run sparse convex clustering
#ifndef LAMBDA_K_PLOT
    cerr << "D = " << D << endl; // # features
    cerr << "N = " << N << endl; // # instances
    cerr << "lambda = " << lambda << endl;
    cerr << "r = " << r << endl;
    int seed = time(NULL);
    srand (seed);
    cerr << "seed = " << seed << endl;
    cerr << "==================================================" << endl; 


    double ** W = mat_init (N, N);
    cvx_clustering (dist_mat, D, N, lambda, W);

    // Output cluster
    // output_objective(clustering_objective (dist_mat, W, N));
    /* Output cluster centroids */
    // output_model (W, N);
    /* Output assignment */
    // output_assignment (W, data, N);

    // output the examplar 
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
                // cout << "member: "<< i << endl;
            }
        }
    }
    // compute the means
    vector< vector<double> > means (nCentroids, vector<double>(D, 0.0));
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
    double means_reg = lambda * nCentroids;
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
            mean_loss +=  df(mean_ins, ins, D) * df(mean_ins, ins, D);
        }
        delete mean_ins;
        cerr << "c=" << centroids[c] << ", mean_loss: " << mean_loss << endl;
        sum_means_loss += mean_loss;
    }
    // output distance
    cerr << "Total_means_loss: " << sum_means_loss << ", Means_reg: " << means_reg << endl;
    cerr << "Total_Error: " << sum_means_loss + means_reg << endl;
    mat_free (W, N, N);
#else
    // =====================================================
    ifstream fin ("lambdas");
    vector<double> lambda_list;
    vector<double> K_list;  // reg / lambda
    vector<bool> fraction;
    while(!fin.eof()){
        double temp;
        fin >> temp;
        if( fin.eof() )break;
        lambda_list.push_back(temp);	
    }
    fin.close();
    for (int i = 0; i < lambda_list.size(); i ++) { 
        double temp_lambda = lambda_list[i];
        cerr << "D = " << D << endl; // # features
        cerr << "N = " << N << endl; // # instances
        cerr << "lambda = " << temp_lambda << endl;
        cerr << "r = " << r << endl;
        cerr << "Screenshot period = " << screenshot_period << endl;
        int seed = time(NULL);
        srand (seed);
        cerr << "seed = " << seed << endl;
        cerr << "==================================================" << endl; 

        double ** wtemp = mat_init (N, N);
        cvx_clustering (dist_mat, fw_max_iter, D, N, temp_lambda, wtemp, ADMM_max_iter, screenshot_period);
        vector<int> centroids;
        double group_lasso = get_reg (wtemp, N, lambda);
        double reg_labmda_ratio = group_lasso/lambda;
        K_list.push_back(reg_labmda_ratio);
        get_all_centroids (wtemp, &centroids, N, N); 
        int nCentroids = centroids.size();
        if (reg_labmda_ratio > nCentroids - 1e-2 && reg_labmda_ratio < nCentroids + 1e-2) 
            fraction.push_back(true);
        else 
            fraction.push_back(false);
        mat_free(wtemp, N, N);
    }
    ofstream lambdaK_out ("lambda_vs_K_dp");
    assert (lambda_list.size() == K_list.size());
    lambdaK_out << "lambda K fractional" << endl;
    for (int i =0 ; i < K_list.size(); i ++) 
        lambdaK_out << lambda_list[i] << " " <<  K_list[i]
            << " " << (fraction[i]?1:0) <<endl;
    lambdaK_out.close();

    // =====================================================
#endif 
    /* reallocation */
    mat_free (dist_mat, N, N);

}
