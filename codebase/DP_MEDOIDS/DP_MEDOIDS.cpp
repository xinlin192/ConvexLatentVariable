#include "DP_MEDOIDS.h"

ofstream objmin_trace("../obj_vs_time_dp/iris/DP-Medoid");
double objmin = 1e300;
double start_time;

/* Note that entries in dist_mat should be squared distance */
double DP_MEDOIDS (double** dist_mat, int N, double lambda, double** W, vector<int>* medoids) {
    // STEP ZERO: validate input
    for (int i = 0; i < N; i ++) 
        for (int j = 0; j < N; j++) 
            assert (dist_mat[i][j]>=0.0);
    // STEP ONE: a. pick up the optimal *global* medoid for initialization
    vector<int> new_medoids (1,0);
    double * sum_col = new double [N];
    mat_sum_col (dist_mat, sum_col, N, N);
    double min_sum = 1e300;
    int min_sum_index = -1;
    for (int i = 0; i < N; i ++) 
       if (sum_col[i] < min_sum) {
           min_sum = sum_col[i];
           min_sum_index = i;
       }
    delete[] sum_col;
    if (min_sum_index >= 0) 
        new_medoids[0] = min_sum_index;
    else {
        cerr << "Init: optimal global medoids not found!" << endl;
        assert (false); 
    }
    cerr << "optimal global medoid: " << new_medoids[0] << endl;
    ofstream init_out ("init_global_medoid");
    double min_val = INF, max_val = -INF;
    for (int i = 0; i < N; i++) {
        double val = dist_mat[i][new_medoids[0]];
        init_out << val << endl;
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
    }
    init_out << "min: " << min_val << endl;
    init_out << "max: " << max_val << endl;
    init_out.close();
    vector<int> last_medoids (new_medoids);
    // OPTIONAL: write randomized medoids to stdout or external files
    double last_cost = min_sum, new_cost = min_sum; 
    vector<int> z (N, 0);
    while (true) {
        // STEP TWO: compute dist square to new medoids d_ic
        vector< vector<int> > members (new_medoids.size(), vector<int>());
        for (int i = 0; i < N; i ++) {
            int nClusters = new_medoids.size(); 
            vector<double> d_ic (nClusters, 0.0);
            for (int c = 0; c < nClusters; c ++) 
                d_ic[c] = dist_mat[i][new_medoids[c]];
            // STEP THREE: create new medoid if min_c d_ic^2 > lambda
            int min_c = -1;
            double min_dic = INF;
            for (int c = 0; c < nClusters; c ++) 
                if (d_ic[c] < min_dic) {
                    min_c = c;
                    min_dic = d_ic[c];
                }
            if (d_ic[min_c] <= lambda) {
                z[i] = min_c;
                members[min_c].push_back(i);
            } else {
                new_medoids.push_back(i);
                z[i] = nClusters;
                members.push_back(vector<int>(1, i));
            }
        }
        // STEP FOUR: M-step elect new K medoids as PAM algorithm
        int nClusters = new_medoids.size();
        for (int c = 0; c < nClusters; c ++) {
            vector<int> cluster_points = members[c];
            int nPoints = cluster_points.size();
            double ** cluster_dist_mat = mat_init(nPoints, nPoints);
            for (int p = 0; p < nPoints; p ++) {
                for (int q = 0; q < nPoints; q ++) {
                    cluster_dist_mat[p][q] = 
                        dist_mat[cluster_points[p]][cluster_points[q]];
                }
            }
            double* squared_dist = new double [nPoints];
            mat_sum_col(cluster_dist_mat, squared_dist, nPoints, nPoints);
            int min_index = -1;
            double min_value = INF;
            for (int j = 0; j < nPoints; j++) 
                if (squared_dist[j] < min_value) {
                    min_index = j;
                    min_value = squared_dist[j];
                }
            assert (min_index < nPoints);
            new_medoids[c] = cluster_points[min_index];
            delete [] squared_dist;
            mat_free (cluster_dist_mat, nPoints, nPoints);
        }

        // STEP THREE: compute cost
        double loss = 0.0; 
        for (int i = 0; i < N; i ++) 
            loss += 0.5 * dist_mat[i][new_medoids[z[i]]];
        double reg = new_medoids.size() * lambda;
        new_cost = loss + reg;
        cout << "loss: " << loss 
             << ", reg: " << reg 
             << ", new_cost: " << new_cost << endl;
        if ( new_cost < objmin ) objmin = new_cost;
        objmin_trace << omp_get_wtime()-start_time << " " << objmin << endl;
        // STEP FOUR: stopping criteria
        if (new_cost == last_cost) {
            cout << "Resulted medoids: ";
            for (int c = 0; c < nClusters; c++) 
                cout << new_medoids[c] << ",";
            cout << endl;
            break;
        } // medoids has been the optimal
        else if (new_cost < last_cost ) {
            last_medoids.clear();
            last_medoids = new_medoids;
            last_cost = new_cost;
        } else {
            assert (false);
        }
    }
    for (int i = 0; i < N; i ++) 
        for (int j = 0; j < N; j++) 
            if (j == new_medoids[z[i]])
                W[i][j] = 1.0;
            else 
                W[i][j] = 0.0;

    // STEP SEVEN: put converged solution to destination W
    for (int i = 0; i < last_medoids.size(); i ++) 
        medoids->push_back(last_medoids[i]);

    return last_cost;
}

int main (int argc, char ** argv) {
    if (argc != 4) {
        cerr << "Usage: DP_MEDOIDS [dataFile] [nRuns] [lambda]" << endl;
        cerr << "Note: dataFile must be scaled to [0,1] in advance." << endl;
        cerr << "Note: nRuns is the number of running to get global optima" << endl;
        exit(-1);
    }

    // parse arguments
    char* dataFile = argv[1];
    int nRuns = atoi(argv[2]);
    double lambda = atof(argv[3]);
    
    objmin_trace << "time objective" << endl;

    int FIX_DIM;
    Parser parser;
    vector<Instance*>* pdata;
    vector<Instance*> data;
    pdata = parser.parseSVM(dataFile, FIX_DIM);
    data = *pdata;

    // explore the data 
    int D = -1;
    int N = data.size(); // data size
    for (int i = 0; i < N; i++) {
        vector< pair<int,double> > * f = &(data[i]->fea);
        int last_index = f->size() - 1;
        if (f->at(last_index).first > D) {
            D = f->at(last_index).first;
        }
    }
    cerr << "D = " << D << endl; // # features
    cerr << "N = " << N << endl; // # instances
    cerr << "nRuns = " << nRuns << endl;
    cerr << "lambda = " << lambda << endl;
    int seed = time(NULL);
    srand (seed);
    cerr << "seed = " << seed << endl;

    // pre-compute distance matrix
    dist_func df = L2norm;
        random_shuffle (data.begin(), data.end());

    // Run sparse convex clustering
    double *** W = new double** [nRuns];
    vector< vector<int> > medoids (nRuns, vector<int>());
    double* objectives = new double [nRuns];
    double min_obj = INF;
    double ** min_w = mat_init (N, N);
    vector<int> min_medoids;
    int min_nMedoids;
    start_time = omp_get_wtime();
    double ** dist_mat = mat_init (N, N);
    compute_dist_mat (data, dist_mat, N, D, df, true); 
    ofstream dmat_out ("dist_mat");
    dmat_out << mat_toString(dist_mat, N, N);
    dmat_out.close();
    for (int i = 0; i < nRuns; i++) {
        W[i] = mat_init (N, N);
        // INVOKE algorithm function
        objectives[i] = DP_MEDOIDS (dist_mat, N, lambda, W[i], &(medoids[i]));
        if (objectives[i] < min_obj) {
            min_obj = objectives[i];
            mat_copy (W[i], min_w, N, N);
            min_nMedoids = medoids[i].size();
            min_medoids.clear();
            for (int j = 0; j < min_nMedoids; j++) {
                min_medoids.push_back(medoids[i][j]);
            }
        } 
        mat_free(W[i], N, N);
    }

    /* Output objective */
    output_objective(min_obj);
    /* Output model */
    output_model (min_w, N);
    /* Output assignments */
    output_assignment (min_w, data, N);
    
    objmin_trace.close();
    /* Deallocation */
    mat_free (min_w, N, N);
    mat_free (dist_mat, N, N);
}
