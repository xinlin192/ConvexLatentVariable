ofstream objmin_trace("../obj_vs_time_dp/iris/DP-Medoid");
double objmin = 1e300;
double start_time;

/* Compute the mutual distance of input instances contained within "data" */
void compute_dist_mat (vector<Instance*>& data, double ** dist_mat, int N, int D, dist_func df, bool isSym) {
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            Instance * xi = data[i];
            Instance * muj = data[j];
            dist_mat[i][j] = df (xi, muj, D);
        }
    }
}

/* Note that entries in dist_mat should be squared distance */
double DP_MEDOIDS (double** dist_mat, int N, double lambda, double** W, vector<int>* medoids) {
    // STEP ZERO: validate input
    for (int i = 0; i < N; i ++) 
        for (int j = 0; j < N; j++) 
            assert (dist_mat[i][j]>=0.0);
    // STEP ONE: a. pick up the optimal *global* medoid for initialization
    vector<int> new_medoids (1,0);
    // new_medoids[0] = random() % N; // by random
    // cout << "randomized medoids: " << new_medoids[0]  << endl;
    double * sum_col = new double [N];
    mat_sum_col (W, sum_col, N, N);
    double max_sum = -1e300;
    int max_sum_index = -1;
    for (int i = 0; i < N; i ++) {
       if (sum_col[i] > max_sum) {
           max_sum = sum_col[i];
           max_sum_index = i;
       }
    }
    if (max_sum_index >= 0) 
        new_medoids[0] = max_sum_index;
    else {
        cerr << "Init: optimal global medoids not found!" << endl;
        assert (false); 
    }
    delete[] sum_col;
    cout << "optimal global medoid: " << new_medoids[0] << endl;
    vector<int> last_medoids (new_medoids);
    // OPTIONAL: write randomized medoids to stdout or external files
    double last_cost = INF, new_cost = INF; // compute last cost
    double ** w = mat_init (N,N);
    mat_zeros(w, N,N);
    for (int i = 0; i < N; i ++) {
        w[i][last_medoids[0]] = 1.0;
    }
    while (true) {
        // STEP TWO: compute dist square to new medoids d_ic
        int last_K = last_medoids.size();
        int K = new_medoids.size(); 
        double ** temp = mat_init (N,K); // temp[i][c] <- d_ic
        for (int i = 0; i < K; i ++) {
            for (int j = 0; j < N; j ++) {
                temp[j][i] = dist_mat[j][new_medoids[i]];
            }
        }
        // STEP THREE: create new cluster medoid if min_c d_ic > lambda
        // for here, we make the condition squared:
        //           min_c d_ic^2 > lambda^2
        double* min_d_ic = new double [N];
        mat_min_row (temp, min_d_ic, N,K);
        bool cont = false;
        vector<int> candiates;
        for (int i = 0; i < N; i ++) {
            if (min_d_ic[i] > lambda) {
                candiates.push_back(i); // add to candidate list
            }
        }
        int nCandidates = candiates.size();
        if (nCandidates > 0) {
            mat_free(temp, N,K);
            int new_med = rand() % nCandidates;
            new_medoids.push_back(candiates[new_med]);
            K = new_medoids.size();
            temp = mat_init (N,K);
            for (int i = 0; i < K; i ++) {
                for (int j = 0; j < N; j ++) {
                    temp[j][i] = dist_mat[j][new_medoids[i]];
                }
            }
        }
        double* min_index = new double [N];
        mat_min_index_row (temp, min_index, N,K);
        mat_zeros(w, N,N);
        for (int i = 0; i < N; i ++) {
            // reassign to one of new medoids
            w[i][new_medoids[min_index[i]]] = 1.0;
        }
        delete[] min_index;
        delete[] min_d_ic;
        mat_free(temp, N,K);

        // STEP THREE: compute cost
        new_cost = 0.5 * mat_frob_dot (w, dist_mat, N,N)+ K * lambda ;

        cout << "new_cost: " << new_cost << endl;
        if( new_cost < objmin ){
            objmin = new_cost;
        }
        objmin_trace << omp_get_wtime()-start_time << " " << objmin << endl;
        

        // STEP FOUR: stopping criteria
        if (new_cost >= last_cost && last_K == K) {
            cout << "CLUSTERING COST: " << last_cost << endl;
            cout << "Resulted medoids: ";
            for (int i = 0; i < last_K; i++) {
                cout << last_medoids[i] << ",";
            }
            cout << endl;
            break;
        } // medoids has been the optimal
        else {
            last_medoids.clear();
            for (int i = 0; i < K; i ++) {
                last_medoids.push_back(new_medoids[i]);
            }
            last_K = K;
            last_cost = new_cost;
        }
        // STEP FIVE: M-step elect new K medoids as PAM algorithm
        for (int i = 0; i < K; i ++) {
            vector<int> cluster_points;
            for (int j = 0; j < N; j ++) {
                if (w[j][last_medoids[i]] > 0) 
                    cluster_points.push_back(j);
            }
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
            int min_index = INTEGER_MAX;
            double min_value = INF;
            for (int j = 0; j < nPoints; j++) {
                if (squared_dist[j] < min_value) {
                    min_index = j;
                    min_value = squared_dist[j];
                }
            }
            assert (min_index < nPoints);
            new_medoids[i] = cluster_points[min_index];
            delete [] squared_dist;
            mat_free (cluster_dist_mat, nPoints, nPoints);
        }
    }

    // STEP SEVEN: put converged solution to destination W
    mat_copy (w, W, N, N);
    for (int i = 0; i < last_medoids.size(); i ++) 
        medoids->push_back(last_medoids[i]);
    // STEP EIGHT: memory recollection
    mat_free (w, N, N);

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
    double ** dist_mat = mat_init (N, N);
    mat_zeros (dist_mat, N, N);
    compute_dist_mat (data, dist_mat, N, D, df, true); 

    // Run sparse convex clustering
    double *** W = new double** [nRuns];
    vector< vector<int> > medoids (nRuns, vector<int>());
    double* objectives = new double [nRuns];
    double min_obj = INF;
    double ** min_w = mat_init (N, N);
    vector<int> min_medoids;
    int min_nMedoids;
    start_time = omp_get_wtime();
    for (int i = 0; i < nRuns; i++) {
        W[i] = mat_init (N, N);
        mat_zeros (W[i], N, N);
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
    }

    /* Output objective */
    output_objective(clustering_objective (dist_mat, min_w, N));
    /* Output model */
    output_model (min_w, N);
    /* Output assignments */
    output_assignment (min_w, data, N);
    
    objmin_trace.close();
    /* Deallocation */
    mat_free (min_w, N, N);
    for (int i = 0; i < nRuns; i++) 
        mat_free (W[i], N, N);
    mat_free (dist_mat, N, N);
}
