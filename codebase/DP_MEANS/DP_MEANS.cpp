#include "DP_MEANS.h"

/* 
    data: N data instance
    assignment: N by 1 vector indicating assignment of one atom
*/
void compute_means (vector<Instance*>& data, vector<int>& assignment, int D, vector< vector<double> >& means) {
    int N = data.size();
    assert (assignment.size() == N);
    // compute number of existing clusters
    int nClusters = means.size();  // for 0-index
    // clear means
    for (int j = 0; j < nClusters; j++) 
        std::fill(means[j].begin(), means[j].end(), 0.0);
    // compute sums and counts
    vector<int> means_count (nClusters, 0);
    for (int i = 0; i < N; i++) {
        int belongsto = assignment[i];
        for (int f = 0; f < data[i]->fea.size(); f++) {
            int j = data[i]->fea[f].first - 1;
            means[belongsto][j] += data[i]->fea[f].second;
        }
        ++ means_count[belongsto];
    }
    // compute means
    for (int c = 0; c < nClusters; c ++) 
        for (int j = 0; j < D; j++) 
            means[c][j] = means[c][j] / means_count[c];
}

double DP_MEDOIDS (vector<Instance*>& data, int N, int D, double lambda, dist_func df) {
    // STEP ZERO: validate input
    assert (data.size() == N);
    // STEP ONE: a. set initial *global* medoid as global mean
    vector< vector<double> > new_means (1, vector<double>(D, 0.0));
    vector<int> assignment (N, 0);
    compute_means (data, assignment, D, new_means);
    vector< vector<double> > last_means = new_means; 
    double last_cost = INF, new_cost = INF; // compute last cost
    while (true) {
        // STEP TWO: compute dist square to new medoids d_ic
        for (int i = 0; i < N; i ++) {
            int nClusters = new_means.size();
            // cout << "nClusters: " << nClusters << endl;
            vector<double> dist_vec (nClusters, 0.0);
            for (int c = 0; c < nClusters; c ++) {
                Instance* mean_ins = new Instance (10000+c);
                for (int f = 0; f < D; f++) 
                    mean_ins->fea.push_back(make_pair(f+1, new_means[c][f]));
                dist_vec[c] = df(mean_ins, data[i], D);
                delete mean_ins;
            }
            int min_index = -1;
            double min_value = INF;
            for (int j = 0; j < nClusters; j++) 
                if (dist_vec[j] < min_value) {
                    min_index = j;
                    min_value = dist_vec[j];
                }
            // cout << "min_value: " << min_value << endl;
            if (min_value <= lambda) 
                assignment[i] = min_index;
            else {
                assignment[i] = nClusters;
                vector<double> tmp_mean (D, 0.0);
                int F = data[i]->fea.size();
                for (int f = 0; f < F; f ++) 
                    tmp_mean[data[i]->fea[f].first-1] = data[i]->fea[f].second;
                new_means.push_back(tmp_mean);
            } 
        }
        compute_means (data, assignment, D, new_means);
        // STEP THREE: compute cost
        int nClusters = new_means.size();
        vector<Instance*> means_ins (nClusters, NULL);
        for (int c = 0; c < nClusters; c++) {
            means_ins[c] = new Instance (5000+c);
            for (int j = 0; j < D; j ++) 
                means_ins[c]->fea.push_back(make_pair(j+1, new_means[c][j]));
        }
        new_cost = 0.0;
        for (int i = 0; i < N; i ++) 
            new_cost += 0.5 * df (means_ins[assignment[i]], data[i], D);
        new_cost += lambda * nClusters;
        for (int c = 0; c < nClusters; c++) delete means_ins[c];
        cout << "CLUSTERING COST: " << new_cost << endl;
        cout << "c: " << nClusters << endl;
        // STEP FOUR: convergence evaluation
        if (new_cost != last_cost) {
            last_cost = new_cost;
            last_means = new_means;
        } else break;
    }
    return last_cost;
}

// entry main function
int main (int argc, char ** argv) {
    if (argc != 3) {
        cerr << "Usage: DP_MEDOIDS [dataFile] [lambda]" << endl;
        cerr << "Note: dataFile must be scaled to [0,1] in advance." << endl;
        cerr << "Note: nRuns is the number of running to get global optima" << endl;
        exit(-1);
    }
    // parse arguments
    char* dataFile = argv[1];
    double lambda = atof(argv[2]);
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
        if (f->at(last_index).first > D) 
            D = f->at(last_index).first;
    }
    cerr << "D = " << D << endl; // # features
    cerr << "N = " << N << endl; // # instances
    cerr << "lambda = " << lambda << endl;
    int seed = time(NULL);
    srand (seed);
    cerr << "seed = " << seed << endl;
    cerr << "==========================================" << endl;
    // pre-compute distance matrix
    dist_func df = L2norm;
    DP_MEDOIDS (data, N, D, lambda, df);

    /* Output objective, model and assignments */
    /*
    output_objective(clustering_objective (dist_mat, min_w, N));
    output_model (min_w, N);
    output_assignment (min_w, data, N);
    */
    /* Deallocation */

}
