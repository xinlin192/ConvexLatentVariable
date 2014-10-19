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
    for (int j = 0; j < nClusters; j++) {
        std::fill(means[j].begin(), means[j].end(), 0.0);
        /*
        for (int i =0; i < means[j].size(); i ++) {
            cout << means[j][i] << " ";
            assert ( means[j][i] == 0.0 );
        }
        cout << endl;
        */
    }
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
            means[c][j] =  means[c][j] / means_count[c];
}

double DP_MEANS (vector<Instance*>& data, int N, int D, double lambda, dist_func df, vector< vector<double> >& means) {
    // STEP ZERO: validate input
    assert (data.size() == N);
    // STEP ONE: a. set initial *global* medoid as global mean
    vector< vector<double> > new_means (1, vector<double>(D, 0.0));
    vector<int> assignment (N, 0);
    vector<int> last_assignment (N, 0);
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
                // d_ic = || x_i - mu_c ||^2
                for (int f = 0; f < D; f++) 
                    mean_ins->fea.push_back(make_pair(f+1, new_means[c][f]));
                double dist_val = df(mean_ins, data[i], D);
                dist_vec[c] = dist_val * dist_val;
                delete mean_ins;
            }
            int min_index = -1;
            double min_value = INF;
            for (int j = 0; j < nClusters; j++) {
                if (dist_vec[j] < min_value) {
                    min_index = j;
                    min_value = dist_vec[j];
                }
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
        double loss = 0.0;
        for (int i = 0; i < N; i ++) {
            double dist_val = df (means_ins[assignment[i]], data[i], D);
            loss +=  dist_val * dist_val;
        }
        double reg = lambda * nClusters;
        new_cost = loss + reg;
        for (int c = 0; c < nClusters; c++) delete means_ins[c];
        cout << "loss: " << loss 
             << ", reg: " << reg 
             << ", new_cost: " << new_cost << endl;
        // STEP FOUR: convergence evaluation
        if (new_cost == last_cost) break;
        else if (new_cost < last_cost) {
            last_cost = new_cost;
            last_means = new_means;
            last_assignment = assignment;
        } else {
            /*
            ofstream la_out ("last_assignment");
            for (int x = 0 ; x < N; x ++)
                la_out << last_assignment[x] <<endl;
            la_out.close();
            ofstream a_out ("assignment");
            for (int x = 0 ; x < N; x ++)
                a_out << assignment[x] <<endl;
            a_out.close();
            ofstream lm_out ("last_means");
            for (int c = 0; c < last_means.size();c++) {
                for (int j = 0; j < D; j++)
                    lm_out << last_means[c][j] << " ";
                lm_out << endl;
            }
            lm_out.close();
            ofstream m_out ("means");
            for (int c = 0; c < new_means.size(); c++) {
                for (int j = 0; j < D; j++)
                    m_out << new_means[c][j] << " ";
                m_out << endl;
            }
            m_out.close();
            */
            assert(false);
        }
    }
    means = last_means;
    return last_cost;
}

// entry main function
int main (int argc, char ** argv) {
    if (argc != 4) {
        cerr << "Usage: DP_MEANS [dataFile] [nRuns] [lambda]" << endl;
        cerr << "Note: dataFile must be scaled to [0,1] in advance." << endl;
        cerr << "Note: nRuns is the number of running to get global optima" << endl;
        exit(-1);
    }
    // parse arguments
    char* dataFile = argv[1];
    int nRuns = atoi(argv[2]);
    double lambda = atof(argv[3]);
    int FIX_DIM;
    Parser parser;
    vector<Instance*>* pdata;
    vector<Instance*> data;
    pdata = parser.parseSVM(dataFile, FIX_DIM);
    data = *pdata;
    // explore the data 
    int D = FIX_DIM;
    int N = data.size(); // data size
    cerr << "D = " << D << endl; // # features
    cerr << "N = " << N << endl; // # instances
    cerr << "lambda = " << lambda << endl;
    int seed = time(NULL);
    srand (seed);
    cerr << "seed = " << seed << endl;
    cerr << "==========================================" << endl;
    // pre-compute distance matrix
    dist_func df = L2norm;
    double min_obj = INF;
    vector<vector<double> > min_means;
    for (int i = 0; i < nRuns; i ++) {
        std::random_shuffle(data.begin(), data.end());
        vector<vector<double> > means;
        double obj = DP_MEANS (data, N, D, lambda, df, means);
        if (obj < min_obj) {
            min_obj = obj;
            min_means = means;
        }
    }
    ofstream obj_out ("opt_objective");
    obj_out << min_obj;
    obj_out.close();

    ofstream model_out ("opt_model");
    model_out << min_means.size() << endl;
    for (int c = 0; c < min_means.size(); c++) {
        for (int j = 0; j < min_means[c].size(); j++)
            model_out << min_means[c][j] << " ";
        model_out << endl;
    }
    model_out.close();
    /* Output objective, model and assignments */
    /*
    output_objective(clustering_objective (dist_mat, min_w, N));
    output_model (min_w, N);
    output_assignment (min_w, data, N);
    */
    /* Deallocation */

}
