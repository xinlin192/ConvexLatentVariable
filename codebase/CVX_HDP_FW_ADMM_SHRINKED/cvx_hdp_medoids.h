#include "../util.h"
#include <time.h>
#include <queue>
#include <cassert>
#include <cstring>

const int FNAME_LEN = 1000;
typedef struct {
	int nWords;
	int nVocs;
	int nDocs;
	vector< pair<int,int> >* doc_lookup; // document start,end word id pair
	vector< pair<int,int> >* word_lookup;  // word id to voc id
	vector< vector<int> >* voc_lookup;  // voc id to word id set
	vector< pair<int,int> >* word_in_doc;
} Lookups ;

class Int_Double_Pair_Dec
{
	public:
		bool operator() (pair<int, double> obj1, pair<int, double> obj2)
		{
			return obj1.second > obj2.second;
		}
};

double lasso_objective (double** z, double lambda, int R_start, int R_end, int C) {
	double lasso = 0.0;
	vector<double> maxn(C, -INF); 
	for (int i = R_start; i < R_end; i ++)
		for (int j = 0; j < C; j ++)
			if ( fabs(z[i][j]) > maxn[j])
				maxn[j] = fabs(z[i][j]);
	for (int j = 0; j < C; j ++) 
		lasso += lambda * maxn[j];
	return lasso;
}

void output_objective (double** dist_mat, double ** W, Lookups * tables, double dummy_rate, vector<double>& lambda) {

	
	char fname[FNAME_LEN];
	sprintf(fname,"lambda_theta_tune_syn/obj-l%.2f-t%.2f",lambda[0],lambda[1]);
	ofstream obj_out (fname);
	int N = tables->nWords;
	int D = tables->nDocs;
	vector< pair<int,int> > doc_lookup = *(tables->doc_lookup);
	double normSum = 0.0;
	for (int i = 0; i < N; i ++) 
		for (int j = 0; j < N; j ++) 
			normSum += W[i][j] * dist_mat[i][j];
	double loss = 0.5 * normSum;
	// STEP TWO: compute dummy loss
	// sum4 = r dot (1 - sum_k w_nk) -> dummy
	double * temp_vec = new double [N];
	mat_sum_row (W, temp_vec, N, N);
	double dummy_penalty = 0.0;
	for (int i = 0; i < N; i ++) 
		dummy_penalty += dummy_rate * max(1 - temp_vec[i], 0.0);
	delete[] temp_vec;
	// STEP THREE: compute group-lasso regularization
	double global_lasso = lasso_objective(W, lambda[0], 0, N, N);
	double sum_local_lasso = 0.0;
	vector<double> local_lasso (D, 0.0);
	for (int d = 0; d < D; d++) {
		local_lasso[d] = lasso_objective(W, lambda[1], doc_lookup[d].first, doc_lookup[d].second, N); 
		sum_local_lasso += local_lasso[d];
	}
	double overall = loss + global_lasso + sum_local_lasso + dummy_penalty;
	obj_out << "loss: " << loss << endl;
	obj_out << "dummy: " << dummy_penalty << endl;
	obj_out << "global_lasso: " << global_lasso << endl;
	obj_out << "sum_local_lasso: " << sum_local_lasso << endl;
	obj_out << "loss+global+local: " <<loss+global_lasso +sum_local_lasso<< endl; 
	obj_out << "overall: " << overall << endl;
	obj_out << "======================================" << endl;
	obj_out << "local lasso dump for each dataset" << endl;
	for (int d = 0; d < D; d++) 
		obj_out << "dataset(" << d << "): " << local_lasso[d] << endl;
	obj_out.close();
}

void output_model (double ** W, Lookups * tables, vector<double>& lambda) {
	int N = tables->nWords;
	int D = tables->nDocs;
	vector< pair<int,int> > doc_lookup = *(tables->doc_lookup); 


	char fname[FNAME_LEN];
	sprintf(fname,"lambda_theta_tune_syn/model-l%.2f-t%.2f",lambda[0],lambda[1]);
	ofstream model_out (fname);

	vector<int> centroids;
	get_all_centroids (W, &centroids, N, N); // contains index of all centroids
	int nCentroids = centroids.size();
	model_out << "nCentroids: " << nCentroids << endl;
	for (int i = 0; i < nCentroids; i ++) 
		model_out << "centroids[" << i <<"]: " << centroids[i] << endl;

	for (int d = 0; d < D; d++) {
		model_out << "======================================" << endl;
		model_out << "d = " << d << endl;
		int nlocalwords =  doc_lookup[d].second-doc_lookup[d].first;
		double ** wtemp = mat_init (nlocalwords, N);
		for (int i = doc_lookup[d].first; i < doc_lookup[d].second; i ++) {
			for (int j = 0; j < N; j ++) {
				wtemp[i-doc_lookup[d].first][j] = W[i][j];
			}
		}

		vector<int> centroids;
		get_all_centroids (wtemp, &centroids, nlocalwords, N); // contains index of all centroids
		int nCentroids = centroids.size();
		model_out << "nCentroids: " << nCentroids << endl;
		for (int i = 0; i < nCentroids; i ++) 
			model_out << "centroids[" << i <<"]: " << centroids[i] << endl;

		mat_free(wtemp, nlocalwords, N);
	}
	model_out.close();
}
void output_assignment (double ** W, Lookups * tables, vector<double>& lambda ) {
	int N = tables->nWords;
	int D = tables->nDocs;
	vector< pair<int,int> > doc_lookup = *(tables->doc_lookup); 
	
	char fname[FNAME_LEN];
	sprintf(fname,"lambda_theta_tune_syn/assign-l%.2f-t%.2f",lambda[0],lambda[1]);
	ofstream asgn_out (fname);
	// get all centroids
	for (int d = 0; d < D; d++) {
		asgn_out << "d = " << d << endl;
		for (int i = doc_lookup[d].first; i < doc_lookup[d].second; i ++) {
			// output identification and its belonging
			asgn_out << "  id=" << i << ", "; 
			for (int j = 0; j < N; j ++) 
				if ( fabs(W[i][j]) > 1e-2 ) 
					asgn_out << j << "(" << W[i][j] << "), ";
			asgn_out << endl;
		}
	}
	asgn_out.close();
}

void get_doc_lookup (vector< Instance* > & data, vector<pair<int,int> >& doc_lookup) {
	// validate the input
	doc_lookup.clear();
	int N = data.size();
	for (int i = 1; i < N; i ++) 
		assert(data[i-1]->label <= data[i]->label);
	// compute list of document info
	int doc_begin = 0;
	int doc_end = 0;
	int last_doc = data[0]->label;
	for (int i = 0; i < N ; i++) {
		int curr_doc = data[i]->label;
		if (curr_doc != last_doc) {
			doc_end = i;
			cerr << "(" << doc_begin << ", " << doc_end << ")" <<endl;
			doc_lookup.push_back(make_pair(doc_begin, doc_end));
			doc_begin = i;
			last_doc = curr_doc;
		} 
	}
	cerr << "(" << doc_begin << ", " << N << ")" <<endl;
	doc_lookup.push_back(make_pair(doc_begin, N));
}
