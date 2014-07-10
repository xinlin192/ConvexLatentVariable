#include "TKF1.h"
#include <cstring>
#include <omp.h>
#include <vector>
#include <iostream>
using namespace std;

void estimate_sub_rate(vector<vector<char> >& anc_seqs, vector<vector<vector<char> > >& seqs_aligned,
			double& sub_rate ){
	
	double subProb; // estimate subProb = \frac{4}{3} * \frac{F}{F+T}
	int mis_match = 0;
	int match = 0;
	int N = seqs_aligned.size();

	char c;
	for(int n=0;n<N;n++){
		for(int i=1;i<seqs_aligned[n].size();i++){
			c = seqs_aligned[n][i][0];
			if( c == '-' )
				continue;

			if( c == anc_seqs[n][i-1] )
				match++;
			else
				mis_match++;
		}
	}

	subProb = (double)mis_match/(mis_match + match) * 4 / 3;
	sub_rate = -log(1-subProb);
}

double logLike_align_avg(vector<vector<vector<char> > >& seqs_aligned, double ins, double del){
	
	int N = seqs_aligned.size();
	double logLike = 0.0;
	for(int n=0;n<N;n++)	
		logLike += logLike_align(seqs_aligned[n], ins, del);
	logLike = logLike / N;

	return logLike;
}

double estimate_indel_rate(vector<vector<vector<char> > >& seqs_aligned,
			double& ins_rate, double& del_rate ){

	int N = seqs_aligned.size();

	//Grid Search parameters
	double prec = 0.003;
	double x_min = 0.1, x_max=0.3;
	double y_min = 0.1, y_max=0.3;
	int I = (y_max-y_min)/prec +1 ;
	int J = (x_max-x_min)/prec +1;
	
	//Grid Search
	double ll_max = -1e300;
	
	double** matrix = new double*[I];
	for(int i=0;i<I;i++){
		matrix[i] = new double[J];
		for(int j=0;j<J;j++){
			matrix[i][j] = 0.0;
		}
	}
	
	#pragma omp parallel for
	for(int i=0;i<I;i++){
		double y = y_min + i*prec;
		for(int j=0;j<J;j++){
			
			double x = x_min + j*prec;

			if( x >= y ){
				continue;
			}

			//compute Loglike at (ins_rate=x, del_rate=y)
			double logLike = logLike_align_avg(seqs_aligned, x, y);
			
			matrix[i][j] = logLike;
			//find max
			if( logLike > ll_max ){
				ll_max = logLike;
				ins_rate = x;
				del_rate = y;
			}
		}
	}

	//Write Grid Matrix to "grid"
	char s[1000];
	sprintf(s,"logLike(insert-rate,deletion-rate),N=%d",N);
	ofstream fout(s);
	for(int i=0;i<I;i++){
		for(int j=0;j<J;j++){
			if(matrix[i][j] != 0.0)
				fout << matrix[i][j] << " ";
			else
				fout << "nan ";
		}
		fout << endl;
	}
	fout.close();

	return ll_max;
}


int main(int argc, char** argv){

	if( argc < 2 ){
		cerr << "./estimate_TKF1 [data_aligned] (model)" << endl;
		exit(0);
	}
	char* dataFile = argv[1];
	char* outputFile;
	if( argc > 2 )
		outputFile = argv[2];
	else
		outputFile = "model";

	//read paired sequence
	vector<vector<vector<char> > > seqs_aligned;
	vector<vector<char> > anc_seqs;
	readPairSeqs(dataFile, anc_seqs, seqs_aligned);
	
	double sub_rate,ins_rate,del_rate;
	estimate_sub_rate(anc_seqs, seqs_aligned, sub_rate);
	double logLike = estimate_indel_rate(seqs_aligned, ins_rate, del_rate);
	cout << "align_loglike=" << logLike << endl;
	cout << "s=" << sub_rate << ", lambda=" << ins_rate << ", mu=" << del_rate << endl;

	ofstream fout(outputFile);
	fout << "sub_rate " << sub_rate << endl;
	fout << "ins_rate " << ins_rate << endl;
	fout << "del_rate " << del_rate << endl;
	fout.close();
}
