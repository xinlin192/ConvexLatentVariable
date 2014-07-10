#include <cstring>
#include <iostream>
#include <fstream>
#include "HMMdec.h"
using namespace std;

void readPairSeqs_unaligned(char* dataFile, vector<vector<char> >& anc_seqs, vector<vector<char> >& des_seqs){

	ifstream fin(dataFile);
	char s[LARGE_NUM];
	int N;
	fin >> N;
	fin.getline(s,LARGE_NUM);
	
	for(int n=0;n<N;n++){

		anc_seqs.push_back(vector<char>());
		des_seqs.push_back(vector<char>());
		
		//anc_seq
		fin.getline(s,LARGE_NUM);
		for(int j=0;j<LARGE_NUM && s[j]!='\0' ;j++)
			anc_seqs[n].push_back(s[j]);
		
		//des_seq
		fin.getline(s,LARGE_NUM);
		for(int j=0;j<LARGE_NUM && s[j]!='\0' ;j++)
			des_seqs[n].push_back(s[j]);
	}
	
	fin.close();	
}

double estimate_TKF1(vector<vector<char> >& anc_seqs, vector<vector<char> >& des_seqs, 
			double& sub_rate, double& ins_rate, double& del_rate ){
	
	int N = anc_seqs.size();

	//Grid Search Parameters
	double prec = 0.01;
	double x_min = 0.15, x_max = 0.3; //sub_rate
	double y_min = 0.15, y_max = 0.25; //ins_rate
	double z_min = 0.15, z_max = 0.25; //del_rate
	
	int I = (z_max-z_min)/prec+1 ;
	int J = (y_max-y_min)/prec+1 ;
	int K = (x_max-x_min)/prec+1 ;
	
	//gird matrix
	double** matrix = new double*[I];
	for(int i=0;i<I;i++){
		matrix[i] = new double[J];
		for(int j=0;j<J;j++)
			matrix[i][j] = 0.0;
	}
	
	//Grid Search
	double ll_max = -1e300;
	//#pragma omp parallel for
	for(int i=0;i<I;i++){
		double z = z_min + i*prec;
		for(int j=0;j<J;j++){
			
			double y = y_min + j*prec;
			if( y >= z ){ //ins_rate >= del_rate
				continue;
			}
			
			for(int k=0;k<K;k++){
				double x = x_min + k*prec;
				//double x = 0.2;
				cerr << x << ", " << y << ", " << z << endl;
				//compute marginal data logLike
				PairHMMState::setupPairHMM(x,y,z); //sub,ins,del
				double logLike = 0.0;
				double ll;
				//#pragma omp parallel for
				for(int n=0;n<N;n++){
					
					ll = forward_backward(des_seqs[n], anc_seqs[n], new PairHMMState(State::BEGIN,0,0));
					//#pragma omp critical
					{
						logLike += ll;
					}
				}
				logLike /= N;
				
				matrix[i][j] = logLike;
				
				if( logLike > ll_max ){
					ll_max = logLike;
					sub_rate = x;
					ins_rate = y;
					del_rate = z;
				}
			}
		}
	}
	
	//print matrix to "grid"
	/*char str[1000];
	sprintf(str,"N=%d",N);
	ofstream fout(str);
	for(int i=0;i<I;i++){
		for(int j=0;j<J;j++){
			if( matrix[i][j] != 0.0 )
				fout << matrix[i][j] << " ";
			else
				fout << "nan ";
		}
		fout << endl;
	}
	fout.close();*/
	
	//delete matrix
	for(int i=0;i<I;i++)
		delete matrix[i];
	delete matrix;
	
	return ll_max;
}

int main(int argc, char** argv){
	
	if( argc < 2 ){
		cerr << "./est_TKF1_HMM  [data_unaligned] (model)" << endl;
		exit(0);
	}

	char* dataFile = argv[1];
	char* outputFile;
	if( argc > 2 )
		outputFile = argv[2];
	else
		outputFile = "model";
	
	//read paired sequences
	vector<vector<char> > anc_seqs;
	vector<vector<char> > des_seqs;
	readPairSeqs_unaligned(dataFile, anc_seqs, des_seqs);
	int N = anc_seqs.size();
	/*
	cerr << "anc:" ;
	for(int i=0;i<anc_seqs[613].size();i++)
		cerr << anc_seqs[0][i] ;
	cerr << endl;
	cerr << "des:" ;
	for(int i=0;i<des_seqs[613].size();i++)
		cerr << des_seqs[0][i];
	cerr << endl;
	*/
	
	//esitimation
	double sub_rate, ins_rate, del_rate;
	double ll_max = estimate_TKF1(anc_seqs, des_seqs, sub_rate, ins_rate, del_rate);
	
	cerr << "logLike=" << ll_max << endl;
	cerr << "s=" << sub_rate << ", lambda=" << ins_rate << ", mu=" << del_rate << endl;
	
	//compute ground truth loglike
	/*PairHMMState::setupPairHMM(0.2,0.416667,0.5); //sub,ins,del
	double logLike = 0.0;
	double ll;
	for(int n=0;n<N;n++){

		ll = forward_backward(des_seqs[n], anc_seqs[n], new PairHMMState(State::BEGIN,0,0));
		logLike += ll;
	}
	logLike /= N;
	cerr << "loglike(ground)=" << logLike << endl;*/

	return 0;
}
