#include <vector>
#include<fstream>
#include<cmath>
#include<stdlib.h>

using namespace std;

const int NUM_CHAR = 4;
const char charSet[] = {'A','T','C','G'};

double rand_u01(){
	return ((double)rand()/RAND_MAX);
}
char sample_char(){
	return charSet[rand() % NUM_CHAR];
}

double compute_alpha(double mu){
	
	return exp(-mu);
}

double compute_beta(double lambda, double mu){
	
	return lambda*(1-exp(lambda-mu)) / (mu-lambda*exp(lambda-mu)) ;
}

double compute_gamma(double lambda, double mu){
	
	return 1.0 - mu*(1-exp(lambda-mu)) / (  (1-exp(-mu)) * (mu-lambda*exp(lambda-mu)) ) ;
}

double compute_subProb(double sub_rate){

	return 1-exp(-sub_rate);
}

double logLike_align( vector< vector<char> >& des_seq_aligned, double ins_rate, double del_rate ){
	
	double alpha = compute_alpha(del_rate);
	double log_alpha = log(alpha);
	double log_1_alpha = log(1-alpha);
	
	double beta = compute_beta(ins_rate,del_rate);
	double log_beta = log(beta);
	double log_1_beta = log(1-beta);

	double gamma = compute_gamma(ins_rate,del_rate);
	double log_gamma = log(gamma);
	double log_1_gamma = log(1-gamma);
	
	double log_ins_del = log(ins_rate/del_rate);
	
	double logLike = des_seq_aligned[0].size()*log(beta) + log(1-beta); //loglike of initial inserts
	for(int i=1;i<des_seq_aligned.size();i++){
		
		logLike += log_ins_del;

		//logLike += log_ins_del;
		vector<char>* frag = &(des_seq_aligned[i]);
		//whether the taxon survived
		if( frag->at(0) == '-' )
			logLike += log_1_alpha;
		else
			logLike += log_alpha;
		//whether exists insertion
		if( frag->size() > 1 ){ 
			logLike += (frag->at(0)=='-')? log_gamma : log_beta ;
			//additional insertions
			if( frag->size() > 2 )
				logLike += (frag->size()-2) * log_beta;
			//stop insert
			logLike += log_1_beta;
		}else{
			//no insert
			logLike += (frag->at(0)=='-')? log_1_gamma : log_1_beta ;
		}
		
	}
	logLike += log(1-ins_rate/del_rate);

	return logLike;
}

void readPairSeqs(char* fname, vector<vector<char> >&  anc_seqs, vector< vector<vector<char> > >& seqs_aligned){
	
	ifstream fin(fname);
	int N;
	fin >> N;
	fin.get();
	
	char c;
	for(int n=0;n<N;n++){
		
		//read ancestral sequence
		vector<char> anc_seq;
		while(1){
			c = fin.get();
			if( c == '\n' ) break;
			
			if( c != ' ' && c!='|' ){
				anc_seq.push_back(c);
			}
		}
		
		//read descendant sequence
		vector<vector<char> > des_seq;
		des_seq.push_back(vector<char>());
		int i = 0;
		while(1){
			c = fin.get();
			if( c == '\n' ) break;
			if( c == '|' ){
			       	i++;
				des_seq.push_back(vector<char>());
				continue;
			}
			
			des_seq[i].push_back(c);
		}
		des_seq.pop_back();

		anc_seqs.push_back(anc_seq);
		seqs_aligned.push_back(des_seq);
	}
	
	fin.close();
}
