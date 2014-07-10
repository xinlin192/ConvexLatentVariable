#include<fstream>
#include<vector>
#include<iostream>
#include"TKF1.h"

using namespace std;

const bool print_align = 1;

int main(int argc, char** argv){
	
	if( argc < 5 ){
		
		cerr << "./TKF1 [N] [length] [sub_rate] [del_rate] (output)" << endl;
		exit(0);
	}
	
	srand(time(NULL));
	//srand(1);
	int N = atoi(argv[1]);
	int L = atoi(argv[2]);
	double sub_rate = atof(argv[3]);
	double del_rate = atof(argv[4]);
	double ins_rate = del_rate*L/(1+L);
	
	double alpha = compute_alpha(del_rate);
	double beta = compute_beta(ins_rate,del_rate);
	double gamma = compute_gamma(ins_rate,del_rate);
	double subProb = compute_subProb(sub_rate);
	
	cout << "sub=" << sub_rate << ", ins=" << ins_rate << ", del=" << del_rate << endl;
	cout << "P(taxon survive)=" << alpha << endl;
	cout << "P(insert|taxon not die)=" << beta << endl;
        cout << "P(insert|taxon die)=" << gamma << endl;
	cout << "P(subsitution)=" << subProb << endl;

	char* outputFile;
	if( argc > 5 ){
		outputFile = argv[5];
	}else{
		outputFile = "output";
	}
	
	ofstream fout(outputFile);
	fout << N << endl ;
	//ancestral sequence (geometric with p = ins_rate/del_rate)
	vector<char> anc_seq;
	for(int i=0;i<L;i++)
		anc_seq.push_back(sample_char());
	
	for(int n=0;n<N;n++){
			
		//descent sequence 
		vector<vector<char> > des_seq;
		//insert before 1st match/delete
		des_seq.push_back(vector<char>());
		while(rand_u01() < beta)
			des_seq[0].push_back(sample_char());
		
		for(int i=0;i<anc_seq.size();i++){
			des_seq.push_back(vector<char>());
			//evolution for the i-th taxon of ancestral seq
			double has_insert_prob;
			if( rand_u01() < alpha ){ //taxon survive
				if(rand_u01()<subProb) des_seq[i+1].push_back(sample_char());
				else			des_seq[i+1].push_back(anc_seq[i]);
				has_insert_prob = beta;
			}else{
				des_seq[i+1].push_back('-');
				has_insert_prob = gamma;
			}
			
			if( rand_u01() < has_insert_prob ){
				
				des_seq[i+1].push_back(sample_char());
				while(rand_u01() < beta)
					des_seq[i+1].push_back(sample_char());
			}
		}
		
		if( print_align ){
			//print on buffer
			vector<char> line1;
			vector<char> line2;
			//initial insertion
			for(int j=0;j<des_seq[0].size();j++){
				line1.push_back(' ');
				line2.push_back(des_seq[0][j]);
			}
			line1.push_back('|');
			line2.push_back('|');

			//for each ancestral taxon
			for(int i=0;i<anc_seq.size();i++){
				line1.push_back(anc_seq[i]);
				line2.push_back(des_seq[i+1][0]);
				for(int j=1;j<des_seq[i+1].size();j++){
					line1.push_back(' ');
					line2.push_back(des_seq[i+1][j]);
				}
				line1.push_back('|');
				line2.push_back('|');
			}

			//pinrt to file
			for(int i=0;i<line1.size();i++)
				fout << line1[i];
			fout << endl;
			for(int i=0;i<line2.size();i++)
				fout << line2[i];
			fout << endl;
		}else{

			//directly print to file
			for(int i=0;i<anc_seq.size();i++)
				fout << anc_seq[i] ;
			fout << endl;
			for(int i=0;i<des_seq.size();i++){
				for(int j=0;j<des_seq[i].size();j++)
					if(des_seq[i][j] != '-')
						fout << des_seq[i][j];
			}
			fout << endl;
		}
	}
	fout.close();
}
