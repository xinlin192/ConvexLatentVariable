#include <cmath>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "State.h"

using namespace std;
int LARGE_NUM = 10000;

void readData( char* file_seq, vector<vector<char> >& data){
	
	ifstream fin(file_seq);
	char s[LARGE_NUM];
	int i=0;
	while( !fin.eof() ){

		data.push_back(vector<char>());
		
		fin.getline(s,LARGE_NUM);
		for(int j=0;j<LARGE_NUM && s[j]!='\0' ;j++)
			data[i].push_back(s[j]);
		
		if( data[i].size() < 1 && fin.eof() ){
			data.pop_back();
			break;
		}

		i++;
	}
	
	fin.close();
}

void readRefSeq( char* file_ref, vector<char>& ref){
	
	
	ifstream fin(file_ref);
	char c;
	int i=0;
	while( 1 ){
		fin >> c;
		if( fin.eof() )break;
		ref.push_back(c);
	}
	
	fin.close();
}


double align_Viterbi( vector<char>& seq, vector<char>& ref, State* start_state, vector<State*>& output ){
	
	vector< SparseDist > max_marg; //state(t) -> log-proby
	vector< SparseTable > max_arg; //state(t+1) -> state(t)
	
	//initialize
	vector<pair<State*,double> > init_marg = start_state->transLogProb(seq,ref);
	max_marg.push_back(SparseDist());
	max_arg.push_back(SparseTable());
	////initial transition probability
	for(int i=0;i<init_marg.size();i++){
		max_marg[0].insert( init_marg[i] );
		max_arg[0].insert( make_pair( init_marg[i].first, start_state ) );
	}
	
	//forward pass
	int t=0;
	while(1){
		
		////times eimission probability
		SparseDist::iterator it;
		vector<State*> del_vec;
		for(it=max_marg[t].begin(); it!=max_marg[t].end(); it++){

			State* s_t = it->first;
			it->second += s_t->emitLogProb(seq,ref);
			if( it->second < -1e10 )
				del_vec.push_back(it->first);
		}
		for(int i=0;i<del_vec.size();i++)
			max_marg[t].erase(del_vec[i]);
		
		
		//compuat max_marg[t+1]
		max_marg.push_back(SparseDist());
		max_arg.push_back(SparseTable());
		
		for(it=max_marg[t].begin();it!=max_marg[t].end();it++){ //for all current states
			
			State* curState = it->first;
			//curState->print(cerr);
			//cerr << ":" << it->second << ", ";
			double logMarg_t = it->second;
			vector<pair<State*,double> > transVec = curState->transLogProb(seq,ref);
			
			for(int j=0;j<transVec.size();j++){ // transit to their next states
				 
				State* nextState = transVec[j].first;
				double nextLogProb = logMarg_t + transVec[j].second;
				
				SparseDist::iterator it2;
				if( (it2=max_marg[t+1].find( nextState )) == max_marg[t+1].end()  ){
					
					max_marg[t+1].insert( make_pair(nextState, nextLogProb) );
					max_arg[t+1].insert( make_pair(nextState, curState) );
				}
				else if( nextLogProb > it2->second ){
					
					it2->second = nextLogProb;
					max_arg[t+1].find( nextState )->second = curState;
					delete nextState;
				}
			}
		}
		
		t = t+1;
		if( max_marg[t].size() == 1 ){ // only END state left
			break;
		}
	}
	int L = max_marg.size();
	double maxLogLike =  max_marg[t].begin()->second;
	
	//backward pass (recover every assignment)
	output.resize( L );
	output[L-1] = max_arg[L-1].begin()->first;
	for(int t=L-2; t>=0; t--){
		
		output[t] = max_arg[t+1].find(output[t+1])->second;
	}
	
	for(int t=L-1; t>=1; t--){
		if( output[t-1]->s_type == State::END )
			output.pop_back();
	}
	
	//delete states
	/*for(int i=0;i<L;i++){
		SparseDist::iterator it;
		for(it=max_marg[i].begin(); it != max_marg[i].end(); it++)
			delete it->first;
	}*/
	
	return maxLogLike;
}

double align_Viterbi( vector<char>& seq, vector<char>& ref, State* start_state  ){
	
	vector<State*> output;
	return align_Viterbi( seq, ref, start_state, output );
}


double forward_backward( vector<char>& seq, vector<char>& ref,  State* start_state){
			//vector<SparseDist>& marg_forward ){//, vector<SparseDist>& marg_backward ){
	
	//initialize
	vector<SparseDist> marg_forward;
	vector<pair<State*,double> > init_marg = start_state->transProb(seq,ref);
	marg_forward.push_back(SparseDist());
	//marg_backward.push_back(SparseDist());
	
	////initial transition probability
	for(int i=0;i<init_marg.size();i++){
		marg_forward[0].insert( init_marg[i] );
	}
	
	//forward pass
	int t=0;
	while(1){
		
		////times eimission probability
		SparseDist::iterator it;
		vector<State*> del_vec;
		for(it=marg_forward[t].begin(); it!=marg_forward[t].end(); it++){

			State* s_t = it->first;
			it->second *= s_t->emitProb(seq,ref);
			if( it->second == 0.0 )
				del_vec.push_back(s_t);
		}
		for(int i=0;i<del_vec.size();i++)
			marg_forward[t].erase(del_vec[i]);
		
		
		//compuat max_marg[t+1]
		marg_forward.push_back(SparseDist());
		
		for(it=marg_forward[t].begin();it!=marg_forward[t].end();it++){ //for all current states
			
			State* curState = it->first;
			double marg_t = it->second;
			vector<pair<State*,double> > transVec = curState->transProb(seq,ref);
			
			for(int j=0;j<transVec.size();j++){ // transit to their next states
				 
				State* nextState = transVec[j].first;
				double nextProb = marg_t * transVec[j].second;
				
				SparseDist::iterator it2;
				if( (it2=marg_forward[t+1].find( nextState )) == marg_forward[t+1].end()  ){
					
					marg_forward[t+1].insert( make_pair(nextState, nextProb) );
				}
				else{
					it2->second += nextProb;
					delete nextState;
				}
			}
		}
		
		t = t+1;
		if( marg_forward[t].size() == 1 ) // only END state left
			break;
	}
	int L = marg_forward.size();
	double margLogLike =  log(marg_forward[L-1].begin()->second);
	
	//delete states
	for(int i=0;i<L;i++){
		SparseDist::iterator it;
		for(it=marg_forward[i].begin(); it != marg_forward[i].end(); it++)
			delete it->first;
	}

	return margLogLike;
}

void writeOutput(char* file_out, vector<vector<char> >& seqs, vector<char>& ref, vector<vector<State*> >& align_out){
	
	int N = seqs.size();
	
	vector<char> output_ref;
	vector<vector<char> > output_mat;
	for(int i=0;i<N;i++)
		output_mat.push_back(vector<char>());
	
	// state by state, output alignments
	int* read_i = new int[N];
	for(int n=0;n<N;n++)
		read_i[n] = 0;
	
	int cur_state_i = 0;
	while( cur_state_i < ref.size() ){
		
		//handle INSERT before current state
		bool has_insert;
		do{
			has_insert = false;
			for(int n=0;n<N;n++){
				State* s = align_out[n][read_i[n]];
				if( s->s_type==State::I ){
					has_insert = true;
					output_mat[n].push_back( seqs[n][s->r_i] );
					read_i[n]++;
				}else{
					output_mat[n].push_back( ' ' );
				}
			}
			output_ref.push_back( ' ' );
			
		}while( has_insert );
		
		//delete the last unnecessary blankspace
		output_ref.pop_back();
		for(int n=0;n<N;n++)
			output_mat[n].pop_back();

		//output one column of alignment (of MATCH/DELETION)
		output_ref.push_back(ref[cur_state_i]);
		for(int n=0;n<N;n++){
			
			State* s = align_out[n][read_i[n]];
			if( s->s_type == State::M )
				output_mat[n].push_back( seqs[n][s->r_i] );
			else if( s->s_type == State::D )
				output_mat[n].push_back( '-' );
			else
				output_mat[n].push_back( '?' );
			
			read_i[n]++;
		}
		cur_state_i++;
	}
	//handle INSERT at the end
	bool has_insert;
	do{
		has_insert = false;
		for(int n=0;n<N;n++){
			State* s = align_out[n][read_i[n]];
			if( s->s_type==State::I ){
				has_insert = true;
				output_mat[n].push_back( seqs[n][s->r_i] );
				read_i[n]++;
			}else{
				output_mat[n].push_back( ' ' );
			}
		}
		output_ref.push_back( ' ' );

	}while( has_insert );
	
	//Output to file
	ofstream fout(file_out);
	////output reference sequence
	for(int i=0;i<output_ref.size();i++)
		fout << output_ref[i] ;
	fout << endl;
	for(int i=0;i<output_ref.size();i++){
		fout << "-";
	}fout << endl;
	////output aligned sequences
	for(int n=0;n<N;n++){
		for(int i=0;i<output_mat[n].size();i++)
			fout << output_mat[n][i];
		fout << endl;
	}
	fout.close();
}
