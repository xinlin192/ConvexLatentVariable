#include "TKF1.h"
#include <map>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;


class State{ //interface
	
	public:
	//M:Match, I:Insert, D:Delete
	enum StateType { M,I,D, BEGIN, END };
	StateType s_type;
	int s_i; //state index
	int r_i; //read index
	
	//Emission Probability P(x[index]|this state), i-th variable
	virtual double emitLogProb(vector<char>& x,vector<char>& ref_seq)=0;
	virtual double emitProb(vector<char>& x, vector<char>& ref_seq)=0;
	virtual vector<pair<State*,double> > transLogProb(vector<char>& x, vector<char>& ref_seq)=0;
	virtual vector<pair<State*,double> > transProb(vector<char>& x, vector<char>& ref_seq)=0;
	virtual bool compareTo(State* const &b)=0;
	virtual void print(ostream& out)=0;
};

void print_dist( vector<pair<State*,double> > &vec, ostream &out ){
	
	out << "[";
	for(int i=0;i<vec.size();i++){
		vec[i].first->print(out);
		out << ":";
		out << vec[i].second;
		out << ",";
	}
	out << "]";
}


/**When the root sequence is known, and the mutation, insertion, deletion probabiliy are unknown/known
 * we use PairHMM to align each sequence to the reference sequence independently and exactly.
 */
class PairHMMState:public State{
	
	public:
	
	static double sub_rate;
	static double del_rate;
	static double ins_rate;
	static double alpha;
	static double beta;
	static double gamma;

	static double subProb;
	
	static double match_prob;
	static double ins_prob;
	static double del_prob;

	static double match_prob_del;
	static double ins_prob_del;
	static double del_prob_del;
	
	static void setupPairHMM(double s, double i, double d){
		sub_rate = s;
		ins_rate = i;
		del_rate = d;

		alpha = compute_alpha(s);
		beta = compute_beta(i,d);
		gamma = compute_gamma(i,d);
		
		subProb = 1-exp(-s);
		
		match_prob = (1-beta)* ins_rate/del_rate*alpha;
		match_prob_del = (1-gamma)* ins_rate/del_rate*alpha;
	
		ins_prob = beta;
		ins_prob_del = gamma;

		del_prob = (1-beta)* ins_rate/del_rate*(1-alpha);
		del_prob_del = (1-gamma)* ins_rate/del_rate*(1-alpha);
	}

	PairHMMState(StateType _s_type, int _s_i, int _r_i){
		s_type = _s_type;
		s_i = _s_i;
		r_i = _r_i;
	}
	
	virtual void print(ostream& out){
		
		out << "(";
		if( s_type == M)
			out << "M" ;
		else if( s_type == I )
			out << "I" ;
		else if( s_type == D )
			out << "D";
		else if( s_type == BEGIN )
			out << "B";
		else if( s_type == END )
			out << "E";
		else
			out << s_type;

		out << "," << s_i << "," << r_i << ")";
	}
	
	virtual bool compareTo(State* const &b_abs){

		PairHMMState* b = (PairHMMState*)b_abs;

		// M > I > D > BEGIN
		if( s_type < b->s_type ){
			return 1;
		}else if( s_type > b->s_type){
			return 0;
		}else{
			if( s_i < b->s_i ){
				return 1;
			}else if( s_i > b->s_i ){
				return 0;
			}else{
				if( r_i < b->r_i )
					return 1;
				else
					return 0;
			}
		}
	}
	
	virtual double emitLogProb(vector<char>& x,vector<char>& ref_seq){
		
		return log(emitProb(x,ref_seq));
	}

	virtual double emitProb(vector<char>& x, vector<char>& ref_seq){ //Assume BEGIN won't be called for this
		
		if( s_i >= ref_seq.size() && r_i < x.size() && s_type != I )
			return 0.0;
		if( r_i >= x.size() && s_i < ref_seq.size() && s_type != D)
			return 0.0;

		if( s_type == D )
			return 1.0;
		if( s_type == I )
			return 0.25;
		if( s_type == END )
			return 1.0;
		
		char c_x = x[r_i];
		char c_ref = ref_seq[s_i];
		if( c_x == c_ref ){
			return 1.0-3.0*subProb/4.0;
		}else{
			return subProb/4.0;
		}
	}
	
	virtual vector<pair<State*,double> > transLogProb(vector<char>& x, vector<char>& ref_seq){
		
		vector<pair<State*, double> > vec = transProb(x,ref_seq);
		for(int i=0;i<vec.size();i++){
			vec[i].second = log( vec[i].second );
		}
		return vec;
	}

	virtual vector<pair<State*,double> > transProb(vector<char>& x, vector<char>& ref_seq){
		
		vector<pair<State*,double> > next_states;
		
		if( s_type == END ){
			next_states.push_back(make_pair(new PairHMMState(END,s_i,r_i),1.0));
			return next_states;
		}
		
		/*
		if( s_i >= ref_seq.size() ){
			if( r_i+1 >= x.size() ){ //transit to END state
				
				if(s_type!=D)
					next_states.push_back( make_pair( 
						new PairHMMState(END,s_i,r_i), (1-ins_rate/del_rate)*(1-beta)
					));
				else
					next_states.push_back( make_pair(
						new PairHMMState(END,s_i,r_i), (1-ins_rate/del_rate)*(1-gamma)
					));
				
				return next_states;

			}else{ //only INSERT is possible to transit to
				
				next_states.push_back( make_pair( 
					new PairHMMState(I,s_i,r_i+1), (s_type!=D)?ins_prob:ins_prob_del
				));
				return next_states;
			}
		}*/
		/*
		if( r_i+1>= x.size() ){ // r_i ends but s_i does not
			
			//only Deletion possible
			next_states.push_back( make_pair( 
				new PairHMMState(D,s_i+1,r_i), (s_type!=D)?del_prob:del_prob_del
			));
			return next_states;
		}*/
		
		if( s_type==BEGIN ){
			next_states.push_back(make_pair( new PairHMMState(M,0,0), match_prob ));
			next_states.push_back(make_pair( new PairHMMState(I,0,0), ins_prob ));
			next_states.push_back(make_pair( new PairHMMState(D,0,0), del_prob ));
		}
		else if( s_type==M ){
			if( s_i+1 >= ref_seq.size() && r_i+1 >= x.size() ){
				next_states.push_back(make_pair(
						new PairHMMState(END,s_i+1,r_i+1), (1-ins_rate/del_rate)*(1-beta)
				));
				return next_states;
			}
			next_states.push_back(make_pair( new PairHMMState(M,s_i+1,r_i+1), match_prob ));
			next_states.push_back(make_pair( new PairHMMState(I,s_i+1,r_i+1), ins_prob ));
			next_states.push_back(make_pair( new PairHMMState(D,s_i+1,r_i+1), del_prob ));
		}
		else if( s_type==I ){
			if( s_i >= ref_seq.size() && r_i+1 >= x.size() ){
				next_states.push_back(make_pair(
						new PairHMMState(END,s_i,r_i+1), (1-ins_rate/del_rate)*(1-beta)
				));
				return next_states;
			}
			next_states.push_back(make_pair( new PairHMMState(M,s_i,r_i+1), match_prob ));
			next_states.push_back(make_pair( new PairHMMState(I,s_i,r_i+1), ins_prob ));
			next_states.push_back(make_pair( new PairHMMState(D,s_i,r_i+1), del_prob ));
		}
		else if( s_type==D ){
			if( s_i+1 >= ref_seq.size() && r_i >= x.size() ){
				next_states.push_back(make_pair(
						new PairHMMState(END,s_i+1,r_i), (1-ins_rate/del_rate)*(1-gamma)
				));
				return next_states;
			}
			next_states.push_back(make_pair( new PairHMMState(M,s_i+1,r_i), match_prob_del ));
			next_states.push_back(make_pair( new PairHMMState(D,s_i+1,r_i), ins_prob_del ));
			next_states.push_back(make_pair( new PairHMMState(I,s_i+1,r_i), del_prob_del ));
		}
		
		return next_states;
	}
};

double PairHMMState::sub_rate;
double PairHMMState::del_rate;
double PairHMMState::ins_rate;
double PairHMMState::alpha;
double PairHMMState::beta;
double PairHMMState::gamma;
double PairHMMState::subProb;
double PairHMMState::match_prob;
double PairHMMState::ins_prob;
double PairHMMState::del_prob;
double PairHMMState::match_prob_del;
double PairHMMState::ins_prob_del;
double PairHMMState::del_prob_del;

struct StateCompare{
	public:
	bool operator()(State* const &a, State* const &b){
		return a->compareTo(b);
	}
};

typedef map<State*,double,StateCompare> SparseDist;
typedef map<State*,State*,StateCompare> SparseTable;

template<class E>
void print_dist( map<State*,E,StateCompare>& m, ostream& out ){
	
	out << "{";
	SparseDist::iterator it;
	for(it=m.begin();it!=m.end();it++){
		it->first->print(out);
		out << ":";
		out << it->second;
		out << ",";
	}
	out << "}" ;
}

void print_states( vector<State*>& vec, ostream& out ){
	
	out << "{";
	vector<State*>::iterator it;
	for(it=vec.begin();it!=vec.end();it++){
		(*it)->print(out);
		out << ",";
	}
	out << "}" ;
}
