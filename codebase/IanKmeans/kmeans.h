#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<vector>
#include<map>
#include<set>
#include "util.h"

using namespace std;

class Instance{
	
	public:
	int id;
	int label;
	vector<pair<int,double> > fea;
	Instance(int _id){
		id = _id;
		_x_sq = -1.0;
	}
	
	double x_sq(){
		
		// Pre-calculated
		if( _x_sq != -1.0 ){
			return _x_sq;
		}
		
		// Calculate and Cache
		_x_sq = 0.0;
		double tmp;
		for(int i=0;i<fea.size();i++){
			tmp = fea[i].second;
			_x_sq += tmp*tmp;
		}

		return _x_sq;
	}

	private:
	double _x_sq;
};
/*
class Cluster{

	public:
	int d;
	vector<double>  mu;
	double mu_sq;
	
	Cluster(int _d){
		
		d = _d;
		for(int i=0;i<d+1;i++)
			mu.push_back( ((double)rand()/RAND_MAX)*2-1 );
		mu[0] = 0.0;
		mu_sq = -1;
	}

	Cluster(vector<double>& _mu){
		
		d = _mu.size()-1;
		for(int i=0;i<_mu.size();i++)
			mu.push_back( _mu[i] );
		mu_sq = -1;
	}
	
	double dist(Instance* ins){
		
		if( mu_sq == -1 ){
			mu_sq = 0.0;
			for(int i=0;i<d+1;i++){
				mu_sq += mu[i]*mu[i];
			}
		}
		
		return  ins->x_sq() + mu_sq - 2*dot(mu,ins->fea);
	}


};*/


class Cluster{

	public:
	int d;
	set<Instance*>  members;
	double*  sum;  //sum = mu * |cluster|, if |cluster|!=0 
		       //      mu,             if |cluster|=0.
	double mu_sq;
	
	Cluster(int _d){
		
		d = _d;
		sum = new double[d];
		for(int i=0;i<d;i++)
			sum[i] = (double)rand()/RAND_MAX;
		sum[0] = 0.0;
		
		mu_sq = -1;
	}
	
	~Cluster(){
		delete[] sum;
	}
	
	double dist(Instance* ins){
		
		if( mu_sq == -1 ){

			mu_sq = dot(sum,sum,d);
			if( members.size() != 0 )
				mu_sq /= (members.size()*members.size());
		}
		
		if( members.size() != 0 )
			return  ins->x_sq() + mu_sq - 2*dot(sum,ins->fea)/members.size();
		else
			return  ins->x_sq() + mu_sq - 2*dot(sum,ins->fea);
	}

	void add(Instance* ins){
		
		if( !members.insert(ins).second )
			cerr << "instance " << ins->id << " already in cluster when added." << endl;

		vadd(sum,ins->fea);
		mu_sq = -1;
	}

	void remove(Instance* ins){
		
		if(members.erase(ins)==0)
			cerr << "instance " << ins->id << "not in cluster when removed." << endl;

		vsub(sum,ins->fea);
		mu_sq = -1;
	}
};

class ScoreComparator{
	
	private:
	vector<double>* score_vec;

	public:
	ScoreComparator(vector<double>* vec){
		score_vec = vec;
	}

	bool operator() (int i,int j) { 
		return score_vec->at(i) < score_vec->at(j); 
	}
};


/** loss = \frac{1}{2} sum_{i,j}{ w_{ij}\| x_i - mu_j \|_2^2 },  w_{iz[i]}=1
 */
double obj_loss( vector<Instance*>& data, int* z, map<int,Cluster*>& clusters ){
	
	double loss = 0.0;
	for(int i=0;i<data.size();i++){
		
		Cluster* cluster = clusters.find( z[i] )->second;
		Instance* ins = data[i];

		loss += cluster->dist( ins )/2;
	}
	
	return loss;
}
