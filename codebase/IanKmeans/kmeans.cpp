#include "kmeans.h"

void read2D(char* fname, vector<Instance*>& data){
	
	ifstream fin(fname);
	int tmp;
	double x,y;
	int id=0;
	while(!fin.eof()){
		
		fin >> tmp >> x >> y;
		if( fin.eof() )break;

		Instance* ins = new Instance(id++);
		ins->fea.push_back(make_pair(1,x));
		ins->fea.push_back(make_pair(2,y));
		
		data.push_back(ins);
	}
	fin.close();
}


void readFixDim(char* fname, vector<Instance*>& data, int D){
	
	ifstream fin(fname);
	int tmp;
	double val;
	int id=0;
	while(!fin.eof()){
		
		fin >> tmp;
		if( fin.eof() )break;
		
		Instance* ins = new Instance(id++);
	       	for(int j=0;j<D;j++){
			fin >> val;
			ins->fea.push_back(make_pair(j+1,val));
		}
		
		data.push_back(ins);
		id++;
	}
	fin.close();
}

void writeResult(char* fname, int* labels, vector<Instance*>& data){

	ofstream fout(fname);
	for(int i=0;i<data.size();i++){
		fout << labels[i] << " ";
	       	for(int j=0;j<data[i]->fea.size();j++)
			fout << data[i]->fea[j].second << " ";
		fout << endl;
	}
	fout.close();
}

void writeModel(char* fname, map<int,Cluster*>& clusters){

	ofstream fout(fname);
	for(map<int,Cluster*>::iterator it = clusters.begin();
		it!=clusters.end(); it++){
		
		fout << it->first << " ";
		Cluster* cluster = it->second;
		
		int size;
		if( (size=cluster->members.size()) != 0 ){

			for(int i=0;i<cluster->d;i++)
				fout << cluster->sum[i] / size  << " ";
		}
		else{
			for(int i=0;i<cluster->d;i++)
				fout << cluster->sum[i] << " ";
		}

		fout << endl;
	}
	
	fout.close();
}

/* z is N by 1 cluster assignment, -1 if not assigned
 * mu is K by D cluster centers
 * E-step first
 */
void kmeans_iter( vector<Instance*>& data, int* z, map<int,Cluster*>& clusters ){
	
	
	int N = data.size();
	int K = clusters.size();
	vector<int> cluster_ids;
	for(map<int,Cluster*>::iterator it=clusters.begin();
		it != clusters.end(); it++){
		cluster_ids.push_back(it->first);
	}
	
	//Main Loop
	int max_iter = 100;
	int iter=0;
	bool z_changed = 1;
	int* z_old = new int[N];
	while( 1 ){
		
		//// E-Step: Find closest cluster z_new_i* for each instance x_i
		z_changed = 0;
		
		for(int i=0;i<data.size();i++){
			
			//compute dist to each cluster, find max
			double dist_min=INF;
			int kmin;
			double dist;
			int csize;
			for(int k=0;k<K;k++){
				
				Cluster* c = clusters.find(cluster_ids[k])->second;
				dist = c->dist( data[i] );
				if( dist < dist_min ){
					dist_min = dist;
					kmin = k;
				}
			}
			z_old[i] = z[i];
			z[i] = cluster_ids[kmin];
			
			if( z[i] != z_old[i] )	
				z_changed = 1;
		}
		
		cout << "iter=" << iter << ", loss=" << obj_loss(data,z,clusters) << endl;
		
		if( iter > max_iter || !z_changed )
			break;
		
		//// M-Step: Adjust each cluster's  sum  and  members based on z (by add/remove)
		for(int i=0;i<N;i++){

			if( z[i] == z_old[i] )
				continue;
			
			if(z_old[i]!=-1){
				Cluster* c_out = clusters.find( z_old[i] )->second;
				c_out->remove( data[i] );
			}

			Cluster* c_in = clusters.find( z[i] )->second;
			c_in->add( data[i] );
		}
		
		iter++;
	}
	
	delete[] z_old;
}

void kmeans( vector<Instance*>& data, int D, int K,   int* z, map<int,Cluster*>& clusters ){
	
	// Random Initialize K clusters ~ uni( [0,1]^D )
	clusters.clear();
	for(int k=0;k<K;k++){
		Cluster* cluster = new Cluster(D); //with 0 being bias
		clusters.insert( make_pair(k, cluster) );
	}
	for(int i=0;i<data.size();i++)
		z[i] = -1;
	kmeans_iter( data, z, clusters);
}

int main(int argc, char** argv){
	
	if( argc < 3 ){
		cerr << "usage: kmeans [dataFile] [K]" << endl;
		cerr << "(assume features scaled to [0,1])" << endl;
		exit(0);
	}
	char* dataFile = argv[1];
	int K = atoi(argv[2]);

	vector<Instance*> data;
	readFixDim( dataFile, data, 13 );
	
	int max_dim = -1;
	for(int i=0;i<data.size();i++){
		
		vector<pair<int,double> >* f = &(data[i]->fea);
		int last_index = f->size()-1;
		if( f->at(last_index).first > max_dim )
			max_dim = f->at(last_index).first;
	}
	int D = max_dim;
	cerr << "D=" << D << endl;
	cerr << "N=" << data.size() << endl;
	cerr << "K=" << K << endl;
	int seed = time(NULL);
	srand(seed);
	cerr << "seed="  << seed << endl;

	/// Run K-means
	map<int,Cluster*> clusters;
	int* z = new int[data.size()];
	kmeans( data, D, K,  z, clusters);
	
	//output results
	writeModel("model.kmeans", clusters);
	writeResult("result.kmeans", z, data );

	return 0;
}
