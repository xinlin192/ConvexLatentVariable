#include<stdlib.h>
#include<vector>
using namespace std;
const double INF = 1e300;

template <class T>
void shuffle( vector<T>& vec ) {

	int l = vec.size();
	int j;
	for(int i=0; i<l; i++)
	{
		j = i+rand()%(l-i);
		swap(vec[i], vec[j]);
	}
}


double dot(vector<double>& v, vector<pair<int,double> >& f) {
	
	double sum = 0.0;
	for(int i=0;i<f.size();i++)
		sum += v[ f[i].first ] * f[i].second;

	return sum;
}

double dot(double* v, vector<pair<int,double> >& f) {
	
	double sum = 0.0;
	for(int i=0;i<f.size();i++)
		sum += v[ f[i].first ] * f[i].second;

	return sum;
}

double dot(double* v, double* w, int d) {
	
	double sum = 0.0;
	for(int i=0;i<d;i++) {
		sum += v[i]*w[i];
	}
	return sum;
}

void vadd(double* v, vector<pair<int,double> >& f) {
	
	for(int i=0;i<f.size();i++)
		v[ f[i].first ] += f[i].second;
}

void vadd(vector<double>& vec, vector<pair<int,double> >& f) {
	
	for(int i=0;i<f.size();i++)
		vec[ f[i].first ] += f[i].second;
}

void vsub(double* v, vector<pair<int,double> >& f) {
	
	for(int i=0;i<f.size();i++)
		v[ f[i].first ] -= f[i].second;
}

int max_index (vector<double>& vec) {

    double vmax = -INF; 
    int imax;
    for (int i = 0; i < vec.size(); i ++) {
        if (vec[i] > vmax) {
            vmax = vec[i];
            imax = i;
        }
    }

    return imax;
}

int max_index (double* arr, int size) {
	
	double vmax= -INF;
	int imax;
	for(int i=0;i<size;i++) {
		if( arr[i] > vmax ) {
			vmax = arr[i];
			imax = i;
		}
	}
	return imax;
}

double flip (double a) {
	
	return (a==0.0)?1.0:0.0;
}

int flip (int a) {
	return (a==0)?1:0;
}
