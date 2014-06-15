/*###############################################################
## MODULE: sparseKmeans.cpp
## VERSION: 1.0 
## SINCE 2014-06-14
## AUTHOR:
##     Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
## DESCRIPTION: 
##      
#################################################################
## Edited by MacVim
## Class Info auto-generated by Snippet 
################################################################*/

#include "sparseClustering.h"


void matadd (double * src1, double * src2, double * dest, int nRows, int nCols) {
    
    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            dest[i][j] = src1[i][j] + src2[i][j];
        }
    }

}

void matsub (double * src1, double * src2, double * dest, int nRows, int nCols) {
    
    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            dest[i][j] = src1[i][j] - src2[i][j];
        }
    }

}

void matdot (double scalar, double * src, double * dest, int nRows, int nCols) {

    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            dest[i][j] = scalar * src[i][j];
        }
    }

}

double L2norm (Instance * ins1, Instance * ins2, int N) {

    double * vec1 = new double [N];
    double * vec2 = new double [N];

    for (int i = 0; i < ins1->fea.size(); i ++) {
        vec1[ ins1->fea[i].first ] = ins1->fea[i].second;
    }
    for (int i = 0; i < ins2->fea.size(); i ++) {
        vec2[ ins2->fea[i].first ] = ins2->fea[i].second;
    }
    
    double norm = 0.0;
    for (int i = 0; i < N; i ++) {
        norm += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }

    return norm;
}

double opt_objective (vector<Instance*>& data, double lambda, int N, double * z) {
    // N is number of entities in "data", and z is N by N.
    // z is current valid solution (average of w_1 and w_2)
    
    // STEP ONE: compute loss function
    // TODO: optimize in terms of symmetry
    double normSum = 0.0;
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            Instance * xi = data[i] ;
            Instance * muj = data[j];
            double dist = L2norm (xi, muj, N);
            normSum += z[i][j] * dist;
        }
    }
    // STEP TWO: compute group-lasso regularization
    double * maxn = new double [N]; 
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            if (z[j][i] > maxn[i])
                maxn[i] = z[j][i];
        }
    }
    double sumk = 0.0;
    for (int i = 0; i < N; i ++) {
        sumk = maxn[i];
    }
    reg = lambda * sumk; 

    return loss + reg;
}

void sparseClustering ( vector<Instance*>& data, int D, int N, double lambda, double * W) {
    double alpha = 0.1;
    // iterative optimization 
    double error = INF;
    while () { 
        // STEP ONE: resolve w_1 and w_2
        frank_wolf ();
        blockwise_closed_form (); 

        // STEP TWO: update z by w_1 and w_2
        double * z = new double [N][N];
        matadd (wone, wtwo, z, N, N);
        // STEP THREE: update the y_1 and y_2 by w_1, w_2 and z
        double * diff = new double [N][N];
        matsub(wone, zone, diffone, N, N);
        matdot(alpha, diffone, diffone, N, N);
        yone = matsub (yone, diffone, yone, N, N);

        double * diff = new double [N][N];
        matsub(wtwo, ztwo, diffone, N, N);
        matdot(alpha, difftwo, difftwo, N, N);
        ytwo = matsub (ytwo, difftwo, ytwo, N, N);

        // STEP FOUR: trace the objective function
        error = opt_objective (data, lambda, N, z);
        cout << "Overall Error: " << error << endl;
    }
}

// entry main function
int main (int argc, char ** argv) {
    
    // exception control
    if (argc < 2) {
        cerr << "Usage: sparseClustering [dataFile] [lambda]" << endl;
        cerr << "Note: dataFile must be scaled to [0,1] in advance." << endl;
        exit(-1);
    }

    // parse arguments
    char * dataFile = argv[1];
    double lambda = atof(argv[2]);

    // read in data
    vector<Instance*> data;
    read2D (dataFile, data);

    // explore the data 
    int dimensions = -1;
    int N = data.size(); // data size
    for (int i = 0; i < N; i++) {
        vector< pair<int,double> > * f = &(data[i]->fea);
        int last_index = f->size() - 1;
        if (f->at(last_index).first > dimensions) {
            dimensions = f->at(last_index).first;
        }
    }

    int D = dimensions;
    cout << "D = " << D << endl;
    cout << "N = " << N << endl;
    cout << "lambda = " << lambda << endl;
    int seed = time(NULL);
    srand(seed);
    cout << "seed = " << seed << endl;

    // Run sparse convex clustering
    map<int, Cluster*> clusters;
    double * W = new double [N][N];

    sparseClustering (data, D, N, lambda, W);
    // Output results

}
