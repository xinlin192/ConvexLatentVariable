/*###############################################################
## MODULE: hdp_medoids.cpp
## VERSION: 2.0 
## SINCE 2014-07-21
## AUTHOR:
##     Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
##
## DESCRIPTION: 
##     Convex relaxation for hdp_medoids resolution
################################################################# 
## Edited by MacVim
## Class Info auto-generated by Snippet 
################################################################*/

#include "hdp_medoids.h"
using namespace std;

/* Compute the mutual distance of input instances contained within "data" */
void compute_dist_mat (Esmat* dist_mat, Lookups* tables, int N, int D) {
    // STEP ZERO: parse input
    vector< pair<int,int> > doc_lookup = *(tables->doc_lookup);
    vector< pair<int,int> > word_lookup = *(tables->word_lookup); 
    // STEP ONE: compute distribution for each document
    vector< map<int, double> > distributions (D, map<int, double>());
    for (int d = 0; d < D; d ++) {
        // a. compute sum of word frequency
        int sumFreq = 0;
        for (int w = doc_lookup[d].first; w < doc_lookup[d].second; w++) {
            sumFreq += word_lookup[w].second;
        }
        // b. compute distribution
        for (int w = doc_lookup[d].first; w < doc_lookup[d].second; w++) {
            int voc_index = word_lookup[w].first;
            double prob = 1.0 * word_lookup[w].second / sumFreq;
            distributions[d].insert(pair<int, double> (voc_index, prob));
        }
    }
    // STEP TWO: compute weight of word within one document
    esmat_zeros(dist_mat);
    for (int j = 0; j < D; j ++) {
        for (int d = 0; d < D; d ++) {
            for (int w = doc_lookup[d].first; w < doc_lookup[d].second; w++) {
                int voc_index = word_lookup[w].first;
                int count_w_d1 = word_lookup[w].second;    
                double prob_w_d2;
                double dist;
                map<int, double>::const_iterator iter;
                iter = distributions[j].find(voc_index);
                if (iter == distributions[j].end()) {
                    prob_w_d2 = 0.0;
                    dist = INF;
                } else {
                    prob_w_d2 = iter->second;
                //   dist(w, d2) =  - count_w_d1 * log( prob_d2(w) )
                    dist = - count_w_d1 * log(prob_w_d2);
                }
                int esmat_index = w + N * j;
#ifdef DIST_MAT_DUMP
                cout << "index: " << esmat_index << ", dist: " << dist << ", ";
                cout << "count_w_d1: " << count_w_d1 << ", prob_w_d2: " << prob_w_d2 << endl;
#endif
                dist_mat->val.push_back(make_pair(esmat_index, dist));
            }
        }
    }
}
/* dummy_penalty = r dot (1 - sum_k w_nk) */
double get_dummy_loss (Esmat* Z) {
    Esmat* temp_vec = esmat_init (Z->nRows, 1);
    esmat_sum_row (Z, temp_vec);
    double dummy_loss = esmat_compute_dummy (temp_vec);
    esmat_free (temp_vec);
    return dummy_loss;
}
/* \lambda_g \sumk \maxn |\wnk| */
double get_global_topic_reg (Esmat* absZ, double lambda) {
    if (absZ->val.size() == 0) {
        return 0.0;
    }
    Esmat* maxn = esmat_init (1, absZ->nCols);
    Esmat* sumk = esmat_init (1, 1);
    esmat_max_over_col (absZ, maxn);
    esmat_sum_row (maxn, sumk);
    double global_topic_reg = -INF;
    if (sumk->val.size() > 0)
        global_topic_reg = lambda * sumk->val[0].second; 
    else 
        global_topic_reg = 0.0;
    esmat_free (sumk);
    esmat_free (maxn);
    return global_topic_reg;
}
/* \lambdal \sum_d \sum_k \underset{n \in d}{\text{max}} |\wnk| */
double get_local_topic_reg (Esmat* absZ, double lambda, vector< pair<int,int> >* doc_lookup) {
    // STEP ONE: initialize sub matrix for each document
    int nDocs = doc_lookup->size();
    vector<Esmat*> sub_absZ (nDocs);
    for (int d = 0; d < nDocs; d ++) {
        sub_absZ[d] = esmat_init (0,0);
    }
    // STEP TWO: separate entire Z to submat Z
    esmat_submat_row (absZ, sub_absZ, doc_lookup);
    // STEP THREE: compute global topic regularizer for each localized doc
    double local_topic_reg = 0.0;
    for (int d = 0; d < nDocs; d ++) {
        // cout << "size[" << d << "]: " << sub_absZ[d]->val.size() << endl;
        local_topic_reg += get_global_topic_reg (sub_absZ[d], lambda); 
    }
    // Final: free resource
    esmat_free_all (sub_absZ);
    return local_topic_reg;
}
double subproblem_objective (int prob_index, Esmat* Y, Esmat* Z, Esmat* W, double RHO, double lambda, Lookups* tables, Esmat* dist_mat) {
    string title = "";
    vector< pair<int,int> >* doc_lookup = tables->doc_lookup;
    vector< pair<int,int> >* word_lookup = tables->word_lookup; 
    vector< vector<int> >* voc_lookup = tables->voc_lookup;
    /*
       esmat_print (W, "W["+ to_string(prob_index) + "]")
       esmat_print (Y, "Y["+ to_string(prob_index) + "]")
       esmat_print (Z, "Z["+ to_string(prob_index) + "]")
       */

    // STEP ONE: compute main term
    double main = -1.0;
    if (prob_index == 1) {
        double loss = esmat_frob_prod (dist_mat, Z);
#ifdef SUBPROBLEM_DUMP
        cout << "loss: " << loss << ", ";
#endif
        title = "dummy";
        // cout << "begin to compute dummy loss" << endl;
        // dummy_penalty = r dot (1 - sum_k w_nk)
        main = get_dummy_loss (W) + loss;
    } else {
        Esmat* absW = esmat_init (W);
        esmat_abs (W, absW);
        if (prob_index == 2) {
            title = "Global_Reg";
            main = get_global_topic_reg (absW, lambda);
        } else if (prob_index == 3) {
            title = "Local_Reg";
            main = get_local_topic_reg (absW, lambda, doc_lookup);
        }
        esmat_free (absW);
    }
    Esmat* w_minus_z = esmat_init ();
    // STEP TWO: compute linear term: linear = y^T dot (w - z) 
    esmat_sub (W, Z, w_minus_z); // temp = w - z
    double linear = esmat_frob_prod (Y, w_minus_z);
    // STEP THREE: compute quadratic term: quadratic = 0.5 * RHO * || w - z ||^2 
    double quadratic = 0.5 * RHO * esmat_frob_norm (w_minus_z);
    total = main + linear + quadratic;
    esmat_free (w_minus_z);
#ifdef SUBPROBLEM_DUMP
    cout << title << ": " << main << ", ";
    cout << "linear: " << linear << ", ";
    cout << "quadratic: " << quadratic << ", ";
    cout << "total: " << total << endl;
#endif
    return total;
}
double original_objective (Esmat* Z, vector<double> LAMBDAs, Lookups* tables) {
    vector< pair<int,int> >* doc_lookup = tables->doc_lookup;
    vector< pair<int,int> >* word_lookup = tables->word_lookup; 
    vector< vector<int> >* voc_lookup = tables->voc_lookup;

    Esmat* absZ = esmat_init (Z);
    esmat_abs (Z, absZ);
    // STEP ONE: compute dummy loss
    double dummy = get_dummy_loss (Z);
    // cout << "dummy =" << dummy << endl;

    // STEP TWO: compute "GLOBAL TOPIC" group-lasso regularization
    double global_topic_reg = get_global_topic_reg (absZ, LAMBDAs[0]);
    esmat_free (absZ); 
    return dummy + global_topic_reg;
}

/* three subproblems that employed group_lasso_solver in different ways */
void global_topic_subproblem (Esmat* Y, Esmat* Z, Esmat* w, double RHO, double lambda) {
    group_lasso_solver (Y, Z, w, RHO, lambda);
}
void local_topic_subproblem (Esmat* Y, Esmat* Z, Esmat* w, double RHO, double lambda, Lookups* tables) {
    vector< pair<int,int> >* doc_lookup = tables->doc_lookup;
    vector< pair<int,int> >* word_lookup = tables->word_lookup; 
    vector< vector<int> >* voc_lookup = tables->voc_lookup;

#ifdef LOCAL_SUBPROBLEM_DUMP
    cout << "[local w input]" << w->nRows << "," << w->nCols << "," << w->val.size() << endl;
    esmat_toString (w);
#endif

    int nDocs = tables->nDocs;
    Esmat* tempW = esmat_init (w);
    vector<Esmat*> subW (nDocs); 
    vector<Esmat*> subY (nDocs);
    vector<Esmat*> subZ (nDocs);
    // STEP ZERO: initialize all submats
    for (int d = 0; d < nDocs; d ++) {
        subW[d] = esmat_init (0,0);
        subY[d] = esmat_init (0,0);
        subZ[d] = esmat_init (0,0);
    }
    // STEP ONE: separate esmat Y, Z, w to multiple small-sized esmat
    // NOTE: all esmat position has to be recomputed
    esmat_submat_row (w, subW, doc_lookup);
    esmat_submat_row (Y, subY, doc_lookup);
    esmat_submat_row (Z, subZ, doc_lookup);

    for (int d = 0; d < nDocs; d ++) {
#ifdef LOCAL_SUBPROBLEM_DUMP
        cout << "subW[d" << d << "]" << esmat_toInfo(subW[d]);
        cout << "subY[d" << d << "]" << esmat_toInfo(subY[d]);
        cout << esmat_toString(subY[d]);
        cout << "subZ[d" << d << "]" << esmat_toInfo(subZ[d]);
        cout << esmat_toString(subZ[d]);
#endif
        // STEP TWO: invoke group_lasso_solver to each individual group
        group_lasso_solver (subY[d], subZ[d], subW[d], RHO, lambda);

#ifdef LOCAL_SUBPROBLEM_DUMP
        cout << "res_subW[d" << d << "]" << esmat_toInfo(subW[d]);
#endif

        // STEP TREE: merge solution of each individual group (place back)
        // NOTE: all esmat position has to be recomputed
        int start_row = (*doc_lookup)[d].first;
        int end_row = (*doc_lookup)[d].second;
        esmat_merge_row (subW[d], start_row, end_row, tempW);
    }
    // realign the mat->val with index-increasing order
    esmat_align (tempW);
    // STEP FIVE: free auxiliary resource
    for (int d = 0; d < nDocs; d ++) {
        esmat_free (subW[d]);
        esmat_free (subY[d]);
        esmat_free (subZ[d]);
    }
#ifdef LOCAL_SUBPROBLEM_DUMP
    cout << "[local w before output]" << w->nRows << "," << w->nCols << "," << w->val.size() << endl;
    cout << "[tempW] " << esmat_toInfo(tempW);
#endif

    // FINAL: update merged solution to w
    esmat_copy (tempW, w);

#ifdef LOCAL_SUBPROBLEM_DUMP
    cout << "[local w output]" << w->nRows << "," << w->nCols << "," << w->val.size() << endl;
#endif
}
void hdp_medoids (Esmat* dist_mat, vector<double> LAMBDAs, Esmat* W, Lookups* tables) {
    // SET MODEL-RELEVANT PARAMETERS 
    assert (LAMBDAs.size() == 2);
    double ALPHA = 1.0;
    double RHO = 1.0;
    int N = tables->nWords;
    int D = tables->nDocs;
    /* A. DECLARE AND INITIALIZE INVOLVED VARIABLES AND MATRICES */
    // 1) randomize initial global cluster medoid 
    
    // 2) randomize initial local cluster medoids for each document j

    // 3) initialize local cluster indicator z_ij = 1 for all j = 1

    // 4) initialize global cluster association v_j1 = 1

    Esmat* w = esmat_init (N, D);

    /* B. Repetition until Convergence */
    while () {
        // 1) preprocess distance

        // 2) if ..,  global and local augmentation
        //  if min_p d_ijp > \lambda_l + \lambda_g

        // 3) otherwise, local assignment or local augmentation
        

        // 4) re-elect medoids for all local clusters

        // 5) if .., global augmentation
        //  if min_p d_jcp > \lambda_g + ... 

        // 6) otherwise, update global association


        // 7) re-elect medoids for global cluster

    }

    /* De-allocation */
    esmat_free (w);

    /*  put converged solution to destinated W*/
    esmat_copy (z, W);
    esmat_free (z);
}

// entry main function
int main (int argc, char ** argv) {

    // EXCEPTION control: illustrate the usage if get input of wrong format
    if (argc < 5) {
        cerr << "Usage: " << endl;
        cerr << "\thdp_medoids [voc_dataFile] [doc_dataFile] [lambda_global] [lambda_local]" << endl;
        exit(-1);
    }

    // PARSE arguments
    string voc_file (argv[1]);
    string doc_file (argv[2]);
    vector<double> LAMBDAs (2, 0.0);
    LAMBDAs[0] = atof(argv[3]); // lambda_document
    LAMBDAs[1] = atof(argv[4]); // lambda_topic

    // preprocess the input dataset
    vector<string> voc_list;
    voc_list_read (voc_file, &voc_list);
    int nVocs = voc_list.size();

    // init lookup_tables
    vector< pair<int,int> > doc_lookup;
    vector< pair<int,int> > word_lookup;
    vector< vector<int> > voc_lookup (nVocs, vector<int>());
    Lookups lookup_tables;
    lookup_tables.doc_lookup = &doc_lookup;
    lookup_tables.word_lookup = &word_lookup;
    lookup_tables.voc_lookup = &voc_lookup;

    document_list_read (doc_file, &lookup_tables);
    lookup_tables.nDocs = lookup_tables.doc_lookup->size();
    lookup_tables.nWords = lookup_tables.word_lookup->size();
    lookup_tables.nVocs = nVocs;

    int seed = time(NULL);
    srand (seed);
    cerr << "###########################################" << endl;
    cerr << "nVocs = " << lookup_tables.nVocs << endl; // # vocabularies
    cerr << "nDocs = " << lookup_tables.nDocs << endl; // # documents
    cerr << "nWords = " << lookup_tables.nWords << endl; // # words
    cerr << "lambda_global = " << LAMBDAs[0] << endl;
    cerr << "lambda_coverage = " << LAMBDAs[1] << endl;
    cerr << "TRIM_THRESHOLD = " << TRIM_THRESHOLD << endl;
    cerr << "seed = " << seed << endl;
    cerr << "###########################################" << endl;

    // Run sparse convex clustering
    int N = lookup_tables.nWords;
    int D = lookup_tables.nDocs;
    Esmat* W = esmat_init (lookup_tables.nWords, lookup_tables.nDocs);
    Esmat* dist_mat = esmat_init (N, D);
    compute_dist_mat (dist_mat, &lookup_tables, N, D);
    
    ofstream dmat_out ("dist_mat");
    dmat_out << esmat_toInfo(dist_mat);
    dmat_out << esmat_toString(dist_mat);
    hdp_medoids (dist_mat, LAMBDAs, W, &lookup_tables);

    /* Output objective */
    output_objective(clustering_objective (dist_mat, W));

    /* Output cluster centroids */
    output_model (W);

    /* Output assignment */
    output_assignment (W, &word_lookup);

    /* reallocation */
    esmat_free (W);
    esmat_free (dist_mat);
}
