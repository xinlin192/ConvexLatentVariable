/*###############################################################
## MODULE: exSparseMat.h
## VERSION: 2.0 
## SINCE: 2014-07-24
## AUTHOR: Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
## DESCRIPTION: 
##     Header file for Extensible Sparse matrix  
#################################################################
## Edited by MacVim
## Class Info auto-generated by Snippet 
################################################################*/

#include<string>
#include<vector>
#include<fstream>
#include<cassert>
#include <algorithm>
#include"math.h"
using namespace std;

/* Global variables */
const double DUMMY_PENALTY_RATE = 1000.0;
const double TRIM_THRESHOLD = 10e-5;

/* Definition of Data structure for Extensible Sparse Matrix (Esmat*) */
typedef struct {
    int nRows; int nCols;
    vector< pair<int, double> >  val;
} Esmat ;

typedef double (* Operation) (double, double); 

/* Prototype for fundamental functions, typically computational frameworks*/
double esmat_unary_operate (Esmat *A, Operation opt);
void esmat_bin_operate (Esmat* A, Esmat* B, Esmat* dest, Operation opt);
void esmat_operate_col (Esmat* A, Esmat* dest, Operation opt, double init_value);
void esmat_operate_row (Esmat* A, Esmat* dest, Operation opt, double init_value);

/* Allocation and De-allocation */
Esmat* esmat_init (int nRows, int nCols);
Esmat* esmat_init (Esmat * A);
Esmat* esmat_init ();
void esmat_init_all (vector<Esmat*>* src);
Esmat* esmat_read (string fname);
bool esmat_equal (Esmat* esmat, double **mat);
void esmat_free (Esmat* src);
void esmat_free_all (vector<Esmat*> src);
void esmat_zeros (Esmat* A);

/* Rearrange one esmat */
void esmat_align (Esmat* mat);
void esmat_copy (Esmat* A, Esmat* D);
void esmat_trim (Esmat* A);

/* submat and merge */
void esmat_submat_row (Esmat* mat, vector<Esmat*> submats, vector< pair<int,int> >* look_up);
void esmat_submat_row (Esmat* mat, vector<Esmat*> submats, vector<int>* word_lookup, vector< vector<int> >* voc_lookup);
void esmat_merge_row (Esmat* submat, int start_index, int end_index, Esmat* mat);
void esmat_merge_row (Esmat* submat, vector<int>* sub_voc_lookup, Esmat* mat);

/* frobenius product and norm */
void esmat_abs (Esmat* A, Esmat* dest);
double esmat_fdot (Esmat* A, Esmat* B);
double esmat_sum (Esmat* A);
double esmat_fnorm (Esmat* A);

/* scalar multiplication */
void esmat_scalar_mult (double scalar, Esmat* A);
void esmat_scalar_mult (double scalar, Esmat* A, Esmat* dest);

/* Auxiliary functions */
bool esmat_isValid (Esmat* A, Esmat* B, int mode);
string esmat_toString (Esmat* A);

/* Compute dummy term */
double esmat_compute_dummy (Esmat* A);

/* Add and Subtract two extensible sparse matrices */
void esmat_add (Esmat* A, Esmat* B, Esmat* dest);
void esmat_sub (Esmat* A, Esmat* B, Esmat* dest);

/* min, max and sum over column and row elements */
void esmat_min_col (Esmat* A, Esmat* dest);
void esmat_max_col (Esmat* A, Esmat* dest);
void esmat_sum_col (Esmat* A, Esmat* dest);
void esmat_min_row (Esmat* A, Esmat* dest);
void esmat_max_row (Esmat* A, Esmat* dest);
void esmat_sum_row (Esmat* A, Esmat* dest); 

/* Powerful function pointers and definition of instantiated operators */
double min (double a, double b) { return a < b ? a : b; }
double max (double a, double b) { return a > b ? a : b; }
double sum (double a, double b) { return a + b; }
double diff (double a, double b) { return a - b; } // non-symmetric 
double times (double a, double b) { return a * b; }
double power2 (double val, double counter) 
    { return counter + val * val; } // non-symmetric
double dummy_penalty (double val, double counter) 
    { return counter + DUMMY_PENALTY_RATE * (1 - val); }
double count (double value, double counter) { 
    // note that for this "count" function, two input is not symmetric
    if (fabs(value) > TRIM_THRESHOLD) counter += 1.0;
    return counter; 
}
bool pair_First_Elem_inc_Comparator (const std::pair<int, double>& firstElem, const std::pair<int, double>& secondElem) {
    // sort pairs by second element with *incremental order*
    return firstElem.first < secondElem.first;
}

/* Allocate one extensible sparse matrix with all zero entries */
Esmat* esmat_init (int nRows, int nCols) {
    Esmat* freshman = new Esmat ();
    freshman->nRows = nRows;
    freshman->nCols = nCols;
    // freshman->val = new vector< pair<int, double> > ();
    return freshman;
}
Esmat* esmat_init (Esmat* src) {
    return esmat_init (src->nRows, src->nCols);
}
Esmat* esmat_init () {
    return esmat_init (0, 0);
}
void esmat_init_all (vector<Esmat*>* src) {
    int N = src->size();
    for (int i = 0; i < N; i++) {
        (*src)[i] = esmat_init();
    }
}
Esmat* esmat_read (string fname) {
   	ifstream fin(fname);
    int nRows, nCols;
    fin >> nRows >> nCols;
    Esmat* result = esmat_init (nRows, nCols);
    double val;
    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            fin >> val;
            if (fabs(val) > TRIM_THRESHOLD) {
                int esmat_index = i + j * nRows;
                result->val.push_back(make_pair(esmat_index, val));
            }
        }
    }
    esmat_align (result);
    fin.close();
    return result;
}
bool esmat_equal (Esmat* esmat, double **mat) {
    bool isequal = true;
    int esmat_size = esmat->val.size();
    int i = 0, esmat_index = esmat->val[i].first, mat_index = 0;
    while (true) {
        esmat_index = esmat->val[i].first;
        if (esmat_index > mat_index) {
            // there is non-zero value in mat but no in esmat, hence not equal
            if (fabs(mat[mat_index % esmat->nRows][mat_index/esmat->nRows]) > TRIM_THRESHOLD) {
                isequal = false; 
                break;
            }
            ++ mat_index;
        } else if (esmat_index == mat_index) {
            // element in the same position is not equal to each other
            if (esmat->val[i].second != mat[mat_index % esmat->nRows][mat_index/esmat->nRows]) {
                isequal = false; 
                break;
            }
            ++ i; ++ mat_index;
        }
    }
    return isequal;
}
/* Deallocate given Esmat */
void esmat_free (Esmat* src) {
    src->val.clear();
    delete src;
}
void esmat_free_all (vector<Esmat*> src) {
    int N = src.size();
    for (int i = 0; i < N; i ++) {
        esmat_free(src[i]);
    }
}
void esmat_zeros (Esmat* A) {
    A->val.clear();
}
/* esmat alignment: place entry of sparse matrix in index-increasing order */
void esmat_align (Esmat* mat) {
    std::sort(mat->val.begin(), mat->val.end(), pair_First_Elem_inc_Comparator);
}

/* submat and merge */
/* pick subset of rows and form a new esmat */
void esmat_submat_row (Esmat* mat, vector<Esmat*> submats, vector< pair<int,int> >* look_up) {
    int nDocs = look_up->size();
    // renew the characteristics of submats
    for (int d = 0; d < nDocs; d ++) {
        int end_index = (*look_up)[d].second;
        int start_index = (*look_up)[d].first;
        submats[d]->nRows = end_index - start_index;
        submats[d]->nCols = mat->nCols;
        submats[d]->val.clear();
    }

    int mat_size = mat->val.size();
    int i = 0, d = 0;
    int start_index = (*look_up)[d].first;
    int end_index = (*look_up)[d].second;
    while (i < mat_size) {
        int mat_index = mat->val[i].first;
        int mat_row_index = mat_index % mat->nRows;
        int mat_col_index = mat_index / mat->nRows;
        // find a corresponding submats
        while (!(start_index <= mat_row_index && mat_row_index < end_index)) {
            ++ d;
            if (d >= nDocs) d = 0;
            start_index = (*look_up)[d].first;
            end_index = (*look_up)[d].second;
        }
        // get value of that entry
        double value = mat->val[i].second;
        // compute corresponding position in submats
        int submat_row_index = mat_row_index - start_index;
        int submat_col_index = mat_col_index;
        int submat_index = submat_row_index + submat_col_index * submats[d]->nRows;
        submats[d]->val.push_back(make_pair(submat_index, value));
        ++ i;
    }
}
// This submat_row version is for coverage_subproblem
void esmat_submat_row (Esmat* mat, vector<Esmat*> submats, vector<int>* word_lookup, vector< vector<int> >* voc_lookup) {
    int nVocs = voc_lookup->size();
    // STEP ONE: renew the characteristics of submats
    for (int v = 0; v < nVocs; v ++) {
        submats[v]->nRows = (*voc_lookup)[v].size();
        submats[v]->nCols = mat->nCols;
        submats[v]->val.clear();
    }
    // STEP TWO:
    int mat_size = mat->val.size();
    for (int i = 0; i < mat_size; i ++) {
        int mat_index = mat->val[i].first;
        int mat_row_index = mat_index % mat->nRows;
        int mat_col_index = mat_index / mat->nRows;

        int v = (*word_lookup)[mat_row_index] - 1;
        vector<int>::iterator it;
        it = find((*voc_lookup)[v].begin(), (*voc_lookup)[v].end(), mat_row_index);
        int submat_row_index = *it;
        int submat_col_index = mat_col_index;
        int submat_index = submat_row_index + submat_col_index * submats[v]->nRows;

        double value = mat->val[i].second;
        submats[v]->val.push_back(make_pair(submat_index, value));
    }
}
/* put submat to specified position of mat */
// this merge_row version is for local_topic_subproblem
void esmat_merge_row (Esmat* submat, int start_index, int end_index, Esmat* mat) {
    assert (start_index < end_index);

    int submat_size = submat->val.size();
    for (int i = 0; i < submat_size; i ++) {
        int submat_index = submat->val[i].first;
        int submat_row_index = submat_index % submat->nRows;
        int submat_col_index = submat_index / submat->nRows;
        // get value of that entry
        double value = submat->val[i].second;
        // compute corresponding position in submat
        int mat_row_index = submat_row_index + start_index;
        int mat_col_index = submat_col_index;
        int mat_index = mat_row_index + mat_col_index * mat->nRows;
        mat->val.push_back(make_pair(mat_index, value));
    }
}
// this merge_row version is for coverage_subproblem
void esmat_merge_row (Esmat* submat, vector<int>* sub_voc_lookup, Esmat* mat) {
    int nWords = sub_voc_lookup->size(); // related to one vocabulary
    assert (submat->nRows == nWords);

    int submat_size = submat->val.size();
    for (int i = 0; i < submat_size; i ++) {
        int submat_index = submat->val[i].first;
        int submat_row_index = submat_index % submat->nRows;
        int submat_col_index = submat_index / submat->nRows;
        
        int mat_row_index = (*sub_voc_lookup)[submat_row_index];
        int mat_col_index = submat_col_index;
        int mat_index = mat_row_index + mat_col_index * mat->nRows;

        double value = submat->val[i].second;
        mat->val.push_back(make_pair(mat_index, value));
    }
}
/* frobenius product */
double esmat_fdot (Esmat* A, Esmat* B) {

    assert (esmat_isValid (A, B, 1));

    int i = 0, j = 0;
    int indexA, indexB;
    int sizeA = A->val.size();
    int sizeB = B->val.size();
    double result = 0.0;

    while (i < sizeA && j < sizeB) {
        indexA = A->val[i].first;
        indexB = B->val[i].first;

        if (indexA < indexB) {
            ++ i; 
        } else if (indexA > indexB) {
            ++ j;
        } else { // equality
            result += A->val[i].second * B->val[j].second;
            ++i; ++j;
        }
    }

    return result;
}
/* abs */
void esmat_abs (Esmat* A, Esmat* dest) {
    dest->nRows = A->nRows;
    dest->nCols = A->nCols;
    dest->val.clear();

    for (int i = 0; i < A->val.size(); i ++) {
        // insert multiplied value to the same position
        dest->val.push_back(make_pair(A->val[i].first, fabs(A->val[i].second))); 
    }
}
/* scalar times a esmat and store on input matrix*/
void esmat_scalar_mult (double scalar, Esmat* A) {
    int sizeA = A->val.size();
    int capacity = A->nRows * A->nCols;
    for (int i = 0; i < sizeA; i ++) {
        assert (A->val[i].first < capacity);
        A->val[i].second *= scalar; 
    }
}

/* scalar times a esmat */
void esmat_scalar_mult (double scalar, Esmat* A, Esmat* dest) {
    dest->nRows = A->nRows; 
    dest->nCols = A->nCols;
    dest->val.clear();
    for (int i = 0; i < A->val.size(); i ++) {
        // insert multiplied value to the same position
        dest->val.push_back(make_pair(A->val[i].first, scalar*A->val[i].second)); 
    }
}

/* Check validity (dim alignment) of input esmat 
 *  mode:
 *    1 - same dim alignment
 *    2 - product alignment
 * */
bool esmat_isValid (Esmat* A, Esmat* B, int mode) {

    bool success = false;

    if (mode == 1) {
        if (A->nRows == B->nRows && A->nCols == B->nCols) 
            success = true;
    } else if (mode == 2) {
        if (A->nCols == B->nRows) 
            success = true;
    }

    return success;
}

string esmat_toString (Esmat* A) {

    assert (A->nRows > 0);
    assert (A->nCols > 0);

    string idx_val_separator = ":";
    string field_seperator = ",";
    string line_separator = "\n";

    vector<string> allStrings (A->nRows, "") ;
    string str = "";

    int sizeA = A->val.size();
    for (int i = 0; i < sizeA; i ++) {
        int overall_idx = A->val[i].first;
        // column major data structure
        int col_idx = overall_idx / A->nRows;
        int row_idx = overall_idx % A->nRows;
        assert (col_idx < A->nCols);
        assert (row_idx < A->nRows);
        // generate newly added string
        string temp = to_string(col_idx) + idx_val_separator + to_string(A->val[i].second);
        // row major string representation
        if (allStrings[row_idx].size() == 0) {
            allStrings[row_idx] = "" + temp;
        } else {
            allStrings[row_idx] += field_seperator + temp;
        }
    }

    for (int i = 0; i < A->nRows; i ++) {
        if (allStrings[i].size() > 0)
            str += allStrings[i] + line_separator;
    }

    return str;
}

/* copy content of Esmat* A to Esmat* D */
void esmat_copy (Esmat* A, Esmat* D) {

    D->nRows = A->nRows;
    D->nCols = A->nCols;
    D->val.clear();

    int sizeA = A->val.size();
    for (int i = 0; i < sizeA; i ++) {
        D->val.push_back(make_pair (A->val[i].first, A->val[i].second));
    }
}

void esmat_trim (Esmat* A) {
    int sizeA = A->val.size();
	for (int i = 0; i < sizeA; i ++) {
        double value = A->val[i].second;
        if ( fabs(value) < TRIM_THRESHOLD ) {
            // remove this index:value pair
            A->val.erase(A->val.begin()+i);
            -- i; -- sizeA;
        }
	}
}

//=========================================================
// The followings are computational frameworks.
// Do not modify them unless you know exactly what you are
// doing. 
// ========================================================
/* Framework of unary operation for Esmat A */
double esmat_unary_operate (Esmat * A, Operation opt) {
    double result = 0.0;
    int sizeA = A->val.size();
    for (int i = 0; i < sizeA; i++) {
        result = opt(A->val[i].second, result);
    }
    return result;
}
/* Framework of binary operation for Esmat A */
void esmat_bin_operate (Esmat* A, Esmat* B, Esmat* dest, Operation opt) {
    assert (esmat_isValid (A, B, 1));

    dest->nRows = A->nRows;
    dest->nCols = A->nCols;
    dest->val.clear();

    int i = 0; int j = 0;
    int indexA, indexB;
    int sizeA = A->val.size();
    int sizeB = B->val.size();

    while (i < sizeA && j < sizeB) {
        indexA = A->val[i].first;
        indexB = B->val[i].first;

        if (indexA < indexB) {
            dest->val.push_back(make_pair(indexA, A->val[i].second));
            ++i; 
        } else if (indexA > indexB) {
            dest->val.push_back(make_pair(indexB, B->val[i].second));
            ++j;
        } else { // equality
            double value = opt(A->val[i].second, B->val[i].second);
            dest->val.push_back(make_pair(indexA, value));
            ++i; ++j;
        }
    }
}
/* Framework for operating over each row of esmat A */
void esmat_operate_row (Esmat* A, Esmat* dest, Operation opt, double init_value) {

    int sizeA = A->val.size();

    // set up dest esmat, here assume it has been initialized
    dest->nRows = A->nRows;
    dest->nCols = 1;
    dest->val.clear();

    vector<double> temp (A->nRows, init_value);
    
    for (int i = 0; i < sizeA; i ++) {
        int row_idx = dest->val[i].first % A->nRows;
        temp[row_idx] = opt(dest->val[i].second, temp[row_idx]);
    }

    for (int j = 0; j < A->nRows; j ++) {
        // trim very small term
        if (fabs(temp[j]) > TRIM_THRESHOLD) 
            dest->val.push_back(make_pair(j, temp[j]));
    }
}
/* Framework for operating over each column of esmat A */
void esmat_operate_col (Esmat* A, Esmat* dest, Operation opt, double init_value) {
    
    int sizeA = A->val.size();

    // set up dest esmat, here assume it has been initialized
    dest->nRows = 1;
    dest->nCols = A->nCols;
    dest->val.clear();

    vector<double> temp (A->nCols, init_value);
    
    for (int i = 0; i < sizeA; i ++) {
        int col_idx = dest->val[i].first / A->nRows;
        temp[col_idx] = opt(dest->val[i].second, temp[col_idx]);
    }

    for (int j = 0; j < A->nCols; j ++) {
        // trim very small term
        if (fabs(temp[j]) > TRIM_THRESHOLD)
            dest->val.push_back(make_pair(j, temp[j]));
    }
}
/* Sum of all element on one matrix */
double esmat_sum (Esmat* A) 
{ return esmat_unary_operate (A, sum); }
/* Frobenius norm of one extensible sparse matrices */
double esmat_fnorm (Esmat* A) 
{ return esmat_unary_operate (A, power2); }
/* Compute dummy term */
double esmat_compute_dummy (Esmat* A) 
{ return esmat_unary_operate (A, dummy_penalty); }

/* Add and Subtract two extensible sparse matrices */
void esmat_add (Esmat* A, Esmat* B, Esmat* dest) 
{ esmat_bin_operate (A, B, dest, sum); }
void esmat_sub (Esmat* A, Esmat* B, Esmat* dest) 
{ esmat_bin_operate (A, B, dest, diff); }

/* min, max and sum over column elements */
void esmat_min_col (Esmat* A, Esmat* dest) 
{ esmat_operate_col (A, dest, min, 10e300); }
void esmat_max_col (Esmat* A, Esmat* dest) 
{ esmat_operate_col (A, dest, max, -10e300); }
void esmat_sum_col (Esmat* A, Esmat* dest) 
{ esmat_operate_col (A, dest, sum, 0.0); }

/* min, max and sum over row elements */
void esmat_min_row (Esmat* A, Esmat* dest) 
{ esmat_operate_row (A, dest, min, 10e300); }
void esmat_max_row (Esmat* A, Esmat* dest) 
{ esmat_operate_row (A, dest, max, -10e300); }
void esmat_sum_row (Esmat* A, Esmat* dest) 
{ esmat_operate_row (A, dest, sum, 0.0); }


