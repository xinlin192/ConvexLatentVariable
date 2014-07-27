/*###############################################################
## MODULE: exSparseMat.cpp
## VERSION: 1.0 
## SINCE 2014-07-25
## AUTHOR Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
## DESCRIPTION: 
##      Test script for Extensible Sparse Matrix
#################################################################
## Edited by MacVim
## Class Info auto-generated by Snippet 
################################################################*/

#include "exSparseMat.h"

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
};
/* Deallocate given Esmat */
void esmat_free (Esmat* src) {
    src->val.clear();
    delete src;
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
void esmat_submat_row (int start_index, int end_index, Esmat* mat, Esmat* submat) {
    assert (start_index < end_index);

    // renew the characteristics of submat
    submat->nRows = end_index - start_index;
    submat->nCols = mat->nCols;
    submat->val.clear();

    int mat_size = mat->val.size();
    for (int i = 0; i < mat_size; i ++) {
        int mat_index = mat->val[i].first;
        int mat_row_index = mat_index % mat->nRows;
        int mat_col_index = mat_index / mat->nRows;
        // if element of that entry is in required submat
        if (start_index <= mat_row_index && mat_row_index < end_index) {
            // get value of that entry
            int value = mat->val[i].second;
            // compute corresponding position in submat
            int submat_row_index = mat_row_index - start_index;
            int submat_col_index = mat_col_index;
            int submat_index = submat_row_index + submat_col_index * submat->nRows;
            submat->val.push_back(make_pair(submat_index, value));
        }
    }
}
/* put submat to specified position of mat */
void esmat_merge_row (Esmat* submat, int start_index, int end_index, Esmat* mat) {
    assert (start_index < end_index);

    int submat_size = submat->val.size();
    for (int i = 0; i < submat_size; i ++) {
        int submat_index = submat->val[i].first;
        int submat_row_index = submat_index % submat->nRows;
        int submat_col_index = submat_index / submat->nRows;
        // get value of that entry
        int value = submat->val[i].second;
        // compute corresponding position in submat
        int mat_row_index = submat_row_index + start_index;
        int mat_col_index = submat_col_index;
        int mat_index = mat_row_index + mat_col_index * mat->nRows;
        mat->val.push_back(make_pair(submat_index, value));
    }
    // realign the mat->val with index-increasing order
    esmat_align (mat);
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

/* regularization for penalizing number of vocabularies used by topics */
/*
mat_nonzero_index_col {
}
*/

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
void esmat_operate_row (Esmat* A, Esmat* dest, Operation opt) {

    int sizeA = A->val.size();

    // set up dest esmat, here assume it has been initialized
    dest->nRows = A->nRows;
    dest->nCols = 1;
    dest->val.clear();

    double init_value = 0.0; 
    /* TODO
    if (opt == *min) init_value = INF;
    else if (opt == *max) init_value = -INF;
    */
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
void esmat_operate_col (Esmat* A, Esmat* dest, Operation opt) {
    
    int sizeA = A->val.size();

    // set up dest esmat, here assume it has been initialized
    dest->nRows = 1;
    dest->nCols = A->nCols;
    dest->val.clear();

    double init_value = 0.0;
    /* TODO
    if (opt == *min) init_value = INF;
    else if (opt == *max) init_value = -INF;
    */
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

