/*###############################################################
## MODULE: mat.cpp
## VERSION: 1.0 
## SINCE 2014-06-16
## AUTHOR Jimmy Lin (xl5224) - JimmyLin@utexas.edu  
## DESCRIPTION: 
##    Matrix data structure and associated operations.
#################################################################
## Edited by MacVim
## Class Info auto-generated by Snippet 
################################################################*/

double ** mat_init (int nRows, int nCols) {

    double ** res = new double * [nRows];
    for (int i = 0; i < nRows; i ++) {
        res[i] = new double [nCols];
    }

    return res;
}

void mat_free (double ** src, int nRows, int nCols) {

    for (int i = 0; i < nRows; i ++) {
        delete src[i];
    }

}

void mat_add (double ** src1, double ** src2, double ** dest, int nRows, int nCols) {

    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            dest[i][j] = src1[i][j] + src2[i][j];
        }
    }

}

void mat_sub (double ** src1, double ** src2, double ** dest, int nRows, int nCols) {
    
    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            dest[i][j] = src1[i][j] - src2[i][j];
        }
    }

}

void mat_times (double ** src1, double ** src2, double ** dest, int nRows, int nCols) {

    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            dest[i][j] = src1[i][j]* src2[i][j];
        }
    }

}

void mat_dot (double scalar, double ** src, double ** dest, int nRows, int nCols) {

    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            dest[i][j] = scalar * src[i][j];
        }
    }

}

// tranpose first matrix and then compute dot product 
void mat_tdot (double ** src1, double ** src2, double ** dest, int nRows, int nCols) {
    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            double temp;
            for (int k = 0; k < nRows; k ++) {
                temp += src1[k][i] * src2[k][i];
            }
            dest[i][j] = temp;
        }
    }
}

void mat_print (double ** src, int nRows, int nCols) {

    string field_seperator = ",";
    string line_separator = "\n";
    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) { 
            cout << src[i][j] << field_seperator;
        }
        cout << line_separator;
    }

}

void mat_copy (double ** src, double ** dest, int nRows, int nCols) {

    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            dest[i][j] = src[i][j];
        }
    }

}

double mat_sum (double ** src, int nRows, int nCols) {

    double sum = 0.0;
    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            sum += src[i][j];
        }
    }

    return sum;
}

void mat_max_col (double ** src, double ** dest, int nRows, int nCols) {

    // we assume that the given dest is all-zero mat

    for (int j = 0; j < nCols; j ++) {
        int max_index = -1;
        int max_value = -INF;
        for (int i = 0; i < nRows; i ++) {
            if (src[i][j] > max_value) {
                max_index = i;
                max_value = src[i][j];
            }
        }
        dest[max_index][j] = 1;
    }

}

void mat_min_col (double ** src, double ** dest, int nRows, int nCols) {

    // we assume that the given dest is all-zero mat

    for (int j = 0; j < nCols; j ++) {
        int min_index = -1;
        int min_value = INF;
        for (int i = 0; i < nRows; i ++) {
            if (src[i][j] < min_value) {
                min_index = i;
                min_value = src[i][j];
            }
        }
        dest[min_index][j] = 1;
    }

}

double mat_norm2 (double ** src, int nRows, int nCols) {
    double sum = 0.0;
    for (int i = 0; i < nRows; i ++) {
        for (int j = 0; j < nCols; j ++) {
            sum += src[i][j] * src[i][j];
        }
    }   
    return sum;
}

// TODO: mat_write and mat_read
