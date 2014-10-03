
#include "cvx_clustering.h"
#include "time.h"
#define INF_INT 60000
/* algorithmic options */ 

/* dumping options */
// #define FRANK_WOLFE_DUMP
// #define EXACT_LINE_SEARCH_DUMP
// #define BLOCKWISE_DUMP
// #define SUBPROBLEM_DUMP 

void blockwise_closed_form (double ** ytwo, double ** ztwo, double ** wtwo, double rho, double lambda, int N) {

    // STEP ONE: compute the optimal solution for truncated problem
    double ** wbar = mat_init (N, N);
    mat_zeros (wbar, N, N);
    mat_dot (rho, ztwo, wbar, N, N); // wbar = rho * z_2
    mat_sub (wbar, ytwo, wbar, N, N); // wbar = rho * z_2 - y_2
    mat_dot (1.0/rho, wbar, wbar, N, N); // wbar = (rho * z_2 - y_2) / rho

    // STEP TWO: find the closed-form solution for second subproblem
    for (int j = 0; j < N; j ++) {
        // 1. bifurcate the set of values
        vector< pair<int,double> > alpha_vec;
        for (int i = 0; i < N; i ++) {
            double value = wbar[i][j];
            alpha_vec.push_back (make_pair(i, abs(value)));
        }
        // 2. sorting
        std::sort (alpha_vec.begin(), alpha_vec.end(), pair_Second_Elem_Comparator);

        // 3. find mstar
        int mstar = 0; // number of elements support the sky
        double separator;
        double max_term = -INF, new_term;
        double sum_alpha = 0.0;
        for (int i = 0; i < N; i ++) {
            sum_alpha += alpha_vec[i].second;
            new_term = (sum_alpha - lambda) / (i + 1.0);
            if ( new_term > max_term ) {
                separator = alpha_vec[i].second;
                max_term = new_term;
                mstar = i + 1;
            }
        }
        cout << "mstar: " << mstar << ", max_term: " << max_term << endl;
        // 4. assign closed-form solution to wtwo
        if( max_term < 0 ){
            for(int i=0;i<N;i++)
                wtwo[i][j] = 0.0;
            continue;
        }
        for (int i = 0; i < N; i ++) {
            // harness vector of pair
            double value = wbar[i][j];
            if ( abs(value) >= separator ) {
                wtwo[i][j] = max_term;
            } else {
                // its ranking is above m*, directly inherit the wbar
                wtwo[i][j] = max(wbar[i][j],0.0);
            }
        }
    }
    // report the #iter and objective function
    /*cout << "[Blockwise] second_subproblem_obj: " << penalty << endl;
      cout << endl;*/

    // STEP THREE: recollect temporary variable - wbar
    mat_free (wbar, N, N);
}

// entry main function
int main (int argc, char ** argv) {
    // exception control: illustrate the usage if get input of wrong format
    if (argc < 6) {
        cerr << "Usage: gs_test [y_input] [z_input] [w_target] [N] [lambda]" << endl;
        cerr << "Note: dataFile must be scaled to [0,1] in advance." << endl;
        exit(-1);
    }

    // parse arguments
    char * yFile = argv[1];
    char * zFile = argv[2];
    char * wFile = argv[3];
    int N = atoi(argv[4]);
    double lambda = atof(argv[5]);

    cerr << "N = " << N << endl;

    double ** y_mat = mat_read (yFile, N, N);
    double ** z_mat = mat_read (zFile, N, N);
    Esmat* y_esmat = mat2esmat (y_mat, N, N);
    Esmat* z_esmat = mat2esmat (z_mat, N, N);

    // double ** w_target = mat_read (wFile, N, N);
    // Esmat* w_target_esmat = mat2esmat (w_target, N, N);
    // esmat_trim (w_target_esmat, 1e-6);
    // ofstream wtargetout ("gs_w2target");
    // wtargetout << esmat_toString (w_target_esmat); 
    // wtargetout.close();

    const int multiple = 1000;
    double ** w_mat = mat_init (N, N);
    clock_t a = clock();
    // /*
    for (int i = 0; i < multiple; i ++) {
    blockwise_closed_form (y_mat, z_mat, w_mat, 1.0, lambda, N);
    }
    // */
    blockwise_closed_form (y_mat, z_mat, w_mat, 1.0, lambda, N);
    cerr << "original: " << clock() - a << endl;
    Esmat* w_target_esmat = mat2esmat (w_mat, N, N);
    esmat_trim (w_target_esmat, 1e-10);
    ofstream wtargetout ("gs_w2target");
    wtargetout << esmat_toString (w_target_esmat); 
    wtargetout.close();

    Esmat* w_esmat = esmat_init (N, N);
    clock_t b = clock ();
    // /*
    for (int i = 0; i < multiple; i ++) {
        group_lasso_solver (y_esmat, z_esmat, w_esmat, 1.0, lambda);
    }
    // */
    group_lasso_solver (y_esmat, z_esmat, w_esmat, 1.0, lambda);
    cerr << "exsparse: " << clock() - b << endl;
    esmat_trim (w_esmat, 1e-10);
    ofstream wout ("gs_w2out");
    wout << esmat_toString (w_esmat); 
    wout.close();

    esmat_free (w_esmat);
    esmat_free (z_esmat);
    esmat_free (y_esmat);
    mat_free (z_mat, N, N);
    mat_free (y_mat, N, N);
}
