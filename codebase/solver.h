#include "util.h"


void frank_wolfe_solver (double ** dist_mat, double ** y, double ** z, double ** w, double rho, int R, int C, int FW_MAX_ITER, set<int>& col_active_set) {
    // cout << "within frank_wolfe_solver" << endl;
    // STEP ONE: compute gradient mat initially
    vector< set< pair<int, double> > > actives (R, set<pair<int,double> >());
    vector< priority_queue< pair<int,double>, vector< pair<int,double> >, Int_Double_Pair_Dec> > pqueues (R, priority_queue< pair<int,double>, vector< pair<int,double> >, Int_Double_Pair_Dec> ());
    for (int i = 0; i < R; i++) {
        for (set<int>::iterator it=col_active_set.begin();it != col_active_set.end(); ++it) {
            int j = *it;
            double grad=0.5*dist_mat[i][j]+y[i][j]+rho*(w[i][j]-z[i][j]); 
            if (w[i][j] > 1e-10) 
                actives[i].insert(make_pair(j,grad));
            else 
                pqueues[i].push(make_pair(j, grad));
        }
    }
    // STEP TWO: iteration solve each row 
    int k = 0;  // iteration number
    vector<bool> is_fw_opt_reached (R, false);
    set<pair<int,double> >::iterator it;
    // cout << "within frank_wolfe_solver: start iteration" << endl;
    while (k < FW_MAX_ITER) { // TODO: change to use portional criteria
        // compute new active atom: can be in active set or not
        vector< pair<int, double> > s (R, pair<int,double>());
        vector<bool> isInActives (R, false);
        for (int i = 0; i < R; i++) {
            if (is_fw_opt_reached[i]) continue;
            if (pqueues[i].size() <= 0) continue;
            pair<int,double> tmp_pair = pqueues[i].top();
            s[i].first = tmp_pair.first;
            s[i].second = tmp_pair.second;
            // cout << s[i].first << ":" << s[i].second << endl;
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) {
                // take the minimal of each row
                if (it->second < s[i].second) {
                    isInActives[i] = true;
                    s[i].first = it->first;
                    s[i].second = it->second;
                }
            }
            // compute gamma: inexact or exact
            double gamma; // step size of line search
#ifdef EXACT_LINE_SEARCH
            double sum1=0.0, sum2=0.0, sum3=0.0, sum4=0.0;
            // gamma* = (sum1 + sum2 + sum3) / sum4, where
            // sum1 = 1/2 sum_n sum_k (w - s)_nk * || x_n - mu_k ||^2
            // sum2 = sum_n sum_k y_nk (w - s)_nk
            // sum3 = - rho * sum_n sum_k  (w - z) (w-s)
            // sum4 = sum_n sum_k rho * (s - w)^2
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) {
                double w_minus_s;
                double w_minus_z = w[i][it->first] - z[i][it->first];
                if (it->first == s[i].first) {
                    w_minus_s = w[i][it->first]-1.0;
                } else {
                    w_minus_s = w[i][it->first];
                }
                sum1 += 0.5 * w_minus_s * (dist_mat[i][it->first] -r);
                sum2 += y[i][it->first] * w_minus_s;
                sum3 += rho * w_minus_s * w_minus_z;
                sum4 += rho * w_minus_s * w_minus_s; 
            }
            if (!isInActives[i]) {
                sum1 += 0.5 * (-1.0) * (dist_mat[i][s[i].first] - r);
                sum2 += y[i][it->first] * (-1.0);
                sum3 += rho * (-1.0) * (w[i][s[i].first]-z[i][s[i].first]);
                sum4 += rho;
            }

            if (fabs(sum4) > 0) {
                gamma = (sum1 + sum2 + sum3) / sum4;
#ifdef EXACT_LINE_SEARCH_DUMP
                cout << "[exact] i=" << i ;
                cout << ",k=" << k;
                cout << ",sum1="<< sum1;
                cout << ",sum2="<< sum2;
                cout << ",sum3="<< sum3;
                cout << ",sum4="<< sum4;
                cout << ",gamma="<< gamma;
                cout << endl;
#endif
                gamma = max(gamma, 0.0);
                gamma = min(gamma, 1.0);
            } else {
                gamma = 0.0;
                is_fw_opt_reached[i] = true;
            }
#else
            gamma = 2.0 / (k+2.0);
#endif
            // update w
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) 
                w[i][it->first] *= (1-gamma);
            w[i][s[i].first] += gamma;
            // update new actives 
            set< pair<int, double> > temp;
            if (!isInActives[i]) {
                actives[i].insert(pqueues[i].top());
                pqueues[i].pop();
            }
            double new_grad;
            for (it=actives[i].begin(); it!=actives[i].end(); ++it) {
                int j = it->first;
                new_grad=0.5*dist_mat[i][j]+y[i][j]+rho*(w[i][j]-z[i][j]); 
                temp.insert (make_pair(it->first, new_grad));
            }
            actives[i].swap(temp);
        }
        k ++;
    }
}

void skyline (double** wout, double**wbar, int R_start, int R_end, int C, double lambda, set<int> col_active_sets) {
    vector< vector< double > > alpha_vec (C, vector<double>());
    vector< int > num_alpha_elem (C, 0);
    for (int i = R_start; i < R_end; i ++) {
        set<int>::iterator it;
        for (it=col_active_sets.begin();it!=col_active_sets.end();++it) {
            int j = *it;
            if (wbar[i][j] > SPARSITY_TOL) {
                alpha_vec[j].push_back (abs(wbar[i][j]));
                ++ num_alpha_elem[j];
            }
        }
    }
    vector<double> max_term (C, -1e50);
    vector<double> separator (C, 0.0);
    int R = R_end - R_start;
    for (int j = 0; j < C; j ++) {
        if (num_alpha_elem[j] == 0) continue;
        // 2. sorting
        std::sort (alpha_vec[j].begin(), alpha_vec[j].end(), double_dec_comp);
        // 3. find mstar
        int mstar = 0; // number of elements support the sky
        double new_term, sum_alpha = 0.0;
        for (int i = 0; i < num_alpha_elem[j]; i ++) {
            sum_alpha += alpha_vec[j][i];
            new_term = (sum_alpha - lambda) / (i + 1.0);
            if ( new_term > max_term[j] ) {
                separator[j] = alpha_vec[j][i];
                max_term[j] = new_term;
                mstar = i;
            }
        }
        if (max_term[j] < 0) 
            max_term[j] = (sum_alpha - lambda) / R;
    }
    for (int i = R_start; i < R_end; i ++) {
        // 4. assign closed-form solution to wout
        set<int>::iterator it;
        for (it=col_active_sets.begin();it!=col_active_sets.end();++it) {
            int j = *it;
            if ( max_term[j] < 0 ) {
                wout[i][j] = 0.0;
                continue;
            }
            double wbar_val = wbar[i][j];
            if ( abs(wbar_val) >= separator[j] ) 
                wout[i][j] = max_term[j];
            else 
                // its ranking is above m*, directly inherit the wbar
                wout[i][j] = max(wbar_val,0.0);
        }
    }
}
void global_group_lasso_solver (double** y, double** z, double** w, double rho, vector<double>& lambda, Lookups *tables, set<int> col_active_sets) {
    int R = tables->nWords;
    int C = tables->nDocs;
    vector< pair<int,int> > word_lookup = *(tables->word_lookup);
    vector< pair<int,int> > doc_lookup = *(tables->doc_lookup);
    int global_lambda = lambda[0];
    int local_lambda = lambda[1];
    double** wbar = mat_init (R, C); mat_zeros (wbar, R, C);
    for (int i = 0; i < R; i ++) {
        for (int j = 0; j < C; j ++) {
            wbar[i][j] = (rho * z[i][j] - y[i][j]) / rho;
        }
    }
    skyline (w, wbar, 0, R, C, global_lambda, col_active_sets);
    mat_free (wbar, R, C);
}
void combined_group_lasso_solver (double** y, double** z, double** w, double rho, vector<double>& lambda, Lookups *tables, set<int> col_active_sets) {
    int R = tables->nWords;
    int C = tables->nDocs;
    vector< pair<int,int> > word_lookup = *(tables->word_lookup);
    vector< pair<int,int> > doc_lookup = *(tables->doc_lookup);
    int global_lambda = lambda[0];
    int local_lambda = lambda[1];

    double** wlocal = mat_init (R, C); mat_zeros (wlocal, R, C);
    double** wbar = mat_init (R, C); mat_zeros (wbar, R, C);
    for (int i = 0; i < R; i ++) {
        for (int j = 0; j < C; j ++) {
            wbar[i][j] = (rho * z[i][j] - y[i][j]) / rho;
        }
    }
    // extend the group lasso solver to both local and global
    for (int d = 0; d < tables->nDocs; d++) {
        int R_start = doc_lookup[d].first;
        int R_end = doc_lookup[d].second;
        skyline (wlocal, wbar, R_start, R_end, C, local_lambda, col_active_sets);
    }
    skyline (w, wlocal, 0, R, C, global_lambda, col_active_sets);

    mat_free (wlocal, R, C);
    mat_free (wbar, R, C);
}


