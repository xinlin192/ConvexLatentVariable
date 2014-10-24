#include <cassert>
#include <queue>
#include <time.h>
#include <omp.h>

#include "../util.h"
class Compare
{
    public:
        bool operator() (pair<int, double> obj1, pair<int, double> obj2)
        {
            return obj1.second > obj2.second;
        }
};

// reg = lambda * sum_k max_n | w_nk |  -> group-lasso
double get_reg (double ** z, int N, double lambda) {
    double * maxn = new double [N]; 
    for (int i = 0; i < N; i ++) 
        maxn[i] = -INF;
    for (int i = 0; i < N; i ++) 
        for (int j = 0; j < N; j ++) 
            if (z[i][j] > maxn[j])
                maxn[j] = z[i][j];
    delete[] maxn;
    double sumk = 0.0;
    for (int i = 0; i < N; i ++) 
        sumk += lambda*maxn[i];
    double group_lasso = sumk; 
    return group_lasso;
}
