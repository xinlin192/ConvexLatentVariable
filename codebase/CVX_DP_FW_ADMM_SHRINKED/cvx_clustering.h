#include <cassert>
#include <queue>
#include <time.h>

#include "../util.h"
class Compare
{
    public:
        bool operator() (pair<int, double> obj1, pair<int, double> obj2)
        {
            return obj1.second > obj2.second;
        }
};

/* Compute the mutual distance of input instances contained within "data" */
void compute_dist_mat (vector<Instance*>& data, double ** dist_mat, int N, int D, dist_func df, bool isSym) {
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            Instance * xi = data[i];
            Instance * muj = data[j];
            dist_mat[i][j] = df (xi, muj, D);
        }
    }
}
