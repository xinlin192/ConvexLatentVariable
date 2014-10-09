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

