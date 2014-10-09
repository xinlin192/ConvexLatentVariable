using namespace std;
#include "DP_MEDOIDS.h"
#include "../util.h"
#include <cassert>
#include <cmath>
#define INTEGER_MAX 3000000

typedef double (* dist_func) (Instance*, Instance*, int); 
