#include "sundials/sundials_types.h"

void JSparseB_colptrs_model0_ma(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 5;
    colptrs[3] = 8;
    colptrs[4] = 10;
    colptrs[5] = 12;
    colptrs[6] = 14;
    colptrs[7] = 17;
}