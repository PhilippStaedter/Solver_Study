#include "sundials/sundials_types.h"

void JSparseB_colptrs_model0_marhl(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 4;
    colptrs[3] = 7;
    colptrs[4] = 12;
    colptrs[5] = 15;
}