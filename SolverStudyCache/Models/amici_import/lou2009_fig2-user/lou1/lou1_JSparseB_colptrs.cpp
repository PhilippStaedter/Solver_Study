#include "sundials/sundials_types.h"

void JSparseB_colptrs_lou1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 6;
    colptrs[2] = 12;
    colptrs[3] = 18;
    colptrs[4] = 23;
    colptrs[5] = 28;
    colptrs[6] = 34;
}