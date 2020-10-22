#include "sundials/sundials_types.h"

void JSparseB_colptrs_model1_perelson2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 6;
    colptrs[3] = 8;
    colptrs[4] = 12;
    colptrs[5] = 15;
}