#include "sundials/sundials_types.h"

void JSparseB_colptrs_model2_balagadde1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 9;
    colptrs[4] = 12;
    colptrs[5] = 12;
    colptrs[6] = 12;
    colptrs[7] = 12;
}