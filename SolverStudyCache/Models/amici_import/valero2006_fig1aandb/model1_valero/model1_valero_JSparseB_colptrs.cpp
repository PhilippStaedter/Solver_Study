#include "sundials/sundials_types.h"

void JSparseB_colptrs_model1_valero(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 5;
    colptrs[3] = 8;
    colptrs[4] = 8;
    colptrs[5] = 9;
    colptrs[6] = 11;
}