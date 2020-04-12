#include "sundials/sundials_types.h"

void dwdx_colptrs_model1_tripathi1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 3;
    colptrs[3] = 7;
    colptrs[4] = 10;
    colptrs[5] = 12;
}