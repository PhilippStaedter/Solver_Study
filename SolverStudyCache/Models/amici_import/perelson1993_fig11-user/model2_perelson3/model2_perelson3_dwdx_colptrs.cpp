#include "sundials/sundials_types.h"

void dwdx_colptrs_model2_perelson3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 5;
    colptrs[3] = 7;
    colptrs[4] = 8;
    colptrs[5] = 11;
}