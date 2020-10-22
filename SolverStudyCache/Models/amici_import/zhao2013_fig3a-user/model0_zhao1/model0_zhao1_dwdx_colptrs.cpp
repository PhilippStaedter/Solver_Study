#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_zhao1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 0;
    colptrs[2] = 4;
    colptrs[3] = 8;
    colptrs[4] = 11;
    colptrs[5] = 14;
}