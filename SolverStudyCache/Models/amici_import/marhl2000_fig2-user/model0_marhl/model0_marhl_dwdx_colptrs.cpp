#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_marhl(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 3;
    colptrs[3] = 4;
    colptrs[4] = 10;
    colptrs[5] = 11;
}