#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_miao2008(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 4;
    colptrs[2] = 8;
    colptrs[3] = 10;
    colptrs[4] = 14;
}