#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_karin5(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 7;
    colptrs[3] = 9;
    colptrs[4] = 11;
    colptrs[5] = 14;
}