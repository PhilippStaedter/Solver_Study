#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_perelson1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 2;
    colptrs[3] = 4;
    colptrs[4] = 5;
}