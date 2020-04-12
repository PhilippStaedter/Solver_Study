#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_kouril8(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 3;
    colptrs[3] = 6;
    colptrs[4] = 6;
    colptrs[5] = 7;
    colptrs[6] = 9;
    colptrs[7] = 10;
    colptrs[8] = 12;
    colptrs[9] = 13;
    colptrs[10] = 14;
    colptrs[11] = 15;
    colptrs[12] = 15;
    colptrs[13] = 15;
}