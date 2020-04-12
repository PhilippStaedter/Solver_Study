#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_borghans2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 8;
    colptrs[5] = 9;
    colptrs[6] = 10;
    colptrs[7] = 11;
    colptrs[8] = 13;
}