#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_borghans3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 8;
    colptrs[5] = 10;
    colptrs[6] = 12;
    colptrs[7] = 14;
    colptrs[8] = 16;
}