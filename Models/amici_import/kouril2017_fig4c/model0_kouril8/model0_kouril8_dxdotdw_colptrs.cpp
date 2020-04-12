#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_kouril8(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 4;
    colptrs[2] = 8;
    colptrs[3] = 10;
    colptrs[4] = 14;
    colptrs[5] = 17;
    colptrs[6] = 18;
}