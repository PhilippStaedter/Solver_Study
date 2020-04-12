#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_bier2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 7;
    colptrs[5] = 8;
    colptrs[6] = 10;
}