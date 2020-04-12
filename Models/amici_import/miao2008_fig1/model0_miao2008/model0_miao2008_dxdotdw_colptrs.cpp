#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_miao2008(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 3;
    colptrs[3] = 5;
    colptrs[4] = 9;
    colptrs[5] = 10;
    colptrs[6] = 12;
    colptrs[7] = 13;
    colptrs[8] = 15;
    colptrs[9] = 16;
}