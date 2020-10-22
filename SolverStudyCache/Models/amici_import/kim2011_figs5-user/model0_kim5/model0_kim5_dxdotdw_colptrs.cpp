#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_kim5(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 3;
    colptrs[4] = 4;
    colptrs[5] = 5;
    colptrs[6] = 6;
    colptrs[7] = 7;
    colptrs[8] = 8;
    colptrs[9] = 9;
}