#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_becker1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 5;
    colptrs[4] = 8;
    colptrs[5] = 10;
    colptrs[6] = 13;
    colptrs[7] = 15;
    colptrs[8] = 17;
}