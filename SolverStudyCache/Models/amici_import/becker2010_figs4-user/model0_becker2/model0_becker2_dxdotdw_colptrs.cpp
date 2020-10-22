#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_becker2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 5;
    colptrs[4] = 8;
    colptrs[5] = 10;
    colptrs[6] = 12;
    colptrs[7] = 14;
    colptrs[8] = 16;
}