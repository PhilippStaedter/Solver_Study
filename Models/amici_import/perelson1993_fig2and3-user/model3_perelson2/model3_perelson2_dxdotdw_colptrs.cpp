#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model3_perelson2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 7;
    colptrs[5] = 9;
    colptrs[6] = 11;
    colptrs[7] = 13;
    colptrs[8] = 15;
    colptrs[9] = 17;
}