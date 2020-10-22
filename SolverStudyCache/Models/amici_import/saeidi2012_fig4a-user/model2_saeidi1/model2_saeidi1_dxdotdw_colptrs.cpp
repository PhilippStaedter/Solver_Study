#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model2_saeidi1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 8;
    colptrs[5] = 11;
    colptrs[6] = 14;
    colptrs[7] = 16;
}