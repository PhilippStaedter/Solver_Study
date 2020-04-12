#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_bachar1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 3;
    colptrs[3] = 4;
    colptrs[4] = 5;
    colptrs[5] = 7;
    colptrs[6] = 9;
    colptrs[7] = 10;
    colptrs[8] = 11;
}