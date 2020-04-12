#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_esteban2013(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 3;
    colptrs[4] = 4;
    colptrs[5] = 5;
    colptrs[6] = 6;
    colptrs[7] = 8;
    colptrs[8] = 9;
    colptrs[9] = 10;
    colptrs[10] = 11;
    colptrs[11] = 12;
    colptrs[12] = 14;
    colptrs[13] = 15;
}