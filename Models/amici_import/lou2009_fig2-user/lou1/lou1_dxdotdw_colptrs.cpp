#include "sundials/sundials_types.h"

void dxdotdw_colptrs_lou1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 3;
    colptrs[3] = 4;
    colptrs[4] = 5;
    colptrs[5] = 7;
    colptrs[6] = 8;
    colptrs[7] = 9;
    colptrs[8] = 10;
    colptrs[9] = 12;
    colptrs[10] = 13;
    colptrs[11] = 14;
    colptrs[12] = 15;
}