#include "sundials/sundials_types.h"

void dxdotdw_colptrs_gardner1_Fig4B(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 3;
    colptrs[3] = 5;
    colptrs[4] = 6;
    colptrs[5] = 7;
    colptrs[6] = 8;
    colptrs[7] = 9;
    colptrs[8] = 10;
    colptrs[9] = 11;
    colptrs[10] = 12;
    colptrs[11] = 13;
    colptrs[12] = 16;
    colptrs[13] = 19;
}