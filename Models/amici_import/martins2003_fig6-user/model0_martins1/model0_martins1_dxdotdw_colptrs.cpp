#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_martins1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 8;
    colptrs[4] = 10;
    colptrs[5] = 12;
    colptrs[6] = 15;
    colptrs[7] = 19;
    colptrs[8] = 22;
    colptrs[9] = 24;
    colptrs[10] = 27;
    colptrs[11] = 30;
    colptrs[12] = 32;
    colptrs[13] = 34;
    colptrs[14] = 37;
    colptrs[15] = 39;
    colptrs[16] = 41;
}