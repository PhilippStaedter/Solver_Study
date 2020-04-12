#include "sundials/sundials_types.h"

void dxdotdw_colptrs_Pritchard2002(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 5;
    colptrs[3] = 7;
    colptrs[4] = 11;
    colptrs[5] = 14;
    colptrs[6] = 16;
    colptrs[7] = 20;
    colptrs[8] = 24;
    colptrs[9] = 26;
    colptrs[10] = 28;
    colptrs[11] = 32;
    colptrs[12] = 34;
    colptrs[13] = 37;
    colptrs[14] = 39;
    colptrs[15] = 42;
    colptrs[16] = 45;
    colptrs[17] = 48;
    colptrs[18] = 51;
    colptrs[19] = 54;
}