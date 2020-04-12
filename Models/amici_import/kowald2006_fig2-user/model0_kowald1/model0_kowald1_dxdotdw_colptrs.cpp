#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_kowald1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 3;
    colptrs[3] = 4;
    colptrs[4] = 7;
    colptrs[5] = 8;
    colptrs[6] = 9;
    colptrs[7] = 11;
    colptrs[8] = 14;
    colptrs[9] = 15;
    colptrs[10] = 18;
    colptrs[11] = 21;
    colptrs[12] = 23;
    colptrs[13] = 24;
    colptrs[14] = 25;
    colptrs[15] = 27;
    colptrs[16] = 29;
    colptrs[17] = 30;
}