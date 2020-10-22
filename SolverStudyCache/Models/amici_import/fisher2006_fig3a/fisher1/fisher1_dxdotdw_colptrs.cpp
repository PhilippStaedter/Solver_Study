#include "sundials/sundials_types.h"

void dxdotdw_colptrs_fisher1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 4;
    colptrs[3] = 7;
    colptrs[4] = 9;
    colptrs[5] = 12;
    colptrs[6] = 15;
    colptrs[7] = 17;
    colptrs[8] = 19;
    colptrs[9] = 21;
    colptrs[10] = 24;
    colptrs[11] = 26;
    colptrs[12] = 28;
    colptrs[13] = 30;
    colptrs[14] = 33;
    colptrs[15] = 35;
    colptrs[16] = 37;
    colptrs[17] = 40;
}