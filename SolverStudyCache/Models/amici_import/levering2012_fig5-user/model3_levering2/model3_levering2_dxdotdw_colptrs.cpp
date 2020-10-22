#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model3_levering2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 4;
    colptrs[2] = 10;
    colptrs[3] = 15;
    colptrs[4] = 18;
    colptrs[5] = 20;
    colptrs[6] = 24;
    colptrs[7] = 26;
    colptrs[8] = 28;
    colptrs[9] = 32;
    colptrs[10] = 35;
    colptrs[11] = 39;
    colptrs[12] = 43;
    colptrs[13] = 45;
    colptrs[14] = 49;
    colptrs[15] = 51;
    colptrs[16] = 54;
    colptrs[17] = 59;
    colptrs[18] = 63;
    colptrs[19] = 67;
    colptrs[20] = 71;
    colptrs[21] = 75;
}