#include "sundials/sundials_types.h"

void dxdotdw_colptrs_Leber2015(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 3;
    colptrs[4] = 4;
    colptrs[5] = 5;
    colptrs[6] = 7;
    colptrs[7] = 9;
    colptrs[8] = 11;
    colptrs[9] = 12;
    colptrs[10] = 13;
    colptrs[11] = 14;
    colptrs[12] = 16;
    colptrs[13] = 18;
    colptrs[14] = 20;
    colptrs[15] = 22;
    colptrs[16] = 24;
    colptrs[17] = 26;
    colptrs[18] = 27;
    colptrs[19] = 28;
    colptrs[20] = 30;
    colptrs[21] = 32;
    colptrs[22] = 34;
    colptrs[23] = 36;
    colptrs[24] = 37;
    colptrs[25] = 39;
    colptrs[26] = 41;
    colptrs[27] = 43;
    colptrs[28] = 44;
    colptrs[29] = 45;
    colptrs[30] = 47;
}