#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_hald(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 5;
    colptrs[3] = 8;
    colptrs[4] = 11;
    colptrs[5] = 15;
    colptrs[6] = 17;
    colptrs[7] = 20;
    colptrs[8] = 24;
    colptrs[9] = 26;
    colptrs[10] = 27;
    colptrs[11] = 31;
    colptrs[12] = 34;
    colptrs[13] = 36;
    colptrs[14] = 40;
    colptrs[15] = 42;
    colptrs[16] = 46;
    colptrs[17] = 49;
    colptrs[18] = 50;
    colptrs[19] = 54;
    colptrs[20] = 56;
    colptrs[21] = 58;
    colptrs[22] = 60;
    colptrs[23] = 62;
    colptrs[24] = 64;
    colptrs[25] = 66;
    colptrs[26] = 70;
    colptrs[27] = 74;
}