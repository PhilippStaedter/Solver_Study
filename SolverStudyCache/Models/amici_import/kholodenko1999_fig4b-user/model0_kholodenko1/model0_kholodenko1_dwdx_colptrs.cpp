#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_kholodenko1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 5;
    colptrs[3] = 9;
    colptrs[4] = 10;
    colptrs[5] = 13;
    colptrs[6] = 14;
    colptrs[7] = 15;
    colptrs[8] = 17;
    colptrs[9] = 19;
    colptrs[10] = 21;
    colptrs[11] = 31;
    colptrs[12] = 33;
    colptrs[13] = 35;
    colptrs[14] = 37;
    colptrs[15] = 40;
    colptrs[16] = 43;
    colptrs[17] = 47;
    colptrs[18] = 49;
    colptrs[19] = 53;
    colptrs[20] = 56;
    colptrs[21] = 59;
    colptrs[22] = 63;
    colptrs[23] = 64;
}