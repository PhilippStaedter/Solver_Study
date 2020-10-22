#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_lee4(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 5;
    colptrs[3] = 8;
    colptrs[4] = 11;
    colptrs[5] = 14;
    colptrs[6] = 17;
    colptrs[7] = 19;
    colptrs[8] = 21;
    colptrs[9] = 24;
    colptrs[10] = 26;
    colptrs[11] = 29;
    colptrs[12] = 32;
    colptrs[13] = 35;
    colptrs[14] = 38;
    colptrs[15] = 41;
    colptrs[16] = 44;
}