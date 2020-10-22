#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_zi1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 7;
    colptrs[4] = 9;
    colptrs[5] = 10;
    colptrs[6] = 12;
    colptrs[7] = 13;
    colptrs[8] = 14;
    colptrs[9] = 16;
    colptrs[10] = 17;
    colptrs[11] = 19;
    colptrs[12] = 22;
    colptrs[13] = 23;
    colptrs[14] = 25;
    colptrs[15] = 28;
    colptrs[16] = 29;
}