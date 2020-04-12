#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_gorlich1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 4;
    colptrs[4] = 6;
    colptrs[5] = 8;
    colptrs[6] = 10;
    colptrs[7] = 11;
    colptrs[8] = 13;
    colptrs[9] = 14;
    colptrs[10] = 16;
    colptrs[11] = 18;
    colptrs[12] = 21;
    colptrs[13] = 23;
}