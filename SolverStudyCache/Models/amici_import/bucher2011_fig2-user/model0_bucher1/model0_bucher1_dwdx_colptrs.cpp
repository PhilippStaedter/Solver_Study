#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_bucher1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 8;
    colptrs[3] = 10;
    colptrs[4] = 11;
    colptrs[5] = 14;
    colptrs[6] = 16;
    colptrs[7] = 17;
    colptrs[8] = 20;
    colptrs[9] = 22;
    colptrs[10] = 23;
    colptrs[11] = 30;
    colptrs[12] = 31;
    colptrs[13] = 32;
    colptrs[14] = 34;
    colptrs[15] = 35;
    colptrs[16] = 36;
    colptrs[17] = 38;
    colptrs[18] = 39;
}