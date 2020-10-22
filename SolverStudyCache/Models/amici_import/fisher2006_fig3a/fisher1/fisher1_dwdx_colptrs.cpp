#include "sundials/sundials_types.h"

void dwdx_colptrs_fisher1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 4;
    colptrs[2] = 8;
    colptrs[3] = 10;
    colptrs[4] = 12;
    colptrs[5] = 14;
    colptrs[6] = 16;
    colptrs[7] = 19;
    colptrs[8] = 22;
    colptrs[9] = 25;
    colptrs[10] = 28;
    colptrs[11] = 31;
    colptrs[12] = 34;
    colptrs[13] = 37;
    colptrs[14] = 40;
}