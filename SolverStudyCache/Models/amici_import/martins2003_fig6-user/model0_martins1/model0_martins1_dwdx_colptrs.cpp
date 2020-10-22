#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_martins1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 0;
    colptrs[2] = 2;
    colptrs[3] = 5;
    colptrs[4] = 8;
    colptrs[5] = 10;
    colptrs[6] = 10;
    colptrs[7] = 10;
    colptrs[8] = 11;
    colptrs[9] = 12;
    colptrs[10] = 12;
    colptrs[11] = 13;
    colptrs[12] = 13;
    colptrs[13] = 15;
    colptrs[14] = 17;
}