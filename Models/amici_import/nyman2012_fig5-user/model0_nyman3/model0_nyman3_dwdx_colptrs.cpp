#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_nyman3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 4;
    colptrs[4] = 5;
    colptrs[5] = 7;
    colptrs[6] = 9;
    colptrs[7] = 12;
    colptrs[8] = 13;
    colptrs[9] = 16;
}