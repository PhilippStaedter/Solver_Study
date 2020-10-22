#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_lee2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 4;
    colptrs[2] = 6;
    colptrs[3] = 8;
    colptrs[4] = 10;
    colptrs[5] = 12;
    colptrs[6] = 13;
    colptrs[7] = 15;
    colptrs[8] = 16;
    colptrs[9] = 16;
}