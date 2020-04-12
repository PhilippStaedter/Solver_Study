#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_lee4(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 6;
    colptrs[2] = 8;
    colptrs[3] = 10;
    colptrs[4] = 10;
    colptrs[5] = 12;
    colptrs[6] = 14;
    colptrs[7] = 17;
    colptrs[8] = 20;
    colptrs[9] = 22;
    colptrs[10] = 23;
    colptrs[11] = 25;
    colptrs[12] = 27;
    colptrs[13] = 28;
    colptrs[14] = 28;
}