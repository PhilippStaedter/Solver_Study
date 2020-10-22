#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_saeidi1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 4;
    colptrs[4] = 6;
    colptrs[5] = 8;
    colptrs[6] = 8;
    colptrs[7] = 10;
    colptrs[8] = 12;
    colptrs[9] = 13;
    colptrs[10] = 13;
}