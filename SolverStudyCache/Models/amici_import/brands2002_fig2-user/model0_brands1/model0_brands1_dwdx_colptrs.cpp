#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_brands1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 1;
    colptrs[3] = 3;
    colptrs[4] = 3;
    colptrs[5] = 3;
    colptrs[6] = 3;
    colptrs[7] = 7;
    colptrs[8] = 10;
    colptrs[9] = 10;
    colptrs[10] = 11;
    colptrs[11] = 13;
}