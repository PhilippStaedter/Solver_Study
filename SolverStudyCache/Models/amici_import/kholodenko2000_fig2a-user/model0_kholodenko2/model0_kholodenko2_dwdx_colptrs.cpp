#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_kholodenko2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 3;
    colptrs[3] = 5;
    colptrs[4] = 6;
    colptrs[5] = 7;
    colptrs[6] = 10;
    colptrs[7] = 12;
    colptrs[8] = 15;
}