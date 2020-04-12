#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_laub1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 4;
    colptrs[3] = 7;
    colptrs[4] = 9;
    colptrs[5] = 12;
    colptrs[6] = 14;
    colptrs[7] = 17;
}