#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_borghans1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 0;
    colptrs[2] = 0;
    colptrs[3] = 3;
    colptrs[4] = 5;
    colptrs[5] = 9;
}