#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_becker2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 3;
    colptrs[3] = 5;
    colptrs[4] = 8;
    colptrs[5] = 8;
    colptrs[6] = 8;
}