#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_fraser1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 6;
    colptrs[3] = 9;
    colptrs[4] = 13;
    colptrs[5] = 17;
    colptrs[6] = 21;
}