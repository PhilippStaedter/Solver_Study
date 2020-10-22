#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_li1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 6;
    colptrs[2] = 11;
    colptrs[3] = 14;
}