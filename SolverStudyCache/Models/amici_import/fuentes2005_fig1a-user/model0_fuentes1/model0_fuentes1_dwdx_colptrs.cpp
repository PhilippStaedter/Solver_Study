#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_fuentes1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 3;
    colptrs[3] = 3;
    colptrs[4] = 5;
}