#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_bier2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 5;
    colptrs[4] = 8;
}