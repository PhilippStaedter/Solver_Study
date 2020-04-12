#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_naresh2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 5;
    colptrs[2] = 5;
    colptrs[3] = 12;
    colptrs[4] = 18;
    colptrs[5] = 22;
}