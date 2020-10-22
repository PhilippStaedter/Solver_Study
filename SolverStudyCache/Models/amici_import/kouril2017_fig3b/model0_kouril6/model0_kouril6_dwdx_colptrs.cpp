#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_kouril6(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 3;
    colptrs[3] = 5;
    colptrs[4] = 6;
    colptrs[5] = 7;
    colptrs[6] = 7;
    colptrs[7] = 7;
}