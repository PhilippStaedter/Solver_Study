#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_kouril7(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 3;
    colptrs[3] = 5;
    colptrs[4] = 5;
    colptrs[5] = 6;
    colptrs[6] = 8;
    colptrs[7] = 9;
    colptrs[8] = 10;
    colptrs[9] = 11;
    colptrs[10] = 12;
    colptrs[11] = 12;
}