#include "sundials/sundials_types.h"

void dwdx_colptrs_model1_kouril4(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 3;
    colptrs[3] = 4;
    colptrs[4] = 5;
    colptrs[5] = 6;
    colptrs[6] = 6;
}