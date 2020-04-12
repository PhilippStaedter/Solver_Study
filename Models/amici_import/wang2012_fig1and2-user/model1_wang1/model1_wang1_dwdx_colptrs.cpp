#include "sundials/sundials_types.h"

void dwdx_colptrs_model1_wang1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 5;
    colptrs[2] = 8;
    colptrs[3] = 11;
}