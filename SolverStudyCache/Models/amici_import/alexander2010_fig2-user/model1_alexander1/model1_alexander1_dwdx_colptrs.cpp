#include "sundials/sundials_types.h"

void dwdx_colptrs_model1_alexander1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 6;
    colptrs[2] = 6;
    colptrs[3] = 9;
    colptrs[4] = 12;
    colptrs[5] = 14;
}