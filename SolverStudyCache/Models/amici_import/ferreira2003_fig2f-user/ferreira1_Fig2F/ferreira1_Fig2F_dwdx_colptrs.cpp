#include "sundials/sundials_types.h"

void dwdx_colptrs_ferreira1_Fig2F(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 2;
    colptrs[3] = 4;
    colptrs[4] = 6;
    colptrs[5] = 8;
    colptrs[6] = 14;
}