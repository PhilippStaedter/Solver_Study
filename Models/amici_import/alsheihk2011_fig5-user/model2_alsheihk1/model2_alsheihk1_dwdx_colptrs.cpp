#include "sundials/sundials_types.h"

void dwdx_colptrs_model2_alsheihk1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 2;
    colptrs[3] = 6;
    colptrs[4] = 9;
    colptrs[5] = 11;
}