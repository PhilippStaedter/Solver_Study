#include "sundials/sundials_types.h"

void dwdx_colptrs_kolodkin1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 3;
    colptrs[3] = 4;
}