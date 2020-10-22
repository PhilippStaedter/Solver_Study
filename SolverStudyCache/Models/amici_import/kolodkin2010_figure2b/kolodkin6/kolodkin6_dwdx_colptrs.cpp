#include "sundials/sundials_types.h"

void dwdx_colptrs_kolodkin6(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 8;
    colptrs[5] = 11;
    colptrs[6] = 13;
    colptrs[7] = 14;
    colptrs[8] = 15;
}