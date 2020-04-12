#include "sundials/sundials_types.h"

void dwdx_colptrs_lou1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 6;
    colptrs[3] = 10;
    colptrs[4] = 14;
    colptrs[5] = 18;
    colptrs[6] = 22;
}