#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_lee3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 2;
    colptrs[3] = 3;
    colptrs[4] = 4;
}