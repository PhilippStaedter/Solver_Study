#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_kirschner(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 7;
}