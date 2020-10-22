#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_bachar1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 0;
    colptrs[2] = 4;
    colptrs[3] = 9;
    colptrs[4] = 12;
}