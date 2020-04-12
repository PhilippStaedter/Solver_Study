#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_esteban2013(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 5;
    colptrs[3] = 8;
    colptrs[4] = 10;
    colptrs[5] = 15;
}