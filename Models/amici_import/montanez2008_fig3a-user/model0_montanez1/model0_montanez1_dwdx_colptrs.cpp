#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_montanez1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 6;
    colptrs[3] = 10;
}