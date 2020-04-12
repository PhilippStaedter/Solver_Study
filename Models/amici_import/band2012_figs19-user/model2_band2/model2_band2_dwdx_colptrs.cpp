#include "sundials/sundials_types.h"

void dwdx_colptrs_model2_band2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
}