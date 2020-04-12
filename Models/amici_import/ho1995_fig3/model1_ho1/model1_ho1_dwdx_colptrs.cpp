#include "sundials/sundials_types.h"

void dwdx_colptrs_model1_ho1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
}