#include "sundials/sundials_types.h"

void dwdx_colptrs_model1_ribeiro1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 3;
}