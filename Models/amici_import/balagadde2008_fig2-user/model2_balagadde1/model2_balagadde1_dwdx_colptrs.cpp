#include "sundials/sundials_types.h"

void dwdx_colptrs_model2_balagadde1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 4;
    colptrs[3] = 8;
    colptrs[4] = 12;
    colptrs[5] = 14;
    colptrs[6] = 14;
    colptrs[7] = 14;
}